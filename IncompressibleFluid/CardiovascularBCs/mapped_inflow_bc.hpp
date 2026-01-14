#pragma once

#include "mfem.hpp"
#include "../../Utilities/kdtree.hpp"
#include "spline.hpp"

using namespace mfem;
using namespace std;
using namespace tk;
using namespace kdtree_v1;

#ifdef MFEM_USE_MPI

// class to approx inflow rate given at discrete points
// using either a fourier series approx of N modes or cubic spline interp
/*
Input expects: 
   file containing (x,y) pairs
   choice of interp
   if fourier, N

file should be formatted as
t1 q1
t2 q2
.....
// IMPORTANT: NO SPACES AT THE END or BEGINNING! file contains only the rows with the data!
where q1 is the flow rate INTO the domain!

Also, provides option to specify one period, the approx will then be periodic about this period
Assumes initial t1 = 0.
*/
class ApproxInflowRate
{
protected:
   ifstream file_reader;
   string filename;
   int interp;

   // arrays for storing file contents
   std::vector<double> t_array, q_array;

   // periodic?
   bool periodic = false;
   double period; // last time val in file

   // related to fourier series approx
   int N;

   // spline interp
   tk::spline s;

public:
   // constructor
   ApproxInflowRate(string &filename_, int interp_, int N_ = 0);

   void PeriodicFlowRate()
   {
      periodic = true;
      if (!t_array.empty())
      {
         period = t_array.back();
      }
      else
      {
         MFEM_ABORT("Internal error: t_array not set");
      }
   }

   // build approx
   void BuildFlowRateApprox();

   // function to eval flow rate
   double flow_rate_func(double t);

   enum interp_type {FOURIER=0,SPLINE=1};

   ~ApproxInflowRate() 
   {

   };
};

// class to map from non-circular to circular cross section
// evaluates r(z)
class MapToCircularCS
{
protected:
   int dim;
   ParMesh *pmesh = nullptr;
   Mesh *mesh = nullptr;

   // geometry parameters
   Vector bbmin, bbmax;
   double r_avg;
   double area;
   Vector z_centroid; // centroid

   // KDTree for finding nearest boundary point and distance
   KDTreeClosestPointFinder *finder = nullptr;
   GridFunction *boundary_nodes = nullptr;
   Vector coords;

   // inflow bdr attribute
   Array<int> inflow_attr;

public:
   // constructor
   MapToCircularCS(int dim_, ParMesh *pmesh_, Mesh *mesh_, Array<int> inflow_attr_);

   // maps z to r
   double MapToCircle(Vector &x_qp);

   // get geometry parameters
   double Get_r_avg()
   {
      return r_avg;
   };

   double GetArea()
   {
      return area;
   };
};

// Class that computes inflow profile on non-circular cross section
// with desired flow rate and form on circular cross section
class MappedInflowCoefficient : public VectorCoefficient
{
protected:
   // // flow rate in cartesian components
   std::function<double (double)> flow_rate;


   // TO DO: interface for user specified radial inflow function!
   int inflow_type;
   // inflow profile in radial coordinates
   // std::function<double (double, double)> radial_inflow;
   // integral normalization
   double integ_inflow = 0.0;
   // flag for time dependent inflow
   bool time_depend;

   // order of FEM approx
   int order;

   // mapping to circular cross section
   MapToCircularCS *map_cs = nullptr;
   double r_avg;
   double area;

   // inflow boundary attribute
   Array<int> inflow_attr;

   // parmesh
   ParMesh *pmesh = nullptr;
   int dim;

   // options related to using a flow rate from a file
   bool flow_rate_from_file = false;
   ApproxInflowRate *q_func = nullptr;

public:
   // constructor for analytical flow_rate
   MappedInflowCoefficient(int dim_, std::function<double (double)> flow_rate_, int inflow_type_, int order_, ParMesh *pmesh_, Mesh *mesh_, Array<int> inflow_attr_, bool time_depend_ = false);

   // constructor for flow rate from a file
   MappedInflowCoefficient(int dim_, bool from_file, int inflow_type_, int order_, ParMesh *pmesh_, Mesh *mesh_, Array<int> inflow_attr_, bool time_depend_ = false);

   // function to set additional options for flow rate from file
   void FlowRateFileOptions(std::string filename_, int approx_type, bool periodic_data, int num_fourier_modes = 0);

   // evaluate profile at quad pt in parent element
   virtual void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override;

   // computes \int_Gamma_in uT(z,t)
   void IntegrateInflow(double t);

   virtual void SetTime(double t) override;

   // evaluate radial inflow profile
   double radial_inflow(double r, double t);

   // evaluates uT(z,t) = uP(r(z),t)
   double cartesian_inflow(Vector &z, double t);

   enum inflowTYPE {Parabolic = 0, Plug = 1};

   ~MappedInflowCoefficient()
   {
      delete q_func;
   }
};

#endif