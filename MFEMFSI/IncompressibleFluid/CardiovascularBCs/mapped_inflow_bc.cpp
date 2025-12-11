#include "mapped_inflow_bc.hpp"

using namespace mfem;
using namespace std;
using namespace tk;
using namespace kdtree_v1;

#ifdef MFEM_USE_MPI

ApproxInflowRate::ApproxInflowRate(string &filename_, int interp_, int N_) : filename(filename_), interp(interp_), N(N_)
{ 
   if (interp == interp_type::FOURIER) MFEM_ASSERT(N > 0, "Error: Specify number of fourier modes > 0!");

   // read in the data from the file
   file_reader.open(filename);
   if(file_reader.fail())
   {
      MFEM_ABORT("Error: flow rate file not found!");
   }

   // NOTE: this will only read properly if file has no extra spaces
   string data;
   do
   {
      // read first col
      double t_curr;
      file_reader >> t_curr;
      if (file_reader.eof()) break;
      t_array.emplace_back(t_curr);

      double q_curr;
      file_reader >> q_curr;
      q_array.emplace_back(q_curr);

   } while (getline(file_reader,data));
   // } while (!file_reader.eof());
   file_reader.close();

   // while (true)
   // {
   //    // read first col
   //    double t_curr;
   //    file_reader >> t_curr;
   //    if (file_reader.eof()) break;
   //    t_array.emplace_back(t_curr);

   //    double q_curr;
   //    file_reader >> q_curr;
   //    q_array.emplace_back(q_curr);
   // }
   // file_reader.close();

   // build the approx
   BuildFlowRateApprox();
}

void ApproxInflowRate::BuildFlowRateApprox()
{
   if (interp == interp_type::FOURIER)
   {
      // least square approx
      MFEM_ABORT("TO DO");
      // Need Eigen library to compute psuedoinverse the rectangular matrix to find coeffs!
   }
   else if (interp == interp_type::SPLINE)
   {
      // setup to match defaults
      s.set_boundary(spline::bd_type::second_deriv,0.0,spline::bd_type::second_deriv,0.0);
      s.set_points(t_array,q_array);
   }
   else
   {
      MFEM_ABORT("Error: Invalid flow rate approx type");
   }
}

double ApproxInflowRate::flow_rate_func(double t)
{
   // if periodic map t to initial period defined in input file
   double t_val = t;
   if (periodic)
   {
      double num_periods = floor(t / period);
      t_val = t - num_periods*period;
   }

   if (interp == interp_type::FOURIER)
   {
      MFEM_ABORT("TO DO");
      return 0.0;
   }
   else if (interp == interp_type::SPLINE)
   {
      return s(t_val);
      // return (*s)(t_val);
   }
   else
   {
      MFEM_ABORT("Error: Invalid flow rate approx type");
      return 0.0;
   }
};

MapToCircularCS::MapToCircularCS(int dim_, ParMesh *pmesh_, Mesh *mesh_, Array<int> inflow_attr_) :
dim(dim_), pmesh(pmesh_), mesh(mesh_), inflow_attr(inflow_attr_)
{
   mesh->GetBoundingBox(bbmin,bbmax);
   
   // create submesh for inflow surface
   SubMesh surface_mesh = SubMesh::CreateFromBoundary(*mesh,inflow_attr);

   // // now extract boundary of inflow surface // DOESNT WORK FOR HIGHER ORDER MESHES
   // Array<int> inflow_bdr_attr(surface_mesh.bdr_attributes.Max());
   // inflow_bdr_attr = 1.0;
   // SubMesh inflow_boundary_mesh = SubMesh::CreateFromBoundary(surface_mesh,inflow_bdr_attr);

   // make refined mesh with uniform points to extract boundary_nodes
   // auto boundary_mesh_refined = inflow_boundary_mesh.MakeRefined(inflow_boundary_mesh,2,BasisType::ClosedUniform); // closed uniform works better
   // boundary_mesh_refined.EnsureNodes();
   // boundary_mesh_refined.SetCurvature(2);
   // boundary_nodes = boundary_mesh_refined.GetNodes(); //
   // auto refined_mesh_uniform = surface_mesh.MakeRefined(surface_mesh,2,BasisType::ClosedUniform);


   auto refined_mesh_uniform = surface_mesh.MakeRefined(surface_mesh,1,BasisType::ClosedUniform);
   refined_mesh_uniform.EnsureNodes();
   refined_mesh_uniform.SetCurvature(2);
   auto uniform_nodes = refined_mesh_uniform.GetNodes();
   auto uniform_mesh_space = refined_mesh_uniform.GetNodalFESpace();
   int ordering = uniform_mesh_space->GetOrdering(); 

   // loop over boundary elements of surface submesh and grab nodal coordinates
   std::vector<double> data;
   Array<int> dofs;
   Array<int> vdofs;
   Vector curr_coords(dim);
   curr_coords = 0.0;
   for (int e = 0; e < uniform_mesh_space->GetNBE(); e++)
   {
      // get bdr element vdofs
      auto bdr_elem = refined_mesh_uniform.GetBdrElement(e);
      // auto dof_trans = uniform_mesh_space->GetBdrElementDofs(e,dofs);
      auto dof_trans = uniform_mesh_space->GetBdrElementVDofs(e,vdofs);

      // get coords from nodes gf
      // vdofs.Print();
      uniform_nodes->GetSubVector(vdofs,curr_coords);
      // curr_coords.Print();

      // // store in data (ordering by vdim)
      int np = vdofs.Size() / dim;
      for (int p = 0; p < np; p++)
      {
         for (int d = 0; d < dim; d++)
         {
            data.emplace_back(curr_coords(p + d*dim));
         }
      }

      // std::cout << "----------\n" << std::endl;
      // std::cout << data[3] << std::endl;
      // std::cout << data[4] << std::endl;
      // std::cout << data[5] << std::endl;
      
   }
   // copy data to coords vector
   {
      coords.SetSize(data.size());
      coords = 0.0;
      for(int i=0; i< data.size(); i++)
      {
         coords(i) = data[i];
      }
   }

   // make refined version with gauss lobatto nodes to compute centroid
   // auto refined_mesh = surface_mesh.MakeRefined(surface_mesh,2,BasisType::GaussLobatto);
   auto refined_mesh = surface_mesh.MakeRefined(surface_mesh,1,BasisType::GaussLobatto);
   refined_mesh.EnsureNodes();
   refined_mesh.SetCurvature(2);

   // compute centroid
   // find center of inlet surface (compute centroid = \int_x x dx / \int_x dx) weighted avg of centroids for each element (weight is area of each elem)
   // loop over elems on proc, compute their contributions to avg and proceed
   double local_area = 0.0;
   double global_area = 0.0;
   Vector center(dim);
   center = 0.0;
   double elem_area = 0.0;
   z_centroid.SetSize(dim);
   z_centroid = 0.0;
   Vector local_sum(dim);
   local_sum = 0.0;
   Vector temp(dim);
   temp = 0.0;
   // for (int e = 0; e < pmesh->GetNE(); e++)
   for (int e = 0; e < refined_mesh.GetNE(); e++)
   {
      // compute element area
      refined_mesh.GetElementCenter(e,center);
      elem_area = refined_mesh.GetElementVolume(e);

      // elem weighted centroid
      temp = center;
      temp *= elem_area;
      local_sum += temp;

      // local area
      local_area += elem_area;
   }

   // communication
   double x_sum = 0.0;
   double y_sum = 0.0;
   double z_sum = 0.0;
   double x_sum_global = 0.0;
   double y_sum_global = 0.0;
   double z_sum_global = 0.0;

   double x_c, y_c, z_c;
   if (dim == 2)
   {
      x_sum = local_sum(0);
      y_sum = local_sum(1);

      // x
      MPI_Allreduce(&x_sum,&x_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      // y
      MPI_Allreduce(&y_sum,&y_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      // area
      MPI_Allreduce(&local_area,&global_area,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      // compute centroid
      x_c = x_sum_global / global_area;
      y_c = y_sum_global / global_area;

      z_centroid(0) = x_c;
      z_centroid(1) = y_c;
   }
   else if (dim == 3)
   {
      x_sum = local_sum(0);
      y_sum = local_sum(1);
      z_sum = local_sum(2);

      // x
      MPI_Allreduce(&x_sum,&x_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      // y
      MPI_Allreduce(&y_sum,&y_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      // z
      MPI_Allreduce(&z_sum,&z_sum_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      // area
      MPI_Allreduce(&local_area,&global_area,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      // compute centroid
      x_c = x_sum_global / global_area;
      y_c = y_sum_global / global_area;
      z_c = z_sum_global / global_area;

      z_centroid(0) = x_c;
      z_centroid(1) = y_c;
      z_centroid(2) = z_c;
   }
   else
   {
      MFEM_ABORT("Error: Computation of centroid in 1D not implemented");
   }

   // if (Mpi::Root())
   // {
   //    std::cout << "x_centroid = " << x_c << std::endl;
   //    std::cout << "y_centroid = " << y_c << std::endl;
   //    std::cout << "z_centroid = " << y_c << std::endl;
   // }

   // compute area using linear form
   // much more accurate area calculation !!!!!
   ConstantCoefficient onecoeff(1.0);
   H1_FECollection fecH1 = H1_FECollection(2, dim);
   FiniteElementSpace pfes(&refined_mesh, &fecH1);
   LinearForm *mass_lf = new LinearForm(&pfes);
   mass_lf->AddDomainIntegrator(new DomainLFIntegrator(onecoeff));
   mass_lf->Assemble();

   GridFunction one_gf(&pfes);
   one_gf.ProjectCoefficient(onecoeff);

   area = mass_lf->operator()(one_gf);

   // compute r_avg
   r_avg = sqrt(area / M_PI);

   // build kdtree
   // probably need list of boundary nodes on each processor, build kdtree on each processor? (YES)
   // finder = new KDTreeClosestPointFinder(*boundary_nodes,dim,bbmin,bbmax);
   finder = new KDTreeClosestPointFinder(coords,Ordering::byVDIM,dim,bbmin,bbmax);

   // clear mesh data
   refined_mesh.Clear();
   refined_mesh_uniform.Clear();
}

double MapToCircularCS::MapToCircle(Vector &x_qp)
{
   // find nearest point on surface and distance
   Vector cp;
   double dist = 0.0;
   int idx;
   finder->FindCPandDist(x_qp,cp,dist,idx);

   //compute r(z)
   double dist_centroid = x_qp.DistanceTo(z_centroid); // distance from centroid
   double r_val = ( dist_centroid / (dist_centroid + dist) ) * r_avg;
   return r_val;
}

MappedInflowCoefficient::MappedInflowCoefficient(int dim_, std::function<double (double)> flow_rate_, int inflow_type_, int order_, ParMesh *pmesh_, Mesh *mesh_, Array<int> inflow_attr_, bool time_depend_) :
VectorCoefficient(dim_), order(order_), time_depend(time_depend_), pmesh(pmesh_), dim(dim_), inflow_attr(inflow_attr_), inflow_type(inflow_type_)
{
   flow_rate = std::move(flow_rate_);
   // radial_inflow = std::move(radial_inflow_);

   map_cs = new MapToCircularCS(dim_,pmesh_,mesh_,inflow_attr);
   r_avg = map_cs->Get_r_avg();
   area = map_cs->GetArea();

   // set time and compute initial normalization integral
   (this)->SetTime(0.0);
}

MappedInflowCoefficient::MappedInflowCoefficient(int dim_, bool from_file, int inflow_type_, int order_, ParMesh *pmesh_, Mesh *mesh_, Array<int> inflow_attr_, bool time_depend_) :
VectorCoefficient(dim_), order(order_), time_depend(time_depend_), pmesh(pmesh_), dim(dim_), inflow_attr(inflow_attr_), inflow_type(inflow_type_), flow_rate_from_file(from_file)
{
   map_cs = new MapToCircularCS(dim_,pmesh_,mesh_,inflow_attr);
   r_avg = map_cs->Get_r_avg();
   area = map_cs->GetArea();

   // set time and compute initial normalization integral
   (this)->SetTime(0.0);
}

void MappedInflowCoefficient::FlowRateFileOptions(std::string filename_, int approx_type, bool periodic_data, int num_fourier_modes)
{
   q_func = new ApproxInflowRate(filename_,approx_type,num_fourier_modes);
   if (periodic_data) q_func->PeriodicFlowRate();
}

// evaluate profile at quad pt in parent element
void MappedInflowCoefficient::Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
{
   // transform integration point to physical space
   Vector x_qp;
   T.Transform(ip,x_qp);

   // evaluate inflow profile
   double inflow_val = cartesian_inflow(x_qp,time);

   // adjust for flow rate
   double flow_rate_val = 0.0;
   if (flow_rate_from_file)
   {
      flow_rate_val = q_func->flow_rate_func(time);
   }
   else
   {
      flow_rate_val = flow_rate(time);
   }

   double normal_velocity = flow_rate_val * inflow_val / integ_inflow;

   // convert velocity in normal, tangential components to cartesian components
   V.SetSize(dim);
   V = 0.0;

   // compute normal in cartesian basis
   Vector nor(dim);
   nor = 0.0;

   MFEM_ASSERT(dim = 3, "Error: Must be 3D mesh for this infow BC");
   CalcOrtho(T.Jacobian(), nor); // normal is in cartesian components of solver

   // normalize nor, see example 18 in general it is not a unit normal
   double normag = 0;
   for (int ii = 0; ii < nor.Size(); ii++)
   {
      normag += nor(ii)*nor(ii);
   }
   normag = sqrt(normag);
   nor *= 1/normag;

   // set output
   // V = nor;
   // V *= -1.0; // since inflow normal points opposite of flow direction
   // V *= normal_velocity;
   V.Set(-normal_velocity,nor);
}

// computes \int_Gamma_in uT(z,t)
void MappedInflowCoefficient::IntegrateInflow(double t)
{
   
   double local_integ = 0.0;
   bool on_surf;

   // loop over boundary elements
   for (int e = 0; e < pmesh -> GetNBE(); e++)
   {
      on_surf = false;

      // check boundary attribute is on surface of interest
      int bdr_attr = pmesh->GetBdrAttribute(e);

      // Note: bdr attributes have same number as specified in GMSH, i.e bdr 1 = bdr attribute 1 in MFEM
      // just when you set the attribute arrays is the bdr attribute need -1;
      for (int ii = 0; ii < inflow_attr.Size(); ii++)
      {
         if(inflow_attr[ii] == 1 && bdr_attr == ii+1) // +1 bc of above comment
         {
            on_surf = true;
         }
      }
      
      // only compute if on surface of interest
      if(on_surf) 
      {
         auto tr = pmesh->GetBdrElementTransformation(e);
         Vector z_qp(pmesh->Dimension());

         // integration rule
         const IntegrationRule *ir = &IntRules.Get(tr->GetGeometryType(), 2*order + 1); // simple integration rule
         for (int j = 0; j < ir->GetNPoints(); j++)
         {
            // setup current integration point
            const IntegrationPoint &ip = ir->IntPoint(j);
            tr->SetIntPoint(&ip);

            // evaluate inflow function at z_qp
            z_qp = 0.0;
            tr->Transform(ip,z_qp);

            // compute integral of inflow velocity profile (potentially time dependent)
            local_integ += ip.weight * cartesian_inflow(z_qp,t) * tr->Weight();
         }
      }
   }

   // communication
   double global_integ = 0.0;
   MPI_Allreduce(&local_integ, &global_integ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   integ_inflow = global_integ;
}

void MappedInflowCoefficient::SetTime(double t)
{
   time = t;

   // if not time dependent then only need to compute \int u^p(r(z)) once
   if (time == 0)
   {
      IntegrateInflow(time);
   }
   else if (time_depend == true)
   {
      // recompute integral normalization
      IntegrateInflow(time);
   }
}

double MappedInflowCoefficient::radial_inflow(double r, double t)
{
   double return_val = 0.0; 
   if (inflow_type == 0) // parabolic
   {
      MFEM_ASSERT(time_depend == false,"Error: for radial parabolic inflow, set time_depend = false");
      return_val =  1 - ((r/r_avg) * (r/r_avg));
   }
   else if (inflow_type == 1) // plug
   {
      return_val = 1.0;
   }
   else
   {
      MFEM_ABORT("Unsupported mapped inflow radial velocity profile");
   } 

   return return_val;
}

// evaluates uT(z,t) = uP(r(z),t)
double MappedInflowCoefficient::cartesian_inflow(Vector &z, double t)
{
   // map z to r
   double r_val = 0.0;
   r_val = map_cs->MapToCircle(z);

   // compute radial inflow
   return radial_inflow(r_val,t);
}

#endif