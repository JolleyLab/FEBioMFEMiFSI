#pragma once

#include "mfem.hpp"

using VecFuncT = void(const mfem::Vector &x, double t, mfem::Vector &u);
using ScalarFuncT = double(const mfem::Vector &x, double t);
using MatFuncT = void(const mfem::Vector &x, double t, mfem::DenseMatrix &gradu);

/// Container for a Dirichlet boundary condition of the velocity field.
class VelDirichletBC_T
{
public:
   VelDirichletBC_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff), coeff_array(nullptr)
   {
      all_comp = true;
      component = -1;
   }

   VelDirichletBC_T(mfem::Array<int> attr, mfem::Coefficient *coeff_[], int component_)
      : attr(attr), coeff_array(coeff_), coeff(nullptr)
   {
      all_comp = false;
      component = component_;
   }

   VelDirichletBC_T(VelDirichletBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      this->coeff_array = obj.coeff_array;
      this->all_comp = obj.all_comp;
      this->component = obj.component;
      obj.coeff = nullptr;
      obj.coeff_array = nullptr;

      this->zero_bdr = obj.zero_bdr;
      this->zero_bdr_tdofs = obj.zero_bdr_tdofs;
   }

   ~VelDirichletBC_T() 
   { 
      delete coeff; 
      if (coeff_array)
      {
         for (int i=0; coeff_array[i] != nullptr; i++)
         {
            delete coeff_array[i];
         }
         delete[] coeff_array; 
      }
   }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff;
   mfem::Coefficient **coeff_array;
   bool all_comp;
   int component = -1;

   bool zero_bdr = false; // if true, the boundary nodes of this bc are set to zero velocity
   mfem::Array<int> zero_bdr_tdofs;
};

/// Container for a Neumann boundary condition of the velocity field.
class VelNeumannBC_T
{
public:
   VelNeumannBC_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff), mcoeff(NULL), scalar_coeff(NULL)
   {}

   VelNeumannBC_T(mfem::Array<int> attr, mfem::MatrixCoefficient *mcoeff)
      : attr(attr), mcoeff(mcoeff), coeff(NULL), scalar_coeff(NULL)
   {}

   VelNeumannBC_T(mfem::Array<int> attr, mfem::Coefficient *scalar_coeff, bool pressure = false)
      : attr(attr), mcoeff(NULL), coeff(NULL), scalar_coeff(scalar_coeff), pressure(pressure)
   {}

   VelNeumannBC_T(VelNeumannBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
      this->mcoeff = obj.mcoeff;
      obj.mcoeff = nullptr;
      this->scalar_coeff = obj.scalar_coeff;
      obj.scalar_coeff = nullptr;
      this->pressure = obj.pressure;
   }

   ~VelNeumannBC_T() 
   { 
      delete coeff;
      delete mcoeff;
      delete scalar_coeff;
   }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff = nullptr;
   mfem::MatrixCoefficient *mcoeff = nullptr;
   mfem::Coefficient *scalar_coeff = nullptr;
   bool pressure = false; // if true, scalar_coeff represents a pressure value
};

/// Container for a Dirichlet boundary condition of the pressure field.
class PresDirichletBC_T
{
public:
   PresDirichletBC_T(mfem::Array<int> attr, mfem::Coefficient *coeff)
      : attr(attr), coeff(coeff)
   {}

   PresDirichletBC_T(PresDirichletBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
   }

   ~PresDirichletBC_T() { delete coeff; }

   mfem::Array<int> attr;
   mfem::Coefficient *coeff;
};

/// Container for a prescribed boundary displacement bc
class PrescBoundaryDispBC_T
{
public:
   PrescBoundaryDispBC_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff)
   {}

   PrescBoundaryDispBC_T(PrescBoundaryDispBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
   }

   ~PrescBoundaryDispBC_T() { delete coeff; }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff;
};

/// Container for an acceleration term.
class AccelTerm_T
{
public:
   AccelTerm_T(mfem::Array<int> attr, mfem::VectorCoefficient *coeff)
      : attr(attr), coeff(coeff)
   {}

   AccelTerm_T(AccelTerm_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // Move the coefficient pointer
      this->coeff = obj.coeff;
      obj.coeff = nullptr;
   }

   ~AccelTerm_T() { delete coeff; }

   mfem::Array<int> attr;
   mfem::VectorCoefficient *coeff;
};

/// Container for resistance bcs
class ResistanceBC_T
{
public:
   #ifdef MFEM_USE_MPI
   ResistanceBC_T(mfem::Array<int> attr, double R_, double p0_, mfem::ParFiniteElementSpace *vfes_)
      : attr(attr), R(R_), p0(p0_), vfes(vfes_)
   {
      flow_rate = 0.0;
      one.constant = 1.0;
      diag_factor = 0.0;

      AssembleVectors();
   }
   #else
   ResistanceBC_T(mfem::Array<int> attr, double R_, double p0_, mfem::FiniteElementSpace *vfes_)
      : attr(attr), R(R_), p0(p0_), vfes(vfes_)
   {
      flow_rate = 0.0;
      one.constant = 1.0;
      diag_factor = 0.0;

      AssembleVectors();
   }
   #endif

   ResistanceBC_T(ResistanceBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // copy over additional data
      this->R = obj.R;
      this->p0 = obj.p0;
      this->flow_rate = obj.flow_rate;
      this->s_form = obj.s_form;
      obj.s_form = nullptr;
      this->s_vec = obj.s_vec;
      this->bdr_true_dofs = obj.bdr_true_dofs;

      this->diag_factor = obj.diag_factor;

      this->vfes = obj.vfes;
      this->temp = obj.temp;
   }

   void AssembleVectors()
   {  
      if (s_form == nullptr)
      {
         #ifdef MFEM_USE_MPI
         s_form = new mfem::ParLinearForm(vfes);
         #else
         s_form = new mfem::LinearForm(vfes);
         #endif
         s_form->AddBoundaryIntegrator(new mfem::VectorBoundaryFluxLFIntegrator(one),attr);
         
         vfes->GetEssentialTrueDofs(attr,bdr_true_dofs);
      }
      s_form->Assemble();
      #ifdef MFEM_USE_MPI
      temp.SetSize(vfes->GetTrueVSize());
      s_form->ParallelAssemble(temp);
      #endif
      
      // now store only boundary dofs
      s_vec.SetSize(bdr_true_dofs.Size());
      for (int i=0; i<bdr_true_dofs.Size(); i++)
      {
         #ifdef MFEM_USE_MPI
         s_vec[i] = temp[bdr_true_dofs[i]];
         #else
         s_vec[i] = (*s_form)[bdr_true_dofs[i]];
         #endif
      }
   };

   mfem::Vector *Get_S_vec()
   {
      return &s_vec;
   }

   mfem::Array<int> attr;
   double R, p0;
   double flow_rate;
   mfem::ConstantCoefficient one;
   mfem::Vector s_vec;
   mfem::Array<int> bdr_true_dofs;
   #ifdef MFEM_USE_MPI
   mfem::ParLinearForm *s_form = nullptr;
   mfem::ParFiniteElementSpace *vfes = nullptr; // not owned
   #else
   mfem::LinearForm *s_form = nullptr;
   mfem::FiniteElementSpace *vfes = nullptr; // not owned
   #endif
   mfem::Vector temp;

   double diag_factor;

   ~ResistanceBC_T()
   {
      delete s_form;
   }
};

/// Container for windkessel bcs
class WindkesselBC_T
{
public:
   #ifdef MFEM_USE_MPI
   WindkesselBC_T(mfem::Array<int> attr, double R_, double Rd_, double C_, mfem::Coefficient *p_0_, std::function<double(double)> p_d_, double init_flow_rate_, 
                  mfem::ParFiniteElementSpace *vfes_)
      : attr(attr), R(R_), Rd(Rd_), C(C_), p_0(p_0_), p_d(std::move(p_d_)), vfes(vfes_)
   {
      init_flow_rate = init_flow_rate_;
      q_np1 = init_flow_rate; // initialize flow rate

      // initialize variables
      q_nm1 = 0.0;
      q_n = 0.0;
      q_npalphaf = 0.0;
      t_npalphaf = 0.0;
      t_np1 = 0.0;
      t_n = 0.0;
      t_nm1 = 0.0;
      dt = 0.0;
      integ_hist_terms = 0.0;
      integ_curr_terms = 0.0;
      hist_terms = 0.0;

      diag_factor = 0.0;

      one.constant = 1.0;
      AssembleSVEC();
      AssemblePVEC();

      // check if p_0 is nonzero
      double p_vec_norm = sqrt(mfem::InnerProduct(vfes->GetComm(),p_vec,p_vec));
      if (p_vec_norm <= 1.0e-12)
      {
         p_flag = false;
      }
   }
   #else
   WindkesselBC_T(mfem::Array<int> attr, double R_, double Rd_, double C_, mfem::Coefficient *p_0_, std::function<double(double)> p_d_, double init_flow_rate_, 
                  mfem::FiniteElementSpace *vfes_)
      : attr(attr), R(R_), Rd(Rd_), C(C_), p_0(p_0_), p_d(std::move(p_d_)), vfes(vfes_)
   {
      init_flow_rate = init_flow_rate_;
      q_np1 = init_flow_rate; // initialize flow rate

      // initialize variables
      q_nm1 = 0.0;
      q_n = 0.0;
      q_npalphaf = 0.0;
      t_npalphaf = 0.0;
      t_np1 = 0.0;
      t_n = 0.0;
      t_nm1 = 0.0;
      dt = 0.0;
      integ_hist_terms = 0.0;
      integ_curr_terms = 0.0;
      hist_terms = 0.0;

      diag_factor = 0.0;

      one.constant = 1.0;
      AssembleSVEC();
      AssemblePVEC();

      // check if p_0 is nonzero
      double p_vec_norm = sqrt(mfem::InnerProduct(p_vec,p_vec));
      if (p_vec_norm <= 1.0e-12)
      {
         p_flag = false;
      }
   }
   #endif

   #ifdef MFEM_USE_MPI
   WindkesselBC_T(mfem::Array<int> attr, double R_, double Rd_, double C_, double p_0_, std::function<double(double)> p_d_, double init_flow_rate_, 
                  mfem::ParFiniteElementSpace *vfes_)
      : attr(attr), R(R_), Rd(Rd_), C(C_), p_0(nullptr), p_d(std::move(p_d_)), vfes(vfes_)
   {
      init_flow_rate = init_flow_rate_;
      q_np1 = init_flow_rate; // initialize flow rate

      // initialize variables
      q_nm1 = 0.0;
      q_n = 0.0;
      q_npalphaf = 0.0;
      t_npalphaf = 0.0;
      t_np1 = 0.0;
      t_n = 0.0;
      t_nm1 = 0.0;
      dt = 0.0;
      integ_hist_terms = 0.0;
      integ_curr_terms = 0.0;
      hist_terms = 0.0;

      diag_factor = 0.0;

      one.constant = 1.0;
      p_0 = new mfem::ConstantCoefficient(p_0_);
      own_p_0 = true;
      AssembleSVEC();
      AssemblePVEC();

      // check if p_0 is nonzero
      double p_vec_norm = sqrt(mfem::InnerProduct(vfes->GetComm(),p_vec,p_vec));
      if (p_vec_norm <= 1.0e-12)
      {
         p_flag = false;
      }
   }
   #else
   WindkesselBC_T(mfem::Array<int> attr, double R_, double Rd_, double C_, double p_0_, std::function<double(double)> p_d_, double init_flow_rate_, 
                  mfem::FiniteElementSpace *vfes_)
      : attr(attr), R(R_), Rd(Rd_), C(C_), p_0(nullptr), p_d(std::move(p_d_)), vfes(vfes_)
   {
      init_flow_rate = init_flow_rate_;
      q_np1 = init_flow_rate; // initialize flow rate

      // initialize variables
      q_nm1 = 0.0;
      q_n = 0.0;
      q_npalphaf = 0.0;
      t_npalphaf = 0.0;
      t_np1 = 0.0;
      t_n = 0.0;
      t_nm1 = 0.0;
      dt = 0.0;
      integ_hist_terms = 0.0;
      integ_curr_terms = 0.0;
      hist_terms = 0.0;

      diag_factor = 0.0;

      one.constant = 1.0;
      p_0 = new mfem::ConstantCoefficient(p_0_);
      own_p_0 = true;
      AssembleSVEC();
      AssemblePVEC();

      // check if p_0 is nonzero
      double p_vec_norm = sqrt(mfem::InnerProduct(p_vec,p_vec));
      if (p_vec_norm <= 1.0e-12)
      {
         p_flag = false;
      }
   }
   #endif

   #ifdef MFEM_USE_MPI
   WindkesselBC_T(mfem::Array<int> attr, double R_, double Rd_, double C_, double p_0_, double p_d_, double init_flow_rate_, 
                  mfem::ParFiniteElementSpace *vfes_)
      : attr(attr), R(R_), Rd(Rd_), C(C_), p_0(nullptr), p_d_val(p_d_), vfes(vfes_)
   {
      init_flow_rate = init_flow_rate_;
      q_np1 = init_flow_rate; // initialize flow rate

      // initialize variables
      q_nm1 = 0.0;
      q_n = 0.0;
      q_npalphaf = 0.0;
      t_npalphaf = 0.0;
      t_np1 = 0.0;
      t_n = 0.0;
      t_nm1 = 0.0;
      dt = 0.0;
      integ_hist_terms = 0.0;
      integ_curr_terms = 0.0;
      hist_terms = 0.0;

      diag_factor = 0.0;

      one.constant = 1.0;
      p_0 = new mfem::ConstantCoefficient(p_0_);
      own_p_0 = true;
      AssembleSVEC();
      AssemblePVEC();

      // check if p_0 is nonzero
      double p_vec_norm = sqrt(mfem::InnerProduct(vfes->GetComm(),p_vec,p_vec));
      if (p_vec_norm <= 1.0e-12)
      {
         p_flag = false;
      }
   }
   #else
   WindkesselBC_T(mfem::Array<int> attr, double R_, double Rd_, double C_, double p_0_, double p_d_, double init_flow_rate_, 
                  mfem::FiniteElementSpace *vfes_)
      : attr(attr), R(R_), Rd(Rd_), C(C_), p_0(nullptr), p_d_val(p_d_), vfes(vfes_)
   {
      init_flow_rate = init_flow_rate_;
      q_np1 = init_flow_rate; // initialize flow rate

      // initialize variables
      q_nm1 = 0.0;
      q_n = 0.0;
      q_npalphaf = 0.0;
      t_npalphaf = 0.0;
      t_np1 = 0.0;
      t_n = 0.0;
      t_nm1 = 0.0;
      dt = 0.0;
      integ_hist_terms = 0.0;
      integ_curr_terms = 0.0;
      hist_terms = 0.0;

      diag_factor = 0.0;

      one.constant = 1.0;
      p_0 = new mfem::ConstantCoefficient(p_0_);
      own_p_0 = true;
      AssembleSVEC();
      AssemblePVEC();

      // check if p_0 is nonzero
      double p_vec_norm = sqrt(mfem::InnerProduct(p_vec,p_vec));
      if (p_vec_norm <= 1.0e-12)
      {
         p_flag = false;
      }
   }
   #endif
   
   void AssembleSVEC()
   {
      if (s_form == nullptr)
      {
         #ifdef MFEM_USE_MPI
         s_form = new mfem::ParLinearForm(vfes);
         #else
         s_form = new mfem::LinearForm(vfes);
         #endif
         s_form->AddBoundaryIntegrator(new mfem::VectorBoundaryFluxLFIntegrator(one),attr);

         vfes->GetEssentialTrueDofs(attr,bdr_true_dofs);
      }
      s_form->Assemble();
      #ifdef MFEM_USE_MPI
      temp.SetSize(vfes->GetTrueVSize());      
      s_form->ParallelAssemble(temp);
      #endif

      // now store only boundary dofs
      s_vec.SetSize(bdr_true_dofs.Size());
      for (int i=0; i<bdr_true_dofs.Size(); i++)
      {
         #ifdef MFEM_USE_MPI
         s_vec[i] = temp[bdr_true_dofs[i]];
         #else
         s_vec[i] = (*s_form)[bdr_true_dofs[i]];
         #endif
      }
   };

   void AssemblePVEC()
   {
      if (p_flag == true)
      {
         // p_0 vec
         if (p_form == nullptr)
         {
            #ifdef MFEM_USE_MPI
            p_form = new mfem::ParLinearForm(vfes);
            #else
            p_form = new mfem::LinearForm(vfes);
            #endif
            p_form->AddBoundaryIntegrator(new mfem::VectorBoundaryFluxLFIntegrator(*p_0),attr);
         }         
         p_form->Assemble();
         #ifdef MFEM_USE_MPI
         temp.SetSize(vfes->GetTrueVSize());
         p_form->ParallelAssemble(temp);
         #endif

         // now store only boundary dofs
         p_vec.SetSize(bdr_true_dofs.Size());
         for (int i=0; i<bdr_true_dofs.Size(); i++)
         {
            #ifdef MFEM_USE_MPI
            p_vec[i] = temp[bdr_true_dofs[i]];
            #else
            p_vec[i] = (*p_form)[bdr_true_dofs[i]];
            #endif
         }
      }      
   }

   mfem::Vector *Get_S_vec()
   {
      return &s_vec;
   }

   mfem::Vector *Get_p_vec()
   {
      return &p_vec;
   }

   WindkesselBC_T(WindkesselBC_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;

      // copy over additional data
      this->R = obj.R;
      this->Rd = obj.Rd;
      this->C = obj.C;
      this->p_0 = obj.p_0;
      this->own_p_0 = obj.own_p_0;
      obj.p_0 = nullptr;
      this->p_d = std::move(obj.p_d); // not sure about this one?

      this->init_flow_rate = obj.init_flow_rate;

      this->q_nm1 = obj.q_nm1;
      this->q_n = obj.q_n;
      this->q_np1 = obj.q_np1;
      this->q_npalphaf = obj.q_npalphaf;

      this->t_nm1 = obj.t_nm1;
      this->t_n = obj.t_n;
      this->t_np1 = obj.t_np1;
      this->t_npalphaf = obj.t_npalphaf;

      this->dt = obj.dt;
      this->integ_curr_terms = obj.integ_curr_terms;
      this->integ_hist_terms = obj.integ_hist_terms;
      this->hist_terms = obj.hist_terms;

      this->s_form = obj.s_form;
      obj.s_form = nullptr;
      this->s_vec = obj.s_vec;

      this->p_form = obj.p_form;
      obj.p_form = nullptr;
      this->p_vec = obj.p_vec;
      this->p_flag = obj.p_flag;

      this->diag_factor = obj.diag_factor;

      this->bdr_true_dofs = obj.bdr_true_dofs;

      this->vfes = obj.vfes;
      this->temp = obj.temp;

      this->p_d_val = obj.p_d_val;
   }

   mfem::Array<int> attr;
   double R, Rd, C;
   mfem::Coefficient *p_0;
   bool own_p_0 = false;
   std::function<double(double)> p_d;
   double p_d_val;
   double init_flow_rate;
   double q_nm1, q_n, q_np1, q_npalphaf;
   double t_nm1, t_n, t_np1;
   double t_npalphaf;
   double dt;
   double integ_hist_terms;
   double integ_curr_terms;
   double hist_terms;

   #ifdef MFEM_USE_MPI
   mfem::ParLinearForm *s_form = nullptr;
   mfem::ParLinearForm *p_form = nullptr;
   mfem::ParFiniteElementSpace *vfes = nullptr; // not owned
   #else
   mfem::LinearForm *s_form = nullptr;
   mfem::LinearForm *p_form = nullptr;
   mfem::FiniteElementSpace *vfes = nullptr; // not owned
   #endif

   mfem::ConstantCoefficient one;
   mfem::Vector s_vec, p_vec;
   bool p_flag = true;   
   mfem::Vector temp;
   double diag_factor;
   mfem::Array<int> bdr_true_dofs;

   ~WindkesselBC_T()
   {
      delete s_form;
      delete p_form;
      if (own_p_0) delete p_0;
   }
};

/// Container for backflow stabilization
class BackflowStabilization_T
{
public:
   BackflowStabilization_T(mfem::Array<int> attr, double beta_)
      : attr(attr), beta(beta_)
   {};

   BackflowStabilization_T(BackflowStabilization_T &&obj)
   {
      // Deep copy the attribute array
      this->attr = obj.attr;
      this->beta = obj.beta;
   };

   mfem::Array<int> attr;
   double beta;
};