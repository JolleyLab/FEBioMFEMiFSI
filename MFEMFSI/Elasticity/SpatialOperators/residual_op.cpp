#include "residual_op.hpp"
#include "../../Software/mfem-master/general/forall.hpp"

using namespace dyn_elast;
using namespace mfem;
using namespace std;
using namespace mat;

void ResidualOpElasticity::Mult(const Vector &k, Vector &y) const
{
   MFEM_ABORT("Not implemented");
}

void ResidualOpElasticity::Mult_PC(const Vector &d, const Vector &v, const Vector &a, const Vector &ft, const Vector &fb, Vector &res) const
{
   ////////////////////////////////
   /////// Compute residual /////// 
   ////////////////////////////////

   if (robin_bc_flag)
   {
      curr_v_gf->SetFromTrueDofs(v);
   }

   // internal work
   Stress_form->Mult(d,res);

   // inertial work
   if (!static_problem) Mmat->AddMult(a,res);

   // external work
   res.Add(-1.0,fb);
   res.Add(-1.0,ft);

   // zero out essential dofs in residual vector
   const int N = ess_tdofs.Size();
   const auto idx = ess_tdofs.Read();
   auto Y_RW = res.ReadWrite();
   mfem::forall(N, [=] MFEM_HOST_DEVICE (int i) { Y_RW[idx[i]] = 0.0; });

   // std::cout << setprecision(16);
   // res.Print(cout);
}

Operator &ResidualOpElasticity::GetGradient(const Vector &k) const
{  
   MFEM_ABORT("not implemented");
   return *Jacobian;
}

Operator &ResidualOpElasticity::GetGradient_PC(const Vector &d, const Vector &v, const Vector &a, const double dt, const double am, const double af) const
{
   // return eliminated jacobian evaluated at current state of the body
   // Tangent Stiffness //

   if (static_problem)
   {
      Jacobian = dynamic_cast<HypreParMatrix *>(&Stress_form->GetGradient(d)); // jacobian not owned here
      Jacobian->EliminateBC(ess_tdofs,Operator::DiagonalPolicy::DIAG_ONE);

      return *Jacobian;
   }
   else
   {
      // inertial and internal tangent stiffness
      HypreParMatrix *stress_tmat = dynamic_cast<HypreParMatrix *>(&Stress_form->GetGradient(d));
      
      if (Jacobian) delete Jacobian;
      Jacobian = Add(am,*Mmat,af,*stress_tmat);

      // external work tangent stiffness (TO DO)
      // NOTE: follower loads are included in stress_tmat above!

      // Eliminate essential dofs
      Jacobian->EliminateBC(ess_tdofs,Operator::DiagonalPolicy::DIAG_ONE);


      // Jacobian->PrintMatlab(std::cout);

      // string file_prefix = "my_mat";
      // dynamic_cast<HypreParMatrix *>(Jacobian)->Print(file_prefix.c_str());

      return *Jacobian;  
   }   
}

ResidualOpElasticity::~ResidualOpElasticity()
{
   if(Jacobian && !static_problem) delete Jacobian;
}