#include "block_prec.hpp"

#include <iomanip>

using namespace mfem;
using namespace std;

JacobianPreconditioner::JacobianPreconditioner(Array<ParFiniteElementSpace *>
                                               &fes,
                                               Array<int> &offsets)
   : Solver(offsets[2]), block_trueOffsets(offsets)
{
   fes.Copy(spaces);

   gamma = 0.00001;

   // form pressure mass matrix
   PressureMassForm = new ParBilinearForm(spaces[1]);
   PressureMassForm->AddDomainIntegrator(new MassIntegrator());
   PressureMassForm->Assemble();
   PressureMassForm->Finalize();
   pressure_mass = PressureMassForm->ParallelAssemble();

   // The mass matrix and preconditioner do not change every Newton cycle, so
   // we only need to define them once
   HypreBoomerAMG *mass_prec_amg = new HypreBoomerAMG();
   mass_prec_amg->SetPrintLevel(0);

   mass_prec = mass_prec_amg;

   CGSolver *mass_pcg_iter = new CGSolver(spaces[1]->GetComm());
   mass_pcg_iter->SetRelTol(1e-12);
   mass_pcg_iter->SetAbsTol(1e-12);
   mass_pcg_iter->SetMaxIter(200);
   mass_pcg_iter->SetPrintLevel(0);
   mass_pcg_iter->SetPreconditioner(*mass_prec);
   mass_pcg_iter->SetOperator(*pressure_mass);
   mass_pcg_iter->iterative_mode = false;

   mass_solver = mass_pcg_iter;

   // The stiffness matrix does change every Newton cycle, so we will define it
   // during SetOperator
}

void JacobianPreconditioner::Mult(const Vector &k, Vector &y) const
{
   // Extract the blocks from the input and output vectors
   Vector disp_in;
   disp_in.MakeRef(const_cast<Vector&>(k), block_trueOffsets[0],
                   block_trueOffsets[1]-block_trueOffsets[0]);
   Vector pres_in;
   pres_in.MakeRef(const_cast<Vector&>(k), block_trueOffsets[1],
                   block_trueOffsets[2]-block_trueOffsets[1]);

   Vector disp_out;
   disp_out.MakeRef(y, block_trueOffsets[0],
                    block_trueOffsets[1]-block_trueOffsets[0]);
   Vector pres_out;
   pres_out.MakeRef(y, block_trueOffsets[1],
                    block_trueOffsets[2]-block_trueOffsets[1]);

   Vector temp(block_trueOffsets[1]-block_trueOffsets[0]);
   Vector temp2(block_trueOffsets[1]-block_trueOffsets[0]);

   // Perform the block elimination for the preconditioner
   mass_solver->Mult(pres_in, pres_out);
   pres_out *= -gamma;

   jacobian->GetBlock(0,1).Mult(pres_out, temp);
   subtract(disp_in, temp, temp2);

   stiff_solver->Mult(temp2, disp_out);

   disp_out.SyncAliasMemory(y);
   pres_out.SyncAliasMemory(y);
}

void JacobianPreconditioner::SetOperator(const Operator &op)
{
   jacobian = (BlockOperator *) &op;

   // Initialize the stiffness preconditioner and solver
   if (stiff_prec == NULL)
   {
      HypreBoomerAMG *stiff_prec_amg = new HypreBoomerAMG();
      stiff_prec_amg->SetPrintLevel(0);

      stiff_prec = stiff_prec_amg;

      GMRESSolver *stiff_pcg_iter = new GMRESSolver(spaces[0]->GetComm());
      stiff_pcg_iter->SetRelTol(1e-8);
      stiff_pcg_iter->SetAbsTol(1e-8);
      stiff_pcg_iter->SetMaxIter(200);
      stiff_pcg_iter->SetPrintLevel(0);
      stiff_pcg_iter->SetPreconditioner(*stiff_prec);
      stiff_pcg_iter->iterative_mode = false;

      stiff_solver = stiff_pcg_iter;
   }

   // At each Newton cycle, compute the new stiffness AMG preconditioner by
   // updating the iterative solver which, in turn, updates its preconditioner
   stiff_solver->SetOperator(jacobian->GetBlock(0,0));
}

JacobianPreconditioner::~JacobianPreconditioner()
{
   delete mass_solver;
   delete mass_prec;
   delete stiff_prec;
   delete stiff_solver;
   delete PressureMassForm;
   delete pressure_mass;
}