#include "ImmersedFSIFindPts.hpp"
#ifdef MFEM_USE_MPI
#include <HYPRE.h>
#include <_hypre_IJ_mv.h>
#endif

// Ignore warnings from the gslib header (GCC version)
#ifdef MFEM_HAVE_GCC_PRAGMA_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

#ifdef MFEM_USE_GSLIB
// External GSLIB header (the MFEM header is gslib.hpp)
namespace gslib
{
#include "gslib.h"
}
#endif

#ifdef MFEM_HAVE_GCC_PRAGMA_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif

using namespace mfem;
using namespace immersed_solid;

#ifdef MFEM_USE_GSLIB
ImmersedFSIFindPts::ImmersedFSIFindPts() : FindPointsGSLIB()
{

}

#ifdef MFEM_USE_MPI
ImmersedFSIFindPts::ImmersedFSIFindPts(MPI_Comm comm_) : mfem::FindPointsGSLIB(comm_)
{

}
#endif

void ImmersedFSIFindPts::SpreadForces(const Vector &point_pos,
                                         const Vector &field_in,
                                         GridFunction &field_out,
                                         int point_pos_ordering)
{
   FindPoints(point_pos, point_pos_ordering);
   Spread(field_in, field_out);

   #ifdef MFEM_USE_MPI
   // communication, sum values across shared dofs and distribute to all shared dofs (shared across proc boundary)
   GroupCommunicator &gcomm = dynamic_cast<ParGridFunction *>(&field_out)->ParFESpace()->GroupComm();
   Array<double> dest_array(field_out.GetData(), field_out.Size());
   gcomm.Reduce<double>(dest_array, GroupCommunicator::Sum);
   gcomm.Bcast(dest_array);
   #endif
}

void ImmersedFSIFindPts::Spread(const Vector &field_in, GridFunction &field_out)
{
   int ncomp   = field_out.VectorDim(),
       nptorig = points_cnt,
       npt     = points_cnt;

   MFEM_ASSERT(field_out.Size() != 0, "field_out has not been sized!");
   field_out = default_interp_value;

   dim = field_out.FESpace()->GetMesh()->Dimension();

   if (gsl_comm->np == 1) // serial
   {
      int curr_elem;
      DofTransformation *doftrans;
      Vector shape;      
      Array<int> vdofs, vdofs_solid;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] == 2) { continue; }
         // set int point for curr solid node in fluid mesh
         IntegrationPoint ip;
         ip.Set2(gsl_mfem_ref.GetData()+index*dim);
         if (dim == 3) { ip.z = gsl_mfem_ref(index*dim + 2); }

         // get fluid element the solid point is in
         curr_elem = (int) gsl_mfem_elem[index];
         doftrans = field_out.FESpace()->GetElementVDofs(curr_elem,vdofs); // should fluid fes!
         const FiniteElement &fe = *field_out.FESpace()->GetFE(curr_elem);
         int nd = fe.GetDof();
         shape.SetSize(nd);

         // evaluate shape funcs at current solid node
         fe.CalcShape(ip,shape);        

         // multiply by nodal solid force and add to output
         for (int a = 0; a < nd; a++)
         {
            // nodal solid force * shape(eval at solid node)
            for (int i = 0; i < ncomp; i++)
            {
               int id = field_out.FESpace()->GetOrdering() == Ordering::byNODES ? index + i*npt : index*ncomp + i;
               field_out[vdofs[a + i*nd]] += shape(a) * field_in(id);
            }
         }         
      }
   }
   else // parallel
   {
      // Determine number of points to be sent
      int nptsend = 0;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] != 2) { nptsend +=1; }
      }

      // Pack data to send via crystal router
      // Add in here solid nodal forces for each solid nodal point
      struct gslib::array *outpt = new gslib::array;
      struct out_pt { double r[3], f[3]; uint index, el, proc; };
      struct out_pt *pt;
      array_init(struct out_pt, outpt, nptsend);
      outpt->n=nptsend;
      pt = (struct out_pt *)outpt->ptr;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] == 2) { continue; }
         pt->index = index;
         pt->proc  = gsl_proc[index];
         pt->el    = gsl_mfem_elem[index];

         for (int d = 0; d < dim; ++d) 
         { 
            pt->r[d]= gsl_mfem_ref(index*dim + d);
         }

         // nodal force in solid
         for (int c = 0; c < ncomp; ++c)
         {            
            // assumes mesh nodes and force are ordered the same way!
            // assumes field in and out are ordered in same way!
            int idx = field_out.FESpace()->GetOrdering() == Ordering::byNODES ? index + c*npt : index*ncomp + c;
            pt->f[c]= field_in(idx);
         }
         ++pt;
      }

      // Transfer data to target MPI ranks
      sarray_transfer(struct out_pt, outpt, proc, 1, cr);

      // Spread the forces
      npt = outpt->n;
      pt = (struct out_pt *)outpt->ptr;

      int curr_elem;
      DofTransformation *doftrans;
      ElementTransformation *eltrans;
      Vector shape, elvect;
      Array<int> vdofs, vdofs_solid;

      // loop over solid nodes sent to curr proc
      for (int index = 0; index < npt; index++)
      {
         // Identify fluid element for solid node A
         curr_elem = (int) pt->el;
         doftrans = field_out.FESpace()->GetElementVDofs(curr_elem,vdofs);
         eltrans = field_out.FESpace()->GetElementTransformation(curr_elem);
         const FiniteElement &fe = *field_out.FESpace()->GetFE(curr_elem);
         shape.SetSize(fe.GetDof());
         elvect.SetSize(fe.GetDof() * ncomp);
         elvect = 0.0;
         
         // evaluate shape funcs at current solid node
         IntegrationPoint ip;
         ip.Set3(&pt->r[0]);
         fe.CalcShape(ip,shape);

         // multiply by nodal solid force and add to output
         int nd = fe.GetDof();
         for (int a = 0; a < fe.GetDof(); a++)
         {
            // nodal solid force * shape(eval at solid node)
            field_out[vdofs[a]] += shape(a) * pt->f[0];
            if (ncomp > 1) field_out[vdofs[a + nd]] += shape(a) * pt->f[1];
            if (dim == 3 && ncomp > 1)
            {
               field_out[vdofs[a + 2*nd]] += shape(a) * pt->f[2];
            }
         }
         ++pt;
      }

      // free memory
      array_free(outpt);
      delete outpt;

   } // parallel
}

void ImmersedFSIFindPts::TransferIndicator(const Vector &point_pos,
                                         GridFunction &field_out,
                                         int point_pos_ordering)
{
   FindPoints(point_pos, point_pos_ordering);
   Transfer(field_out);

   #ifdef MFEM_USE_MPI
   // communication, get largest value across shared dofs and distribute to all shared dofs (shared across proc boundary)
   GroupCommunicator &gcomm = dynamic_cast<ParGridFunction *>(&field_out)->ParFESpace()->GroupComm();
   Array<double> dest_array(field_out.GetData(), field_out.Size());
   gcomm.Reduce<double>(dest_array, GroupCommunicator::Max); // either 1 or 0
   gcomm.Bcast(dest_array);
   #endif
}

void ImmersedFSIFindPts::Transfer(GridFunction &field_out)
{
   int ncomp   = field_out.VectorDim(),
       nptorig = points_cnt,
       npt     = points_cnt;

   MFEM_ASSERT(field_out.Size() != 0, "field_out has not been sized!");
   field_out = default_interp_value;

   dim = field_out.FESpace()->GetMesh()->Dimension();

   if (gsl_comm->np == 1) // serial
   {
      int curr_elem;   
      Array<int> vdofs;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] == 2) { continue; }

         // get fluid element the solid point is in
         curr_elem = (int) gsl_mfem_elem[index];
         field_out.FESpace()->GetElementVDofs(curr_elem,vdofs);   
         
         // set fluid dofs on the element to 1
         field_out.SetSubVector(vdofs,1.0);  
      }
   }
   else // parallel
   {
      // Determine number of points to be sent
      int nptsend = 0;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] != 2) { nptsend +=1; }
      }

      // Pack data to send via crystal router
      // Add in here solid nodal forces for each solid nodal point
      struct gslib::array *outpt = new gslib::array;
      struct out_pt { uint index, el, proc; };
      struct out_pt *pt;
      array_init(struct out_pt, outpt, nptsend);
      outpt->n=nptsend;
      pt = (struct out_pt *)outpt->ptr;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] == 2) { continue; }
         pt->index = index;
         pt->proc  = gsl_proc[index];
         pt->el    = gsl_mfem_elem[index];

         ++pt;
      }

      // Transfer data to target MPI ranks
      sarray_transfer(struct out_pt, outpt, proc, 1, cr);

      npt = outpt->n;
      pt = (struct out_pt *)outpt->ptr;

      int curr_elem;
      Array<int> vdofs;
      // loop over solid nodes sent to curr proc
      for (int index = 0; index < npt; index++)
      {
         // Identify fluid element for solid node A
         curr_elem = (int) pt->el;
         field_out.FESpace()->GetElementVDofs(curr_elem,vdofs);

         // set fluid dofs on the element to 1
         field_out.SetSubVector(vdofs,1.0);

         ++pt;
      }

      // free memory
      array_free(outpt);
      delete outpt;

   } // parallel
}

void ImmersedFSIFindPts::AssembleCouplingMatrix(const Vector &solid_nodes, FiniteElementSpace *fes_in, FiniteElementSpace *fes_out, int orderingtype_)
{
   // find solid nodes on fluid grid
   FindPoints(solid_nodes, orderingtype_);

   #ifdef MFEM_USE_MPI
   // C = fluid dofs x solid dofs
   // fes_in = fluid fes, fes_out = solid fes
   if (C) delete C;
   
   // Setup HypreIJ matrix
   HYPRE_BigInt *row_starts = dynamic_cast<ParFiniteElementSpace*>(fes_in)->GetTrueDofOffsets();
   int row_start = row_starts[0];
   int row_end = row_starts[1];

   HYPRE_BigInt *col_starts = dynamic_cast<ParFiniteElementSpace*>(fes_out)->GetTrueDofOffsets();
   int col_start = col_starts[0];
   int col_end = col_starts[1];

   HYPRE_IJMatrixCreate(dynamic_cast<ParFiniteElementSpace*>(fes_in)->GetComm(),row_start,row_end-1,col_start,col_end-1,&ij_matrix);
   HYPRE_IJMatrixSetObjectType(ij_matrix,HYPRE_PARCSR);
   
   // // TODO: Specify max rowsize or diag/offd exact size for performance
   // HYPRE_IJMatrixSetRowSizes(ij_matrix, max_rowsize);
   
   HYPRE_IJMatrixInitialize(ij_matrix);

   AssembleC(fes_in,fes_out,solid_nodes,orderingtype_);

   HYPRE_Int ret = HYPRE_IJMatrixAssemble(ij_matrix);
   ret = HYPRE_IJMatrixGetObject(ij_matrix, (void**)&pcsr_matrix);
   C = new HypreParMatrix((hypre_ParCSRMatrix*)(pcsr_matrix), true);

   if (CT) delete CT;
   CT = C->Transpose();

   #else
   // C = fluid dofs x solid dofs
   // reform sparsity pattern
   if (C) delete C;
   C = new SparseMatrix(fes_in->GetVSize(),fes_out->GetVSize());
   *C = 0.0;

   AssembleC(fes_in,fes_out);

   // Finalize
   C->Finalize();

   // Transpose the coupling matrix
   if (CT) delete CT;
   CT = Transpose(*C);
   #endif
}

#ifndef MFEM_USE_MPI
void ImmersedFSIFindPts::AssembleC(FiniteElementSpace *fes_in, FiniteElementSpace *fes_out)
{
   int ncomp   = fes_out->GetVDim(),
       nptorig = points_cnt,
       npt     = points_cnt;

   const int skip_zeros = 0;

   dim = fes_out->GetMesh()->Dimension();

   MFEM_ASSERT(gsl_comm->np == 1, "Wrong code path");

   int curr_elem;
   Vector shape;
   DenseMatrix elmat;
   Array<int> vdofs, vdofs_solid;
   for (int index = 0; index < npt; index++)
   {
      if (gsl_code[index] == 2) { continue; }
      // set int point for curr solid node in fluid mesh
      IntegrationPoint ip;
      ip.Set2(gsl_mfem_ref.GetData()+index*dim);
      if (dim == 3) { ip.z = gsl_mfem_ref(index*dim + 2); }

      // get fluid element the solid point is in
      curr_elem = (int) gsl_mfem_elem[index];
      fes_in->GetElementVDofs(curr_elem,vdofs);
      const FiniteElement &fe = *fes_in->GetFE(curr_elem);
      int nd = fe.GetDof();
      shape.SetSize(nd);         

      // evaluate shape funcs at current solid node
      fe.CalcShape(ip,shape);   

      // form local matrix
      elmat.SetSize(ncomp*nd,ncomp);
      elmat = 0.0;
      for (int a = 0; a < nd; a++)
      {
         elmat(a,0) = shape(a);
         if (ncomp > 1) elmat(a + nd,1) = shape(a);
         if (ncomp > 2) elmat(a + 2*nd,2) = shape(a);
      }

      // solid dofs
      vdofs_solid.SetSize(ncomp);
      vdofs_solid[0] = fes_out->DofToVDof(index,0);
      if (ncomp > 1) vdofs_solid[1] = fes_out->DofToVDof(index,1);
      if (ncomp > 2) vdofs_solid[2] = fes_out->DofToVDof(index,2);
      
      // Cmat->AddSubMatrix(vdofs_solid,vdofs,elmat,skip_zeros);
      C->SetSubMatrix(vdofs,vdofs_solid,elmat,skip_zeros);
   }
}
#else
void ImmersedFSIFindPts::AssembleC(FiniteElementSpace *fes_in, FiniteElementSpace *fes_out, const Vector &solid_nodes, int ordering_type)
{
   int ncomp   = fes_in->GetVDim(),
       nptorig = points_cnt,
       npt     = points_cnt;

   const int skip_zeros = 0;

   dim = fes_in->GetMesh()->Dimension();

   if (gsl_comm->np == 1)
   {
      int curr_elem;
      DofTransformation *doftrans;
      Vector shape;
      DenseMatrix elmat;
      Array<int> vdofs, vdofs_solid;
      auto par_fes_out = dynamic_cast<ParFiniteElementSpace*>(fes_out);
      int nn_solid = solid_nodes.Size() / ncomp;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] == 2) { continue; }
         // set int point for curr solid node in fluid mesh
         IntegrationPoint ip;
         ip.Set2(gsl_mfem_ref.GetData()+index*dim);
         if (dim == 3) { ip.z = gsl_mfem_ref(index*dim + 2); }

         // get fluid element the solid point is in
         curr_elem = (int) gsl_mfem_elem[index];
         doftrans = fes_in->GetElementVDofs(curr_elem,vdofs);
         const FiniteElement &fe = *fes_in->GetFE(curr_elem);
         int nd_in = fe.GetDof();
         shape.SetSize(nd_in);

         // evaluate shape funcs at current solid node
         fe.CalcShape(ip,shape);

         // form local matrix
         elmat.SetSize(ncomp*nd_in,ncomp);
         elmat = 0.0;
         for (int a = 0; a < nd_in; a++)
         {
            elmat(a,0) = shape(a);
            if (ncomp > 1) elmat(a + nd_in,1) = shape(a);
            if (ncomp > 2) elmat(a + 2*nd_in,2) = shape(a);
         }

         // pack output global tdofs
         int idx, idy, idz; // local tdof numbers = global tdof numbers for 1 proc
         vdofs_solid.SetSize(ncomp);
         idx = ordering_type == Ordering::byNODES ? index + 0*nn_solid : index*ncomp + 0;
         // vdofs_solid[0] = ltdof_to_ldof[idx];
         vdofs_solid[0] = par_fes_out->GetGlobalTDofNumber(ltdof_to_ldof[idx]);
         if (ncomp > 1)
         {
            idy = ordering_type == Ordering::byNODES ? index + 1*nn_solid : index*ncomp + 1;
            // vdofs_solid[1] = ltdof_to_ldof[idy];
            vdofs_solid[1] = par_fes_out->GetGlobalTDofNumber(ltdof_to_ldof[idy]);
         } 
         if (ncomp > 2) 
         {
            idz = ordering_type == Ordering::byNODES ? index + 2*nn_solid : index*ncomp + 2;
            // vdofs_solid[2] = ltdof_to_ldof[idz];
            vdofs_solid[2] = par_fes_out->GetGlobalTDofNumber(ltdof_to_ldof[idz]);
         }

         // assemble into IJ matrix
         int single_col[1];
         single_col[0] = 1;
         for (int i=0; i<nd_in*ncomp; ++i)
         {
            // global row number in the IJ matrix
            int global_row[1];
            if (dynamic_cast<ParFiniteElementSpace*>(fes_in))
            {
               global_row[0] = dynamic_cast<ParFiniteElementSpace*>(fes_in)->GetGlobalTDofNumber(vdofs[i]);
            }
            else
            {
               global_row[0] = vdofs[i];
            }
            for (int j=0; j<ncomp; ++j)
            {
               double val[1];
               val[0] = elmat(i,j);
               int global_col = vdofs_solid[j]; // global column number in the IJ matrix
               HYPRE_IJMatrixSetValues(ij_matrix, 1, single_col, global_row, &global_col, val);
               // HYPRE_IJMatrixAddToValues(ij_matrix, 1, single_col, global_row, &global_col, val);
            }
         }
      }
   }
   else
   {
      // Determine number of points to be sent
      int nptsend = 0;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] != 2) { nptsend +=1; }
      }

      // Pack data to send via crystal router
      struct gslib::array *outpt = new gslib::array;
      struct out_pt { double r[3]; int id[3]; uint index, el, proc; };
      struct out_pt *pt;
      array_init(struct out_pt, outpt, nptsend);
      outpt->n=nptsend;
      pt = (struct out_pt *)outpt->ptr;
      int nn_solid = solid_nodes.Size() / ncomp;
      auto par_fes_out = dynamic_cast<ParFiniteElementSpace*>(fes_out);
      int count_diff_proc = 0;
      for (int index = 0; index < npt; index++)
      {
         if (gsl_code[index] == 2) { continue; }
         pt->index = index;
         pt->proc  = gsl_proc[index];
         pt->el    = gsl_mfem_elem[index];

         if (pt->proc != par_fes_out->GetMyRank()) { count_diff_proc += 1; }

         for (int d = 0; d < dim; ++d) 
         { 
            pt->r[d]= gsl_mfem_ref(index*dim + d);
         }

         // pack output global tdofs
         int idx, idy, idz; // local tdof numbers
         idx = ordering_type == Ordering::byNODES ? index + 0*nn_solid : index*ncomp + 0;
         pt->id[0] = par_fes_out->GetGlobalTDofNumber(ltdof_to_ldof[idx]);
         if (ncomp > 1)
         {
            idy = ordering_type == Ordering::byNODES ? index + 1*nn_solid : index*ncomp + 1;
            pt->id[1] = par_fes_out->GetGlobalTDofNumber(ltdof_to_ldof[idy]);
         } 
         if (ncomp > 2) 
         {
            idz = ordering_type == Ordering::byNODES ? index + 2*nn_solid : index*ncomp + 2;
            pt->id[2] = par_fes_out->GetGlobalTDofNumber(ltdof_to_ldof[idz]);
         }
         
         ++pt;
      }

      // Transfer data to target MPI ranks
      sarray_transfer(struct out_pt, outpt, proc, 1, cr);

      npt = outpt->n;
      pt = (struct out_pt *)outpt->ptr;

      // store data for dofs not owned by this rank
      struct send_dof_data { int glob_row, glob_col[3]; double val[3]; uint proc; };
      Array<send_dof_data> send_dof_array(0);

      int curr_elem;
      Vector shape;
      DenseMatrix elmat;
      Array<int> vdofs;
      int single_col[1];
      single_col[0] = 1;
      for (int index = 0; index < npt; index++)
      {
         // get fluid element the solid point is in
         curr_elem = (int) pt->el;
         fes_in->GetElementVDofs(curr_elem,vdofs);
         const FiniteElement &fe = *fes_in->GetFE(curr_elem);
         int nd_in = fe.GetDof();
         shape.SetSize(nd_in);        

         // evaluate shape funcs at current solid node
         // set int point for curr solid node in fluid mesh
         IntegrationPoint ip;
         ip.Set3(&pt->r[0]);
         fe.CalcShape(ip,shape);     

         // form local matrix
         elmat.SetSize(ncomp*nd_in,ncomp);
         elmat = 0.0;
         for (int a = 0; a < nd_in; a++)
         {
            elmat(a,0) = shape(a);
            if (ncomp > 1) elmat(a + nd_in,1) = shape(a);
            if (ncomp > 2) elmat(a + 2*nd_in,2) = shape(a);
         }

         // assemble locally owned dofs into IJ matrix
         for (int i=0; i<nd_in*ncomp; ++i)
         {
            int local_dof = vdofs[i];
            int global_row[1];
            global_row[0] = dynamic_cast<ParFiniteElementSpace*>(fes_in)->GetGlobalTDofNumber(local_dof);

            // Set the value of the matrix from the owning rank (required by HYPRE_IJMatrixSetValues)
            if (dynamic_cast<ParFiniteElementSpace*>(fes_in)->GetLocalTDofNumber(local_dof) == -1) 
            {  
               // dof owned by another rank
               send_dof_data curr_dof;
               curr_dof.glob_row = global_row[0];
               curr_dof.proc = my_rank_gf(local_dof); // master rank for this dof
               curr_dof.glob_col[0] = pt->id[0];
               if (ncomp > 1) curr_dof.glob_col[1] = pt->id[1];
               if (ncomp > 2) curr_dof.glob_col[2] = pt->id[2];
               curr_dof.val[0] = elmat(i,0);
               if (ncomp > 1) curr_dof.val[1] = elmat(i,1);
               if (ncomp > 2) curr_dof.val[2] = elmat(i,2);
               
               send_dof_array.Append(curr_dof);
            }
            else 
            {
               // dof owned by this rank
               for (int j=0; j<ncomp; ++j)
               {
                  double val[1];
                  val[0] = elmat(i,j);
                  int global_col = pt->id[j]; // global column number in the IJ matrix

                  HYPRE_IJMatrixSetValues(ij_matrix, 1, single_col, global_row, &global_col, val);
                  // HYPRE_IJMatrixAddToValues(ij_matrix, 1, single_col, global_row, &global_col, val);
               }
            }            
         }

         ++pt;
      }

      // free memory
      array_free(outpt);
      delete outpt;

      // Pack data to send via crystal router for dofs not owned by this rank
      struct gslib::array *sendpt = new gslib::array;
      struct send_pt { int glob_row, glob_col[3]; double val[3]; uint proc; };
      struct send_pt *sdpt;
      array_init(struct send_pt, sendpt, send_dof_array.Size());
      sendpt->n=send_dof_array.Size();
      sdpt = (struct send_pt *)sendpt->ptr;

      for (int index = 0; index < send_dof_array.Size(); ++index)
      {
         sdpt->glob_row = send_dof_array[index].glob_row;
         sdpt->proc = send_dof_array[index].proc;
         sdpt->glob_col[0] = send_dof_array[index].glob_col[0];
         if (ncomp > 1) sdpt->glob_col[1] = send_dof_array[index].glob_col[1];
         if (ncomp > 2) sdpt->glob_col[2] = send_dof_array[index].glob_col[2];
         sdpt->val[0] = send_dof_array[index].val[0];
         if (ncomp > 1) sdpt->val[1] = send_dof_array[index].val[1];
         if (ncomp > 2) sdpt->val[2] = send_dof_array[index].val[2];

         ++sdpt;
      }

      // Transfer data for dofs not owned by this rank to the master rank
      sarray_transfer(struct send_pt, sendpt, proc, 1, cr);

      npt = sendpt->n;
      sdpt = (struct send_pt *)sendpt->ptr;

      // assemble into IJ matrix on the master rank
      for (int index = 0; index < npt; ++index)
      {
         // mfem::out << "Setting global row = " << sdpt->glob_row << std::endl;
         for (int j=0; j<ncomp; ++j)
         {
            HYPRE_IJMatrixSetValues(ij_matrix, 1, single_col, &sdpt->glob_row, &sdpt->glob_col[j], &sdpt->val[j]);
         }
         ++sdpt;
      }     

      // free memory
      array_free(sendpt);
      delete sendpt;      
   }
}
#endif

#ifdef MFEM_USE_MPI
void ImmersedFSIFindPts::SetupLocalDofToMasterRankMap(ParFiniteElementSpace *fes)
{
   // Per MFEM paper https://doi.org/10.1016/j.camwa.2020.06.009
   // the master rank is always the lowest rank that shares a given dof

   my_rank_gf.SetSpace(fes);
   my_rank_gf = fes->GetMyRank(); // set all dofs to my rank

   // communication
   GroupCommunicator &gcomm = fes->GroupComm();
   Array<double> dest_array(my_rank_gf.GetData(), my_rank_gf.Size());
   gcomm.Reduce<double>(dest_array, GroupCommunicator::Min);
   gcomm.Bcast(dest_array);

   // determine max row size for the hypre ij matrix
   DetermineMaxRowSize(fes);
}

void ImmersedFSIFindPts::DetermineMaxRowSize(ParFiniteElementSpace *fes)
{
   auto my_mesh = fes->GetMesh();
   auto my_fec = fes->FEColl();
   Array<Geometry::Type> geom_types;
   my_mesh->GetGeometries(my_mesh->SpaceDimension(),geom_types);
   int max_order = fes->GetMaxElementOrder();

   // Max dofs per element
   int max_dofs = 0;
   for (auto &geom : geom_types)
   {
      int ndofs = my_fec->GetNumDof(geom,max_order);
      if (ndofs > max_dofs) { max_dofs = ndofs; }
   }

   // Not sure here, but corner of hex shared by at most 8 elements, setting to 20 just in case
   max_rowsize[0] = max_dofs * 20; // max dofs per element * max elements sharing a dof (todo for second part)
}

void ImmersedFSIFindPts::SetupLocalTDofToLDofSolidMap(ParFiniteElementSpace *solid_fes)
{
   for (int i=0; i<solid_fes->GetVSize(); ++i)
   {
      int local_tdof = solid_fes->GetLocalTDofNumber(i);
      if (local_tdof == -1) { continue; } // dof not owned by this rank
      ltdof_to_ldof[local_tdof] = i;
   }
}
#endif

ImmersedFSIFindPts::~ImmersedFSIFindPts()
{
   delete C;
   delete CT;
}
#endif



