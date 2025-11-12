#include "linear_form_submesh.hpp"

using namespace mfem;

void LinearFormSubMesh::Assemble()
{
   Array<int> vdofs;
   ElementTransformation *eltrans;
   DofTransformation *doftrans;
   Vector elemvect;

   Vector::operator=(0.0);

   // The above operation is executed on device because of UseDevice().
   // The first use of AddElementVector() below will move it back to host
   // because both 'vdofs' and 'elemvect' are on host.

   if (fast_assembly && ext) { return ext->Assemble(); }

   if (domain_integs.Size())
   {
      for (int k = 0; k < domain_integs.Size(); k++)
      {
         if (domain_integs_marker[k] != NULL)
         {
            MFEM_VERIFY(domain_integs_marker[k]->Size() ==
                        (fes->GetMesh()->attributes.Size() ?
                           fes->GetMesh()->attributes.Max() : 0),
                        "invalid element marker for domain linear form "
                        "integrator #" << k << ", counting from zero");
         }
      }

      for (int i = 0; i < fes -> GetNE(); i++)
      {
         int elem_attr = fes->GetMesh()->GetAttribute(i);
         for (int k = 0; k < domain_integs.Size(); k++)
         {
            const Array<int> * const markers = domain_integs_marker[k];
            if ( markers == NULL || (*markers)[elem_attr-1] == 1 )
            {
               doftrans = fes -> GetElementVDofs (i, vdofs);
               eltrans = fes -> GetElementTransformation (i);
               domain_integs[k]->AssembleRHSElementVect(*fes->GetFE(i),
                                                         *eltrans, elemvect);
               if (doftrans)
               {
                  doftrans->TransformDual(elemvect);
               }
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }
   AssembleDelta();

   if (boundary_integs.Size())
   {
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker;
      if (global_bdr_attr_size == -1)
      {
         bdr_attr_marker.SetSize(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      }
      else
      {
         bdr_attr_marker.SetSize(global_bdr_attr_size);
      }

      bdr_attr_marker = 0;
      for (int k = 0; k < boundary_integs.Size(); k++)
      {
         if (boundary_integs_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *boundary_integs_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (int i = 0; i < fes -> GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }
         doftrans = fes -> GetBdrElementVDofs (i, vdofs);
         eltrans = fes -> GetBdrElementTransformation (i);
         for (int k=0; k < boundary_integs.Size(); k++)
         {
            if (boundary_integs_marker[k] &&
                  (*boundary_integs_marker[k])[bdr_attr-1] == 0) { continue; }

            boundary_integs[k]->AssembleRHSElementVect(*fes->GetBE(i),
                                                         *eltrans, elemvect);

            if (doftrans)
            {
               doftrans->TransformDual(elemvect);
            }
            AddElementVector (vdofs, elemvect);
         }
      }
   }
   if (boundary_face_integs.Size())
   {
      FaceElementTransformations *tr;
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker;
      if (global_bdr_attr_size == -1)
      {
         bdr_attr_marker.SetSize(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      }
      else
      {
         bdr_attr_marker.SetSize(global_bdr_attr_size);
      }
      
      bdr_attr_marker = 0;
      for (int k = 0; k < boundary_face_integs.Size(); k++)
      {
         if (boundary_face_integs_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *boundary_face_integs_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary face integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (int i = 0; i < mesh->GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

         tr = mesh->GetBdrFaceTransformations(i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            for (int k = 0; k < boundary_face_integs.Size(); k++)
            {
               if (boundary_face_integs_marker[k] &&
                     (*boundary_face_integs_marker[k])[bdr_attr-1] == 0)
               { continue; }

               boundary_face_integs[k]->
               AssembleRHSElementVect(*fes->GetFE(tr->Elem1No),
                                       *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }

   if (interior_face_integs.Size())
   {
      Mesh *mesh = fes->GetMesh();

      for (int k = 0; k < interior_face_integs.Size(); k++)
      {
         for (int i = 0; i < mesh->GetNumFaces(); i++)
         {
            FaceElementTransformations *tr = NULL;
            tr = mesh->GetInteriorFaceTransformations (i);
            if (tr != NULL)
            {
               fes -> GetElementVDofs (tr -> Elem1No, vdofs);
               Array<int> vdofs2;
               fes -> GetElementVDofs (tr -> Elem2No, vdofs2);
               vdofs.Append(vdofs2);
               interior_face_integs[k]->
               AssembleRHSElementVect(*fes->GetFE(tr->Elem1No),
                                       *fes->GetFE(tr->Elem2No),
                                       *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }
}

#ifdef MFEM_USE_MPI
void ParLinearFormSubMesh::Assemble()
{
   Array<int> vdofs;
   ElementTransformation *eltrans;
   DofTransformation *doftrans;
   Vector elemvect;

   Vector::operator=(0.0);

   // The above operation is executed on device because of UseDevice().
   // The first use of AddElementVector() below will move it back to host
   // because both 'vdofs' and 'elemvect' are on host.

   if (fast_assembly && ext) { return ext->Assemble(); }

   if (domain_integs.Size())
   {
      for (int k = 0; k < domain_integs.Size(); k++)
      {
         if (domain_integs_marker[k] != NULL)
         {
            MFEM_VERIFY(domain_integs_marker[k]->Size() ==
                        (fes->GetMesh()->attributes.Size() ?
                           fes->GetMesh()->attributes.Max() : 0),
                        "invalid element marker for domain linear form "
                        "integrator #" << k << ", counting from zero");
         }
      }

      for (int i = 0; i < fes -> GetNE(); i++)
      {
         int elem_attr = fes->GetMesh()->GetAttribute(i);
         for (int k = 0; k < domain_integs.Size(); k++)
         {
            const Array<int> * const markers = domain_integs_marker[k];
            if ( markers == NULL || (*markers)[elem_attr-1] == 1 )
            {
               doftrans = fes -> GetElementVDofs (i, vdofs);
               eltrans = fes -> GetElementTransformation (i);
               domain_integs[k]->AssembleRHSElementVect(*fes->GetFE(i),
                                                         *eltrans, elemvect);
               if (doftrans)
               {
                  doftrans->TransformDual(elemvect);
               }
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }
   AssembleDelta();

   if (boundary_integs.Size())
   {
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker;
      if (global_bdr_attr_size == -1)
      {
         bdr_attr_marker.SetSize(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      }
      else
      {
         bdr_attr_marker.SetSize(global_bdr_attr_size);
      }

      bdr_attr_marker = 0;
      for (int k = 0; k < boundary_integs.Size(); k++)
      {
         if (boundary_integs_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *boundary_integs_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (int i = 0; i < fes -> GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }
         doftrans = fes -> GetBdrElementVDofs (i, vdofs);
         eltrans = fes -> GetBdrElementTransformation (i);
         for (int k=0; k < boundary_integs.Size(); k++)
         {
            if (boundary_integs_marker[k] &&
                  (*boundary_integs_marker[k])[bdr_attr-1] == 0) { continue; }

            boundary_integs[k]->AssembleRHSElementVect(*fes->GetBE(i),
                                                         *eltrans, elemvect);

            if (doftrans)
            {
               doftrans->TransformDual(elemvect);
            }
            AddElementVector (vdofs, elemvect);
         }
      }
   }
   if (boundary_face_integs.Size())
   {
      FaceElementTransformations *tr;
      Mesh *mesh = fes->GetMesh();

      // Which boundary attributes need to be processed?
      Array<int> bdr_attr_marker;
      if (global_bdr_attr_size == -1)
      {
         bdr_attr_marker.SetSize(mesh->bdr_attributes.Size() ?
                                 mesh->bdr_attributes.Max() : 0);
      }
      else
      {
         bdr_attr_marker.SetSize(global_bdr_attr_size);
      }
      
      bdr_attr_marker = 0;
      for (int k = 0; k < boundary_face_integs.Size(); k++)
      {
         if (boundary_face_integs_marker[k] == NULL)
         {
            bdr_attr_marker = 1;
            break;
         }
         Array<int> &bdr_marker = *boundary_face_integs_marker[k];
         MFEM_ASSERT(bdr_marker.Size() == bdr_attr_marker.Size(),
                     "invalid boundary marker for boundary face integrator #"
                     << k << ", counting from zero");
         for (int i = 0; i < bdr_attr_marker.Size(); i++)
         {
            bdr_attr_marker[i] |= bdr_marker[i];
         }
      }

      for (int i = 0; i < mesh->GetNBE(); i++)
      {
         const int bdr_attr = mesh->GetBdrAttribute(i);
         if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

         tr = mesh->GetBdrFaceTransformations(i);
         if (tr != NULL)
         {
            fes -> GetElementVDofs (tr -> Elem1No, vdofs);
            for (int k = 0; k < boundary_face_integs.Size(); k++)
            {
               if (boundary_face_integs_marker[k] &&
                     (*boundary_face_integs_marker[k])[bdr_attr-1] == 0)
               { continue; }

               boundary_face_integs[k]->
               AssembleRHSElementVect(*fes->GetFE(tr->Elem1No),
                                       *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }

   if (interior_face_integs.Size())
   {
      Mesh *mesh = fes->GetMesh();

      for (int k = 0; k < interior_face_integs.Size(); k++)
      {
         for (int i = 0; i < mesh->GetNumFaces(); i++)
         {
            FaceElementTransformations *tr = NULL;
            tr = mesh->GetInteriorFaceTransformations (i);
            if (tr != NULL)
            {
               fes -> GetElementVDofs (tr -> Elem1No, vdofs);
               Array<int> vdofs2;
               fes -> GetElementVDofs (tr -> Elem2No, vdofs2);
               vdofs.Append(vdofs2);
               interior_face_integs[k]->
               AssembleRHSElementVect(*fes->GetFE(tr->Elem1No),
                                       *fes->GetFE(tr->Elem2No),
                                       *tr, elemvect);
               AddElementVector (vdofs, elemvect);
            }
         }
      }
   }

   // parallel part
   if (interior_face_integs.Size())
   {
      pfes->ExchangeFaceNbrData();
      AssembleSharedFaces();
   }
}
#endif