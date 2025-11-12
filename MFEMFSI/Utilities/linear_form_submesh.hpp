#pragma once

#include "mfem.hpp"

class LinearFormSubMesh : public mfem::LinearForm
{
protected:
   /// Size of bdr attribute array to use. Useful when dealing with submeshes
   int global_bdr_attr_size = -1;

public:
   /// Constructor
   LinearFormSubMesh() : mfem::LinearForm() {};
   LinearFormSubMesh(mfem::FiniteElementSpace *f) : mfem::LinearForm(f) {};

   void SetGlobalBdrAttrSize(int val_) { global_bdr_attr_size = val_;};

   void Assemble();
};

#ifdef MFEM_USE_MPI
class ParLinearFormSubMesh : public mfem::ParLinearForm
{
protected:
   /// Size of bdr attribute array to use. Useful when dealing with submeshes
   int global_bdr_attr_size = -1;

public:
   /// Constructor
   ParLinearFormSubMesh() : mfem::ParLinearForm() {};
   ParLinearFormSubMesh(mfem::ParFiniteElementSpace *f) : mfem::ParLinearForm(f) {};

   void SetGlobalBdrAttrSize(int val_) { global_bdr_attr_size = val_;};

   void Assemble();
};
#endif