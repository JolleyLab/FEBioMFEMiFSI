#pragma once
#include <FECore/FEModelParam.h>
#include <FECore/FEMaterialPoint.h>
#include <mfem.hpp>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/FELinearSystem.h>
#include <FECore/CompactUnSymmMatrix.h>

class MFEMImmersedElasticDomain;

class FEScalarCoefficient : public mfem::Coefficient
{
private:
	FEModel* m_fem = nullptr;
	FEParamDouble& m_val;

public:
	FEScalarCoefficient(FEModel* fem, FEParamDouble& val) : m_fem(fem), m_val(val) {}

	double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;
};

class FEVectorCoefficient : public mfem::VectorCoefficient
{
private:
	FEModel* m_fem = nullptr;
	FEParamDouble& m_vx;
	FEParamDouble& m_vy;
	FEParamDouble& m_vz;

public:
	FEVectorCoefficient(FEModel* fem, FEParamDouble& vx, FEParamDouble& vy, FEParamDouble& vz) : m_fem(fem), m_vx(vx), m_vy(vy), m_vz(vz), VectorCoefficient(3) {}

	void Eval(mfem::Vector& V, mfem::ElementTransformation& T,
		const mfem::IntegrationPoint& ip) override;
};

class CustomStreamBuffer : public std::streambuf 
{
public:
   CustomStreamBuffer() {};
   CustomStreamBuffer(std::streambuf* buf, FEModel *fe_model_) : originalBuffer(buf), fe_model(fe_model_) {};

   FEModel* GetFEModel() const { return fe_model; };

   std::streambuf* originalBuffer = nullptr;
   FEModel* fe_model = nullptr;
   
protected:
   virtual int overflow(int c) override {
      if (c != EOF) {
         char ch = static_cast<char>(c);
         buffer << ch;
         if (ch == '\n') {
               std::string str = buffer.str();
               // Escape percent symbols
               size_t pos = 0;
               while ((pos = str.find('%', pos)) != std::string::npos) {
                  str.insert(pos, "%");
                  pos += 2; // Move past the escaped percent symbol
               }
               feLog("%s", str.c_str());
               buffer.str("");
               buffer.clear();
         }
      }
      return c;
   }

private:
   std::ostringstream buffer;
};

class MFEMGlobalVector : public FEGlobalVector
{
public:
   // constructor
   MFEMGlobalVector(FEModel* fem, int size, std::vector<double>& R, std::vector<double>& Fr);
   // destructor
   ~MFEMGlobalVector() noexcept override;

   //! Assemble the element vector into this global vector
	void Assemble(std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe, bool bdom = false) override;

	//! Assemble into this global vector
	void Assemble(std::vector<int>& lm, std::vector<double>& fe) override;

	//! assemble a nodel value
	void Assemble(int node, int dof, double f) override;

   // Assemble into global mfem vector
   void Assemble(mfem::Array<int> &vdofs, std::vector<double>& fe);
   void Assemble(mfem::Array<int> &vdofs, mfem::Vector& fe);

   mfem::Vector &GetMfemVector() { return mfem_R; }

public:
   mfem::Vector mfem_R;
};

class MFEMGlobalMatrix : public FEGlobalMatrix
{
public:
   MFEMGlobalMatrix(CRSSparseMatrix* pK, bool del = true);

   void CreateProfile(FEModel* pfem, int neq, bool reshape);
   void BuildMatrixProfile(FEModel* pfem, bool breset);

   mfem::SparseMatrix *GetMfemMatrix() { return mfem_K; }

private:
   mfem::SparseMatrix *mfem_K = nullptr;
};
