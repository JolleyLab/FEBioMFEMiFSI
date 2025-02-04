#pragma once
#include <FECore/FEModelParam.h>
#include <FECore/FEMaterialPoint.h>
#include <mfem.hpp>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <FECore/log.h>
#include <FECore/FEModel.h>

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
               feLog("%s", buffer.str().c_str());
               buffer.str("");
               buffer.clear();
         }
      }
      return c;
   }

private:
   std::ostringstream buffer;
};