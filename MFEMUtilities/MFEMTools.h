#pragma once
#include <cmath>
#include <type_traits>
#include <FECore/FEModelParam.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <mfem.hpp>

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
	FEVectorCoefficient(FEModel* fem, FEParamDouble& vx, FEParamDouble& vy, FEParamDouble& vz) : m_fem(fem), m_vx(vx), m_vy(vy), m_vz(vz), mfem::VectorCoefficient(3) {}

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

// Truncate a floating-point value to a fixed number of decimal places.
// This performs truncation (not rounding): values are cut off beyond 'places'.
// - Leaves NaN/Inf unchanged.
// Example: TruncateToDigits(3.14159, 2) -> 3.14
template <typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, T>::type
TruncateToDigits(T value, int places)
{
   if (!std::isfinite(value)) { return value; }
   if (places <= 0) { return std::trunc(value); }

   // Limit places to avoid overflow in pow for extreme inputs
   const int max_places = 18; // safe for double; float effectively smaller
   const int p = (places > max_places) ? max_places : places;

   const T factor = static_cast<T>(std::pow(10.0, static_cast<double>(p)));
   const T scaled = value * factor;
   const T truncated = std::trunc(scaled);
   return truncated / factor;
}

template float TruncateToDigits<float>(float value, int places);
template double TruncateToDigits<double>(double value, int places);