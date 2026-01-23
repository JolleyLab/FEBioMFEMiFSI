#include "MFEMNodalAverageOutputs.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

namespace {

template<typename Quantity>
bool ParseQuantityToken(const std::string &token,
                            const std::pair<const char*, Quantity> *table,
                            size_t count,
                            Quantity &out)
{
   const std::string norm = NormalizeNodalAverageToken(token);
   if (norm.empty()) return false;

   for (size_t i = 0; i < count; ++i)
   {
      if (norm == table[i].first)
      {
         out = table[i].second;
         return true;
      }
   }
   return false;
}

template<typename Quantity>
std::vector<Quantity> ParseQuantityList(const std::string &list,
                                            const std::pair<const char*, Quantity> *table,
                                            size_t count)
{
   std::vector<Quantity> result;
   std::stringstream ss(list);
   std::string token;
   while (std::getline(ss, token, ','))
   {
      Quantity q;
      if (ParseQuantityToken(token, table, count, q))
      {
         if (std::find(result.begin(), result.end(), q) == result.end())
         {
            result.push_back(q);
         }
      }
   }
   return result;
}

} // namespace

std::string NormalizeNodalAverageToken(const std::string &token)
{
   std::string norm;
   norm.reserve(token.size());
   for (char c : token)
   {
      if (c == ' ' || c == '\t') continue;
      if (c == '_' || c == '-') continue;
      norm.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
   }
   return norm;
}

bool ParseStressQuantityToken(const std::string &token, nodal_avg::StressQuantity &out)
{
   static const std::pair<const char*, nodal_avg::StressQuantity> table[] = {
      {"vonmises", nodal_avg::StressQuantity::VonMises},
      {"vm", nodal_avg::StressQuantity::VonMises},
      {"sxx", nodal_avg::StressQuantity::CauchyXX},
      {"cauchyxx", nodal_avg::StressQuantity::CauchyXX},
      {"syy", nodal_avg::StressQuantity::CauchyYY},
      {"cauchyyy", nodal_avg::StressQuantity::CauchyYY},
      {"szz", nodal_avg::StressQuantity::CauchyZZ},
      {"cauchyzz", nodal_avg::StressQuantity::CauchyZZ},
      {"sxy", nodal_avg::StressQuantity::CauchyXY},
      {"cauchyxy", nodal_avg::StressQuantity::CauchyXY},
      {"syz", nodal_avg::StressQuantity::CauchyYZ},
      {"cauchyyz", nodal_avg::StressQuantity::CauchyYZ},
      {"sxz", nodal_avg::StressQuantity::CauchyXZ},
      {"cauchyxz", nodal_avg::StressQuantity::CauchyXZ},
      {"principal1", nodal_avg::StressQuantity::Principal1},
      {"p1", nodal_avg::StressQuantity::Principal1},
      {"maxprincipal", nodal_avg::StressQuantity::Principal1},
      {"principal2", nodal_avg::StressQuantity::Principal2},
      {"p2", nodal_avg::StressQuantity::Principal2},
      {"midprincipal", nodal_avg::StressQuantity::Principal2},
      {"principal3", nodal_avg::StressQuantity::Principal3},
      {"p3", nodal_avg::StressQuantity::Principal3},
      {"minprincipal", nodal_avg::StressQuantity::Principal3}
   };
   return ParseQuantityToken(token, table, sizeof(table)/sizeof(table[0]), out);
}

bool ParseStrainQuantityToken(const std::string &token, nodal_avg::StrainQuantity &out)
{
   static const std::pair<const char*, nodal_avg::StrainQuantity> table[] = {
      {"vonmises", nodal_avg::StrainQuantity::VonMises},
      {"vm", nodal_avg::StrainQuantity::VonMises},
      {"exx", nodal_avg::StrainQuantity::Exx},
      {"eyy", nodal_avg::StrainQuantity::Eyy},
      {"ezz", nodal_avg::StrainQuantity::Ezz},
      {"exy", nodal_avg::StrainQuantity::Exy},
      {"eyz", nodal_avg::StrainQuantity::Eyz},
      {"exz", nodal_avg::StrainQuantity::Exz},
      {"principal1", nodal_avg::StrainQuantity::Principal1},
      {"p1", nodal_avg::StrainQuantity::Principal1},
      {"maxprincipal", nodal_avg::StrainQuantity::Principal1},
      {"principal2", nodal_avg::StrainQuantity::Principal2},
      {"p2", nodal_avg::StrainQuantity::Principal2},
      {"midprincipal", nodal_avg::StrainQuantity::Principal2},
      {"principal3", nodal_avg::StrainQuantity::Principal3},
      {"p3", nodal_avg::StrainQuantity::Principal3},
      {"minprincipal", nodal_avg::StrainQuantity::Principal3}
   };
   return ParseQuantityToken(token, table, sizeof(table)/sizeof(table[0]), out);
}

std::vector<nodal_avg::StressQuantity> ParseStressQuantityList(const std::string &list)
{
   static const std::pair<const char*, nodal_avg::StressQuantity> table[] = {
      {"vonmises", nodal_avg::StressQuantity::VonMises},
      {"vm", nodal_avg::StressQuantity::VonMises},
      {"sxx", nodal_avg::StressQuantity::CauchyXX},
      {"cauchyxx", nodal_avg::StressQuantity::CauchyXX},
      {"syy", nodal_avg::StressQuantity::CauchyYY},
      {"cauchyyy", nodal_avg::StressQuantity::CauchyYY},
      {"szz", nodal_avg::StressQuantity::CauchyZZ},
      {"cauchyzz", nodal_avg::StressQuantity::CauchyZZ},
      {"sxy", nodal_avg::StressQuantity::CauchyXY},
      {"cauchyxy", nodal_avg::StressQuantity::CauchyXY},
      {"syz", nodal_avg::StressQuantity::CauchyYZ},
      {"cauchyyz", nodal_avg::StressQuantity::CauchyYZ},
      {"sxz", nodal_avg::StressQuantity::CauchyXZ},
      {"cauchyxz", nodal_avg::StressQuantity::CauchyXZ},
      {"principal1", nodal_avg::StressQuantity::Principal1},
      {"p1", nodal_avg::StressQuantity::Principal1},
      {"maxprincipal", nodal_avg::StressQuantity::Principal1},
      {"principal2", nodal_avg::StressQuantity::Principal2},
      {"p2", nodal_avg::StressQuantity::Principal2},
      {"midprincipal", nodal_avg::StressQuantity::Principal2},
      {"principal3", nodal_avg::StressQuantity::Principal3},
      {"p3", nodal_avg::StressQuantity::Principal3},
      {"minprincipal", nodal_avg::StressQuantity::Principal3}
   };
   return ParseQuantityList(list, table, sizeof(table)/sizeof(table[0]));
}

std::vector<nodal_avg::StrainQuantity> ParseStrainQuantityList(const std::string &list)
{
   static const std::pair<const char*, nodal_avg::StrainQuantity> table[] = {
      {"vonmises", nodal_avg::StrainQuantity::VonMises},
      {"vm", nodal_avg::StrainQuantity::VonMises},
      {"exx", nodal_avg::StrainQuantity::Exx},
      {"eyy", nodal_avg::StrainQuantity::Eyy},
      {"ezz", nodal_avg::StrainQuantity::Ezz},
      {"exy", nodal_avg::StrainQuantity::Exy},
      {"eyz", nodal_avg::StrainQuantity::Eyz},
      {"exz", nodal_avg::StrainQuantity::Exz},
      {"principal1", nodal_avg::StrainQuantity::Principal1},
      {"p1", nodal_avg::StrainQuantity::Principal1},
      {"maxprincipal", nodal_avg::StrainQuantity::Principal1},
      {"principal2", nodal_avg::StrainQuantity::Principal2},
      {"p2", nodal_avg::StrainQuantity::Principal2},
      {"midprincipal", nodal_avg::StrainQuantity::Principal2},
      {"principal3", nodal_avg::StrainQuantity::Principal3},
      {"p3", nodal_avg::StrainQuantity::Principal3},
      {"minprincipal", nodal_avg::StrainQuantity::Principal3}
   };
   return ParseQuantityList(list, table, sizeof(table)/sizeof(table[0]));
}

std::string StressQuantityLabel(nodal_avg::StressQuantity q)
{
   switch (q)
   {
      case nodal_avg::StressQuantity::VonMises: return "VonMisesStress";
      case nodal_avg::StressQuantity::CauchyXX: return "Stress_xx";
      case nodal_avg::StressQuantity::CauchyYY: return "Stress_yy";
      case nodal_avg::StressQuantity::CauchyZZ: return "Stress_zz";
      case nodal_avg::StressQuantity::CauchyXY: return "Stress_xy";
      case nodal_avg::StressQuantity::CauchyYZ: return "Stress_yz";
      case nodal_avg::StressQuantity::CauchyXZ: return "Stress_xz";
      case nodal_avg::StressQuantity::Principal1: return "PrincipalStress_1";
      case nodal_avg::StressQuantity::Principal2: return "PrincipalStress_2";
      case nodal_avg::StressQuantity::Principal3: return "PrincipalStress_3";
      default: return "Stress";
   }
}

std::string StrainQuantityLabel(nodal_avg::StrainQuantity q)
{
   switch (q)
   {
      case nodal_avg::StrainQuantity::VonMises: return "VonMisesStrain";
      case nodal_avg::StrainQuantity::Exx: return "Strain_xx";
      case nodal_avg::StrainQuantity::Eyy: return "Strain_yy";
      case nodal_avg::StrainQuantity::Ezz: return "Strain_zz";
      case nodal_avg::StrainQuantity::Exy: return "Strain_xy";
      case nodal_avg::StrainQuantity::Eyz: return "Strain_yz";
      case nodal_avg::StrainQuantity::Exz: return "Strain_xz";
      case nodal_avg::StrainQuantity::Principal1: return "PrincipalStrain_1";
      case nodal_avg::StrainQuantity::Principal2: return "PrincipalStrain_2";
      case nodal_avg::StrainQuantity::Principal3: return "PrincipalStrain_3";
      default: return "Strain";
   }
}

namespace nodal_avg
{

double StressMeasureValue(const mat3ds &s, StressQuantity measure)
{
   switch (measure)
   {
      case StressQuantity::VonMises: return s.effective_norm();
      case StressQuantity::CauchyXX: return s.xx();
      case StressQuantity::CauchyYY: return s.yy();
      case StressQuantity::CauchyZZ: return s.zz();
      case StressQuantity::CauchyXY: return s.xy();
      case StressQuantity::CauchyYZ: return s.yz();
      case StressQuantity::CauchyXZ: return s.xz();
      case StressQuantity::Principal1:
      case StressQuantity::Principal2:
      case StressQuantity::Principal3:
      {
         double l1,l2,l3;
         SymmetricEigenvalues(s.xx(), s.yy(), s.zz(), s.xy(), s.xz(), s.yz(), l1,l2,l3);
         if (measure == StressQuantity::Principal1) return l1;
         if (measure == StressQuantity::Principal2) return l2;
         return l3;
      }
      default: return 0.0;
   }
}

double StrainMeasureValue(const mat3ds &E, StrainQuantity measure)
{
   switch (measure)
   {
      case StrainQuantity::VonMises: return E.effective_norm();
      case StrainQuantity::Exx: return E.xx();
      case StrainQuantity::Eyy: return E.yy();
      case StrainQuantity::Ezz: return E.zz();
      case StrainQuantity::Exy: return E.xy();
      case StrainQuantity::Eyz: return E.yz();
      case StrainQuantity::Exz: return E.xz();
      case StrainQuantity::Principal1:
      case StrainQuantity::Principal2:
      case StrainQuantity::Principal3:
      {
         double l1,l2,l3;
         SymmetricEigenvalues(E.xx(), E.yy(), E.zz(), E.xy(), E.xz(), E.yz(), l1,l2,l3);
         if (measure == StrainQuantity::Principal1) return l1;
         if (measure == StrainQuantity::Principal2) return l2;
         return l3;
      }
      default: return 0.0;
   }
}

} // namespace nodal_avg
