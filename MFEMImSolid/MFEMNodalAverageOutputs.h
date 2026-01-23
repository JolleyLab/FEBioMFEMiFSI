#pragma once

#include <string>
#include <vector>
#include <FECore/mat3d.h>

#include "../MFEMFSI/ImmersedFSI/NodalAverageOutputs.h"

namespace nodal_avg
{
   // FEBio mat3ds overloads for stress/strain measures
   double StressMeasureValue(const mat3ds &s, StressQuantity measure);
   double StrainMeasureValue(const mat3ds &E, StrainQuantity measure);
}

std::string NormalizeNodalAverageToken(const std::string &token);
bool ParseStressQuantityToken(const std::string &token, nodal_avg::StressQuantity &out);
bool ParseStrainQuantityToken(const std::string &token, nodal_avg::StrainQuantity &out);
std::vector<nodal_avg::StressQuantity> ParseStressQuantityList(const std::string &list);
std::vector<nodal_avg::StrainQuantity> ParseStrainQuantityList(const std::string &list);
std::string StressQuantityLabel(nodal_avg::StressQuantity q);
std::string StrainQuantityLabel(nodal_avg::StrainQuantity q);