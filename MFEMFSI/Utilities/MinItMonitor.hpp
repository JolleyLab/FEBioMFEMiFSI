#pragma once
#include "mfem.hpp"

class MinItMonitor : public mfem::IterativeSolverMonitor
{
private:
   int min_it = 1;

public:
   MinItMonitor(int min_it_) : min_it(min_it_) { };

   void MonitorResidual(int it, mfem::real_t norm, const mfem::Vector &r, bool final) override;
};