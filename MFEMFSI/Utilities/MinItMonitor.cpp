#include "MinItMonitor.hpp"

using namespace mfem;

void MinItMonitor::MonitorResidual(int it, real_t norm, const Vector &r, bool final)
{
   if (it < min_it)
   {
      this->converged = false;
   }
   else
   {
      return;
   }
}