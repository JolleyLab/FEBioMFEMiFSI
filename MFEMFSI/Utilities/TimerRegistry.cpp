#include "TimerRegistry.hpp"
#include "mfem.hpp"

std::vector<double> TimerRegistry::CollectSeconds() const
{
   std::vector<double> secs;
   secs.reserve(entries.size());
   for (const auto& e : entries)
   {
      secs.push_back(e.enabled && e.sw ? e.sw->RealTime() : 0.0);
   }
   return secs;
}