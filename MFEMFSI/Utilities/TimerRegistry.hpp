#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cassert>

namespace mfem { class StopWatch; }

class TimerRegistry
{
public:
   struct Entry
   {
      std::string key;          // unique id (e.g., "timestep")
      std::string label;        // printed name (e.g., "TIMESTEP")
      mfem::StopWatch* sw;      // stopwatch pointer
      std::string percent_of;   // key used for % calc (e.g., "timestep")
      std::string parent;       // optional parent key for grouping/indent
      double restart_time = 0.0;// collected time prior to restart
      int indent = 0;           // printed indent level
      bool enabled = true;      // allow toggling at runtime
      bool show_percent = true;  // whether to show percent column
   };

   void Clear()
   {
      entries.clear();
      index.clear();
   }

   // Add an entry; returns index
   int Add(const std::string& key,
           const std::string& label,
           mfem::StopWatch* sw,
           const std::string& percent_of_key,
           const std::string& parent_key = "",
           int indent = 0,
           bool show_percent = true,
           bool enabled = true)
   {
      assert(index.count(key) == 0 && "Duplicate timer key");
         // Initialize fields in the same order as declared in Entry so that
         // restart_time is not accidentally set from the indent parameter.
         entries.push_back({key,            // key
             label,          // label
             sw,             // stopwatch
             percent_of_key, // percent_of
             parent_key,     // parent
             0.0,            // restart_time (can be set later)
             indent,         // indent
             enabled,        // enabled
             show_percent}); // show_percent
      int idx = static_cast<int>(entries.size() - 1);
      index[key] = idx;
      return idx;
   }

   void SetEnabled(const std::string& key, bool enabled)
   {
      auto it = index.find(key);
      if (it != index.end()) entries[it->second].enabled = enabled;
   }

   // Set the stored restart_time for a given timer key. This value is
   // added to the live stopwatch time in CollectSeconds(), allowing us
   // to carry timing information across restarts.
   void SetRestartTime(const std::string& key, double restart_time)
   {
      auto it = index.find(key);
      if (it == index.end()) { return; }
      entries[it->second].restart_time = restart_time;
   }

   const std::vector<Entry>& Entries() const { return entries; }

   // Collect raw seconds in same order as Entries()
   std::vector<double> CollectSeconds() const;

   // Utility to look up a time by key from a times vector (same order as Entries)
   double TimeOf(const std::string& key, const std::vector<double>& times) const
   {
      auto it = index.find(key);
      if (it == index.end()) return 0.0;
      return times[it->second];
   }

   // return number of timers
   int NumTimers() const { return static_cast<int>(entries.size()); }

private:
   std::vector<Entry> entries;
   std::unordered_map<std::string,int> index;
};