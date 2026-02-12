#pragma once

#include <functional>
#include <memory>
#include <vector>

namespace ifs_solver
{
   class IncompressibleFluidSolver;
}

// Interface for individual time step estimators
class MFEMTimeStepEstimator 
{
public:
   virtual ~MFEMTimeStepEstimator() = default;
   virtual double ComputeCandidate(double current_dt) const = 0;

   // Per-estimator scaling factor (CFL-like). Defaults to 1.0.
   void SetFactor(double c) { factor_ = (c > 0.0 ? c : 1.0); }
   double GetFactor() const { return factor_; }

protected:
   double factor_ = 1.0;
};

// Convective CFL-based estimator
class MFEMConvectiveCFLTimeStepEstimator : public MFEMTimeStepEstimator
{
public:
   MFEMConvectiveCFLTimeStepEstimator(ifs_solver::IncompressibleFluidSolver* solver_, double cfl)
      : solver(solver_) { SetFactor(cfl); }
   
   double ComputeCandidate(double current_dt) const override;

private:
   ifs_solver::IncompressibleFluidSolver* solver;
};

// A estimator that wraps a function and applies a factor.
class MFEMFunctionTimeStepEstimator : public MFEMTimeStepEstimator 
{
public:
   explicit MFEMFunctionTimeStepEstimator(std::function<double(double)> fn, double c = 1.0)
      : fn_(std::move(fn)) { SetFactor(c); }

   double ComputeCandidate(double current_dt) const override 
   {
      if (!fn_) return current_dt;
      double raw = fn_(current_dt);
      return raw;
   }

private:
   std::function<double(double)> fn_;
};

class MFEMTimeStepController
{
public:
   MFEMTimeStepController();
   ~MFEMTimeStepController();

   // Aggregates all registered estimators and returns the minimum candidate.
   // If no estimators are registered, returns the input dt.
   virtual double ComputeTimeStep(double dt);

   // Register an estimator object (strategy pattern).
   void AddEstimator(std::unique_ptr<MFEMTimeStepEstimator> estimator);

   // Register a simple functional rule returning a candidate dt with per-estimator factor C.
   void AddEstimatorFn(std::function<double(double)> estimator_fn, double C = 1.0);

   // Remove all estimators.
   void ClearEstimators();

   // Optional safety factor applied to the final dt (<= 1.0 typical).
   void SetSafetyFactor(double factor);
   double GetSafetyFactor() const;

   // Optional clamping to [dt_min, dt_max]. Use non-positive to disable.
   void SetClamp(double dt_min, double dt_max);

private:
   double safety_factor_ = 1.0;
   double clamp_min_ = 0.0;
   double clamp_max_ = 0.0;

   std::vector<std::unique_ptr<MFEMTimeStepEstimator>> estimators_;
};
