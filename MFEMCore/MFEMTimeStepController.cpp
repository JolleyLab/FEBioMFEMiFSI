#include "MFEMTimeStepController.h"

#include "../MFEMFSI/IncompressibleFluid/IncompressibleFluidSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

double MFEMConvectiveCFLTimeStepEstimator::ComputeCandidate(double current_dt) const
{
   return solver->estimate_convective_cfl_dt();
}

MFEMTimeStepController::MFEMTimeStepController() {}

MFEMTimeStepController::~MFEMTimeStepController()
{
	ClearEstimators();
	std::vector<std::unique_ptr<MFEMTimeStepEstimator>>().swap(estimators_);
}

double MFEMTimeStepController::ComputeTimeStep(double dt)
{
	double candidate = std::numeric_limits<double>::infinity();
	for (const auto &est : estimators_)
   {
		if (!est) continue;
		double v = est->ComputeCandidate(dt);
		v *= est->GetFactor();
		if (std::isfinite(v) && v > 0.0) 
      {
			candidate = std::min(candidate, v);
		}
	}

	// If no valid candidates, fall back to input dt
	if (!std::isfinite(candidate)) 
   {
		candidate = dt;
	}

	// Apply optional clamping
	if (clamp_min_ > 0.0) candidate = std::max(candidate, clamp_min_);
	if (clamp_max_ > 0.0) candidate = std::min(candidate, clamp_max_);

	// Apply optional safety factor
	candidate *= safety_factor_;

	return std::min(candidate,dt);
}

void MFEMTimeStepController::AddEstimator(std::unique_ptr<MFEMTimeStepEstimator> estimator)
{
	if (estimator) estimators_.emplace_back(std::move(estimator));
}

void MFEMTimeStepController::AddEstimatorFn(std::function<double(double)> estimator_fn, double C)
{
	if (!estimator_fn) return;
	auto est = std::make_unique<MFEMFunctionTimeStepEstimator>(std::move(estimator_fn), C);
	estimators_.emplace_back(std::move(est));
}

void MFEMTimeStepController::ClearEstimators()
{
	estimators_.clear();
}

void MFEMTimeStepController::SetSafetyFactor(double factor)
{
	safety_factor_ = (factor > 0.0 ? factor : 1.0);
}

double MFEMTimeStepController::GetSafetyFactor() const
{
	return safety_factor_;
}

void MFEMTimeStepController::SetClamp(double dt_min, double dt_max)
{
	clamp_min_ = dt_min;
	clamp_max_ = dt_max;
}