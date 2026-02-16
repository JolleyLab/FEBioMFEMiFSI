#pragma once

#include <FECore/FECallback.h>

// run post-solve actions when the model is solved
bool timing_callback(FEModel* fem, unsigned int nreason, void* pd);

// Adjust time step size
bool time_step_callback(FEModel* fem, unsigned int nreason, void* pd);