#pragma once

namespace ifs_solver
{
// viscous term type
enum ViscousTermType {LaplaceForm = 0, StressForm = 1};

// Dimensional flag
enum dim_flag {NonDimensional = 0, Dimensional = 1};

// Nonlinear form type
enum NonlinearFormType {Convective = 0, EMAC = 1};

// Preconditioner type
enum PreconditionerType {SIMPLE = 0, SCR = 1, aSCR = 2, ILU = 3, Jacobi = 4};

enum StabilizationType 
{
   NONE = 0,
   SUPG = 1,
   SUPG_PSPG = 2,
   VMS = 3
};

// WSS averaging type
enum WSSAvgType {Standard = 0, LumpedL2 = 1};

// Time Stepping Method
enum TimeSteppingMethod {GeneralizedAlpha = 0, BackwardEuler = 1};

// Mapped Inflow BC type
enum inflowTYPE {Parabolic = 0, Plug = 1};

} // namespace ifs_solver