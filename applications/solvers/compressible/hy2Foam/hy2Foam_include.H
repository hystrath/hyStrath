#include <map>

// --- Default libraries
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "fvOptions.H"

// --- Smoothing for Local-Time-Stepping
#include "fvcSmooth.H"

// --- High-temperature chemistry library
#include "rho2HTCModel.H"

// --- Relaxation processes
#include "relaxationTimeModel.H"
#include "relaxationTimeModelVV.H"
#include "relaxationTimeModelHE.H"
#include "relaxationTimeModeleV.H"

// --- Multi-species transport models
#include "mixingRule.H"
#include "multiSpeciesTransportModel.H"

// --- Rarefied gas library
#include "rarefactionParameter.H"

// --- Numerical flux functions
#include "numerics/directionInterpolate.H"
#include "numerics/AUSM_functions.H"
#include "numerics/KNP-KT.H"
#include "numerics/KFVS.H"
#include "numerics/LxF.H"
#include "numerics/ROE.H"

// --- Others
#include "wallFvPatch.H"
