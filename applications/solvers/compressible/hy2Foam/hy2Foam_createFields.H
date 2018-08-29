#include "createTimeControls.H"
#include "createThermoFields.H"   
#include "readTimeControls.H"

// --- Initialisation of reacting and two-temperature fields
#include "createReactingFields.H"
#include "createVibrationalFields.H"

#include "createFvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "numerics/readFluxScheme.H"

// --- Write fields in the 0 folder
#include "write/write0.H" 

dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
/*dimensionedScalar vv_zero("vv_zero", dimArea, 0.0);
dimensionedVector Vv_one("vv_one", dimArea, vector::one);
dimensionedScalar vv_one("vv_one", dimArea, 1.0);
dimensionedScalar vvv_zero("vvv_zero", dimMass/dimTime, 0.0);
dimensionedScalar vvvv_zero("vvvv_zero", dimless, 0.0);*/

// --- Upwind interpolation of primitive fields on faces
#include "numerics/upwindInterpolation.H"
#include "numerics/fluxesCalculation.H"

// --- Time control and local time stepping
#include "numerics/compressibleCourantNo.H"
#include "LTS/setInitialrDeltaT.H"
#include "setInitialDeltaT.H"
