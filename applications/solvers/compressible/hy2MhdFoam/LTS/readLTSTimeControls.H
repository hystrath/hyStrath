// Smoothing parameter (0-1) when smoothing iterations > 0
scalar rDeltaTSmoothingCoeff
(
    runTime.controlDict().lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 1)
);

// Damping coefficient (1-0)
scalar rDeltaTDampingCoeff
(
    runTime.controlDict().lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1)
);

// Maximum change in cell temperature per iteration (relative to previous value)
scalar alphaTemp(runTime.controlDict().lookupOrDefault("alphaTemp", 0.05));
