
p.dimensionedInternalField() = rho.dimensionedInternalField() / psi.dimensionedInternalField(); 
p.correctBoundaryConditions();
rho.boundaryField() = psi.boundaryField()*p.boundaryField();
