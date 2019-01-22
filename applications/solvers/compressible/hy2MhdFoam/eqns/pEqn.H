
p.ref() = rho() / psi(); 
p.correctBoundaryConditions();
rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField();
