thermalRelaxModelVT().correct();
if(noSolvedVibEq > 1)
{
    thermalRelaxModelVV().correct();
}

if(containsElectrons)
{
    thermalRelaxModelHE().correct();
    thermalRelaxModeleV().correct();
}

// --- Solve vibrational energy equations, INVISCID
forAll(rhoEvk, speciei) // NEW VINCENT 05/08/2016
{
    fvScalarMatrix rhoEvkInviscid
    (
        fvm::ddt(rhoEvk[speciei]) 
      + fvc::div(phiEvk[speciei]) 
      - thermalRelaxModelVT().QVT(speciei)
      - reaction->Scv(speciei)
    );
    
    if(noSolvedVibEq > 1)
    {
        rhoEvkInviscid -= thermalRelaxModelVV().QVV(speciei);
    }
    
    if(containsElectrons)
    { 
        rhoEvkInviscid -= thermalRelaxModeleV().QeV(speciei)
            + reaction->Seiir(speciei);
            
        if(composition.vibTempAssociativity("e-") == speciei) 
        {
            rhoEvkInviscid += (fvc::grad(composition.pP("e-")) & U)
                - thermalRelaxModelHE().QHE();
        }   
    }
    
    rhoEvkInviscid.solve
    (
        mesh.solver("rhoEv_" + composition.species()[speciei])
    );

    forAll(Y[speciei], celli)
    {
        if(Y[speciei][celli] > 1.0e-6)
        {
            hevk[speciei][celli] = rhoEvk[speciei][celli]/pD[speciei][celli];
        }
        else
        {
            hevk[speciei][celli] = 0.0;
        }
    } 
       
    hevk[speciei].correctBoundaryConditions();
    rhoEvk[speciei].boundaryField() = pD[speciei].boundaryField()*hevk[speciei].boundaryField();
}

thermo.correctNonMoleculeHEVK();

rhoEv *= 0.0;
forAll(hevk, speciei)
{
    rhoEv += pD[speciei]*hevk[speciei]; // NEW VINCENT 05/08/2016
}

ev = rhoEv/rho;

