thermo.correct();
thermo.correct2T();
transport().correct();

thermalRelaxModelVT().correct();
thermalRelaxModelVV().correct();
thermalRelaxModelHE().correct();

dimensionedScalar epsilon("SMALL", dimMass/dimVolume, Foam::SMALL);

// --- Solve vibrational energy equations, INVISCID
forAll (rhoEvk_mode, speciei)
{
    forAll (rhoEvk_mode[speciei], vibMode)
    {
        rhoEvk_mode[speciei][vibMode] = pD[speciei]*hevk_mode[speciei][vibMode];

        if (composition.particleType(speciei) == 2)
        {
            if (containsElectrons)
            {
                solve
                (
                    fvm::ddt(rhoEvk_mode[speciei][vibMode])
                  //+ fvc::div(phiEvk_mode[speciei][vibMode])  // TODO create this in declaration/uppwind declarations
                  ==
                    thermalRelaxModelVT().QVT_mode(speciei, vibMode)
                  //+ thermalRelaxModelVV().QVV(speciei, vibMode)
                  + thermalRelaxModelHE().QHE()
                  //+ reaction->Scv(speciei) // TODO what to do here ?
                  - composition.pP()[Y.size()-1]*fvc::div(U)
                );
            }
            else
            {
                solve
                (
                    fvm::ddt(rhoEvk_mode[speciei][vibMode])
                  //+ fvc::div(phiEvk_mode[speciei][vibMode]) // TODO create this in declaration/uppwind declarations
                  ==
                    thermalRelaxModelVT().QVT_mode(speciei, vibMode)
                  //+ thermalRelaxModelVV().QVV(speciei, vibMode)
                  + thermalRelaxModelHE().QHE()
                  //+ reaction->Scv(speciei) // TODO what to do here ?
                );
            }
        }

        if (composition.particleType(speciei) == 2)
        {
            forAll(pD[speciei], celli)
            {
                if (Y[speciei][celli] > 1e-21)
                {
                    hevk_mode[speciei][vibMode][celli] = rhoEvk_mode[speciei][vibMode][celli]/pD[speciei][celli];
                }
                else
                {
                    hevk_mode[speciei][vibMode][celli] = 0.0;
                }
            }
        }

        hevk_mode[speciei][vibMode].correctBoundaryConditions(); //TODO correct bdry cdts
    }
}

thermo.correct();
thermo.correct2T();
transport().correct();

forAll (rhoEvk_mode, speciei)
{
    forAll (rhoEvk_mode[speciei], vibMode)
    {
        rhoEvk_mode[speciei][vibMode].boundaryFieldRef() = pD[speciei].boundaryField()*hevk_mode[speciei][vibMode].boundaryField();
    }
}

/*// --- Solve vibrational energy equations, VISCOUS
if (!inviscid)
{
    forAll (hevk, speciei)
    {
        if (composition.particleType(speciei) == 2)
        {
            solve
            (
                fvm::ddt(pD[speciei], hevk[speciei]) - fvc::ddt(pD[speciei], hevk[speciei])
              - fvm::laplacian(speciesDiffusion().rhoD(speciei), hevk[speciei], "laplacian(D,ev)")
              - fvm::laplacian(transport().alphave(speciei), hevk[speciei])
            );

        }

        rhoEvk[speciei] = pD[speciei]*hevk[speciei];
    }

    thermo.correct();
    thermo.correct2T();
    transport().correct();
}*/


// --- Calculation of the total vibrational energy
forAll (rhoEvk_mode, speciei)
{
    rhoEvk[speciei] = 0.0*rhoEvk_mode[0][0];
    forAll (rhoEvk_mode[speciei], vibMode)
    {
      rhoEvk[speciei] += rhoEvk_mode[speciei][vibMode];
      rhoEvk[speciei].boundaryFieldRef() += rhoEvk_mode[speciei][vibMode].boundaryField();
    }
}


rhoEv = rhoEvk[0]*0.0;
forAll (rhoEvk, speciei)
{
  rhoEv += rhoEvk[speciei];
  rhoEv.boundaryFieldRef() += rhoEvk[speciei].boundaryField();
}

ev = rhoEv/rho;
ev.boundaryFieldRef() = rhoEv.boundaryField()/rho.boundaryField();

thermo.correct();
thermo.correct2T();
transport().correct();
