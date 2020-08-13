//scalar timeg = runTime.elapsedCpuTime();
thermalRelaxModelVT().correct();
//Info<< "VTcorrect time" << tab << runTime.elapsedCpuTime() - timeg << " s" << endl;
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
forAll(rhoEvk, speciei)
{
    // NEW VINCENT 29/12/2016 *************************************************
    DynamicList<label> removedCellsFromMatrix;
    DynamicList<scalar> valuesToImpose;

    forAll(Y[speciei], celli)
    {
        if(Y[speciei][celli] < Foam::rho2ReactionThermo::miniYforSolvingEvEqn)
        {
            removedCellsFromMatrix.append(celli);
        }
    }

    forAll(removedCellsFromMatrix, i)
    {
        valuesToImpose.append(0.0);
    }
    // END NEW VINCENT 29/12/2016 *********************************************

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
            + reaction->Siir(speciei);

        if(composition.vibTempAssociativity("e-") == speciei)
        {
            rhoEvkInviscid += (fvc::grad(composition.pP("e-")) & U)
                - thermalRelaxModelHE().QHE();
        }
    }

    rhoEvkInviscid.setValues(removedCellsFromMatrix, valuesToImpose); // NEW VINCENT 29/12/2016

    rhoEvkInviscid.solve
    (
        mesh.solver("rhoEv_" + composition.species()[speciei])
    );

    forAll(Y[speciei], celli)
    {
        if(Y[speciei][celli] > Foam::rho2ReactionThermo::miniYforSolvingEvEqn)
        {
            hevk[speciei][celli] = rhoEvk[speciei][celli]/pD[speciei][celli];
        }
        /*else // DELETED VINCENT 09/02/2017
        {
            hevk[speciei][celli] = 0.0; //composition.Hvels(speciei, 1.0, 210.0);
        }*/
    }

    hevk[speciei].correctBoundaryConditions();
    rhoEvk[speciei].boundaryFieldRef() = pD[speciei].boundaryField()*hevk[speciei].boundaryField();
}

thermo.correctHEVK();

rhoEv *= 0.0;
forAll(hevk, speciei)
{
    rhoEv += pD[speciei]*hevk[speciei];
}

ev = rhoEv/rho;

