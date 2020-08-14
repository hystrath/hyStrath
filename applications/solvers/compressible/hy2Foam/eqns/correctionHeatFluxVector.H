
if(hyDev)
{
    /*
    *heatFluxVector = -turbulence->kappaEff()*fvc::grad(thermo.T());

    scalarField Ycoord = mesh.C().component(1);
    vectorField& myHFV = heatFluxVector->internalField();

    forAll(Y, speciei)
    {
        *heatFluxVector -= transport().kappave(speciei)*fvc::grad(composition.Tv()[speciei]);
    }

    volScalarField hFV = heatFluxVector->component(1);
    volVectorField derivhFV = fvc::grad(hFV);
    scalarField derivhFV1 = derivhFV.component(1);

    forAll(*heatFluxVector, celli)
    {
        scalar xs = 1.0;
        scalar mfp = 0.015;
        scalar extentLayer = 5.0*mfp;
        scalar x1 = xs - extentLayer - 0.0025;
        scalar y1 = hFV[897];; //3.07989; //hFV[897];
        scalar a = derivhFV1[205]; //(y1-y2)/(x1-x2);
        scalar gprime = 0; //max(0, min(a/1.05, a));
        scalar b = y1 - a*x1;
        scalar c = (a-gprime)/(2.0*(x1-xs));
        scalar d = -2.0*c*xs + gprime;
        scalar e = -c*Foam::sqr(x1) + (a-d)*x1 + b;

        if(Ycoord[celli] > x1)
        {
            myHFV[celli].y() = c*Foam::sqr(Ycoord[celli]) + d*Ycoord[celli] + e;
        }


        xs = 0.0;
        extentLayer = 0.1;
        x1 = 1-0.9625;
        y1 = hFV[38]; //-3.07989; //hFV[102];
        b = y1 - a*x1;
        c = (a-gprime)/(2.0*(x1-xs));
        d = -2.0*c*xs + gprime;
        e = -c*Foam::pow(x1, 2.0) + (a-d)*x1 + b;

        if(Ycoord[celli] < x1)
        {
            myHFV[celli].y() = c*Foam::sqr(Ycoord[celli]) + d*Ycoord[celli] + e;
        }
    }

    *heatFluxVector += turbulence->kappaEff()*fvc::grad(thermo.T());
    forAll(Y, speciei)
    {
        *heatFluxVector += transport().kappave(speciei)*fvc::grad(composition.Tv()[speciei]);
    }
    */
}

