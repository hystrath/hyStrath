/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Fick.H"
#include "fvm.H"
#include <ctime>

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::Fick<ThermoType>::updateCoefficients()
{     
    DijModel_().update();

    forAll(species(), speciei)
    {
        volScalarField tmpSum = 0 / Dij(0,0);

        forAll(species(), speciej)
        {
            if (speciej != speciei and thermo_.composition().particleType(speciej) != 0)
            {     
                tmpSum += thermo_.composition().X(speciej) / Dij(speciei, speciej);
            }
        }
        
        volScalarField& Xi = thermo_.composition().X(speciei);
        
        D_[speciei] = thermo_.rho()*(1.0 - Xi) 
            / (tmpSum + dimensionedScalar("VSMALL", dimTime/dimArea, Foam::VSMALL));   

        
        forAll(D_[speciei], celli)
        {
            if (1.0 - Xi[celli] < miniXs_)
            {
                D_[speciei][celli] = 0;
            }
        }
        
        forAll(D_[speciei].boundaryField(), patchi)
        {
            forAll(D_[speciei].boundaryField()[patchi], facei)
            {
                if (1.0 - Xi.boundaryField()[patchi][facei] < miniXs_)
                {
                    D_[speciei].boundaryField()[patchi][facei] = 0;
                }  
            }
        }
    }
} 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::Fick<ThermoType>::Fick
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    multiSpeciesTransportModel(thermo, turbulence),
    
    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),
    
    miniXs_(1.0e-4)
{    
    D_.setSize(species().size());
    
    forAll(species(), speciei)
    {
        D_.set
        (
            speciei, 
            new volScalarField
            (
                IOobject
                (
                    "rhoD_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("D", dimMass/dimLength/dimTime, 0.0)
            )
        );
    } 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::Fick<ThermoType>::correct()
{
    updateCoefficients();

    forAll(species(), speciei)
    {
        calculateJ(speciei);
    }
    
    calculateSumDiffusiveFluxes();
}

    
/*template<class ThermoType>
Foam::scalar Foam::Fick<ThermoType>::correct
(
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields
)
{
    updateCoefficients();

    scalar maxResidual = 0;
    scalar eqnResidual = 1;

    volScalarField yt = 0.0*thermo_.composition().Y(0);
    surfaceScalarField nt = turbulence_.phi();
    
    forAll(this->D_, i)
    {  
        volScalarField& yi = thermo_.composition().Y(i);
        surfaceScalarField& spMassFluxi = spMassFlux_[i];

        tmp<fv::convectionScheme<scalar> > mvConvection
        (
            fv::convectionScheme<scalar>::New
            (
                mesh_,
                fields,
                turbulence_.phi(),
                mesh_.divScheme("div(phi,Yi_h)")
            )
        );

        if (mesh_.relaxField("Yi"))//Mohsen
        {
            yi.storePrevIter();
        }
            
        tmp<fvScalarMatrix> yEqn
        (   
            fvm::ddt(thermo_.rho(), yi)
//           + fvm::div(turbulence_.phi(), yi, "div(phi,Yi_h)")
          + mvConvection->fvmDiv(turbulence_.phi(), yi)
          - fvm::laplacian(D_[i],yi, "laplacian(D,Yi)")
          ==
            Sy_[i]
        );

        eqnResidual = solve(yEqn() , mesh_.solver("Yi")).initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

        if (mesh_.relaxField("Yi"))//Mohsen
        {
	          yi.relax(mesh_.fieldRelaxationFactor("Yi"));//Mohsen
        }

        yi.max(0.0);
//         yi.min(1.0);

        spMassFluxi = yEqn().flux();

        spMassFluxt -= spMassFluxi;
        yt += yi;  
    }
     
    // Calculate inert species
    // CONSIDERED FOR DELETION VINCENT 29/01/2016
    volScalarField& yInert = thermo_.composition().Y()[inertIndex_];
    yInert = 1 - yt;
    forAll(yInert.boundaryField(), patchi)
    {
        forAll(yInert.boundaryField()[patchi], facei)
        {
            yInert.boundaryField()[patchi][facei] = 1 - yt.boundaryField()[patchi][facei];
        }
    }
    yInert.max(0.0);
    spMassFlux_[inertIndex_] = spMassFluxt;
    // END CONSIDERED FOR DELETION VINCENT 29/01/2016
          
    updateMolarFractions();

    return maxResidual;
}*/

template<class ThermoType>
bool Foam::Fick<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}
   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
