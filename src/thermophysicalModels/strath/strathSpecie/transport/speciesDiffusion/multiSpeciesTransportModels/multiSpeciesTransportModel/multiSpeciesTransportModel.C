/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "multiSpeciesTransportModel.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(multiSpeciesTransportModel, 0);
    defineRunTimeSelectionTable(multiSpeciesTransportModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSpeciesTransportModel::multiSpeciesTransportModel
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    IOdictionary
    (
        thermo.transportDictionary()
    ),
    
    mesh_(thermo.Tt().mesh()),     
    thermo_(thermo),
    turbulence_(turbulence),
    
    spMassFlux_(species().size()),
    JnonCorrected_(species().size()),
    
    sumDiffusiveFluxes_
    (
        IOobject
        (
            "sumDiffusiveFluxes",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("sumDiffusiveFluxes", dimMass/dimArea/dimTime, vector::zero)
    ),
    
    useNonCorrected_(subDict("transportModels").lookupOrDefault<bool>("nonCorrectedDiffusiveFluxes", false))
{  
    const word dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).substr
        (
            fileName(thermo.lookup("foamChemistryThermoFile")).find("constant/") + 9
        )
    );
        
    DijModel_.set
    (
        new diffusivityModel
        (
            IOdictionary::name(),
            dictThermoPhy,
            thermo.p(), 
            thermo.Tt(), 
            species()
         )
    );
    
    forAll(species(), speciei)
    {
        spMassFlux_.set
        (
            speciei,
            new surfaceScalarField
            (
                IOobject
                (
                    "massFlux_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("massFlux", dimMass/dimTime, 0.0)
            )
        );
        
        JnonCorrected_.set
        (
            speciei,
            new volVectorField
            (
                IOobject
                (
                    "JnonCorrected_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("JnonCorrected", dimMass/dimArea/dimTime, vector::zero)
            )
        );
    } 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiSpeciesTransportModel::calculateJ
(
    const label i
)
{
    if(thermo_.composition().particleType(i) != 0)
    {
        JnonCorrected_[i] = -rhoD(i)*fvc::grad(thermo_.composition().Y(i));
    }
    else
    {
        volVectorField sum = JnonCorrected_[0]*thermo_.composition().particleCharge(0)/W(0);
        for(label specier=1 ; specier < species().size(); specier++)
        {
            if(thermo_.composition().particleType(specier) != 0)
            {
                sum += JnonCorrected_[specier]*thermo_.composition().particleCharge(specier)/W(specier);
            }
        }
        JnonCorrected_[i] = W(i)*sum;
    }
}


void Foam::multiSpeciesTransportModel::calculateSumDiffusiveFluxes()
{
    // Uses the non-corrected diffusive fluxes for the calculation
    // of the diffusive fluxes
    sumDiffusiveFluxes_ = JnonCorrected_[0];
    
    for(label speciej=1 ; speciej < species().size(); speciej++)
    {
        if(thermo_.composition().particleType(speciej) != 0)
        {
            sumDiffusiveFluxes_ += JnonCorrected_[speciej];
        }
    }
}


Foam::volVectorField
Foam::multiSpeciesTransportModel::Jcorrected(const label i) const
{
    if(thermo_.composition().particleType(i) != 0 and not useNonCorrected_)
    {
        return J(i) - thermo_.composition().Y(i)*sumDiffusiveFluxes();
    }
    else
    {
        return J(i);
    }
}


Foam::volVectorField
Foam::multiSpeciesTransportModel::multiSpeciesHeatSource() const
{
    volVectorField multiSpeciesHeatSource
    (
        IOobject
        (
            "multiSpeciesHeatSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("multiSpeciesHeatSource", dimEnergy/dimArea/dimTime, vector::zero)
    );
    
    forAll(species(), speciej)
    {
        if(thermo_.composition().particleType(speciej) != 0)
        {
            const volScalarField pCells = thermo_.p();
            const volScalarField TtCells = thermo_.Tt();
            const volScalarField TvCells = thermo_.composition().Tv(speciej);
            
            // Initialisation of the volScalarField with the right units
            volScalarField hsj = thermo_.composition().e();
          
            forAll(hsj, celli)
            {
                hsj[celli] = hs(speciej, pCells[celli], TtCells[celli], TvCells[celli]);
            }
            
            forAll(hsj.boundaryField(), patchi)
            {
                const fvPatchScalarField& pp = thermo_.p().boundaryField()[patchi];
                const fvPatchScalarField& pTt = thermo_.Tt().boundaryField()[patchi];
                const fvPatchScalarField& pTv = thermo_.composition().Tv(speciej).boundaryField()[patchi];

                fvPatchScalarField& phsj = hsj.boundaryField()[patchi];

                forAll(pTt, facei)
                {
                    phsj[facei] = hs(speciej, pp[facei], pTt[facei], pTv[facei]);
                }
            }
            
            if(useNonCorrected_)
            {
                multiSpeciesHeatSource += JnonCorrected_[speciej]*hsj;
            }
            else
            {
                multiSpeciesHeatSource += Jcorrected(speciej)*hsj;
            }
        }
    }
    
    return multiSpeciesHeatSource;
}


void Foam::multiSpeciesTransportModel::getSpeciesMassFlux
(
    const label i, 
    const surfaceScalarField& flux
)
{
    spMassFlux_[i] = flux;
}


bool Foam::multiSpeciesTransportModel::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
