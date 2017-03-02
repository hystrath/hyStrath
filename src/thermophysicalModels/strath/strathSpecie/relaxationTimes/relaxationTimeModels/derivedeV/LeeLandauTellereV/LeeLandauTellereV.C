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

#include "LeeLandauTellereV.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::LeeLandauTellereV<ThermoType>::updateCoefficients()
{     
    taueViModel_().update();
} 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::LeeLandauTellereV<ThermoType>::LeeLandauTellereV
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    relaxationTimeModeleV(thermo, turbulence),
    
    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{    
    taueV_.setSize(solvedVibEqSpecies().size());
    
    forAll(taueV_, speciei)
    {
        taueV_.set
        (
            speciei, 
            new volScalarField
            (
                IOobject
                (
                    "taueV_" + solvedVibEqSpecies()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("taueV", dimTime, 0.0)
            )
        );
    } 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::LeeLandauTellereV<ThermoType>::correct()
{
    updateCoefficients(); 
    
    const volScalarField& ee = thermo_.composition().hevel("e-");
    const volScalarField& Xe = thermo_.composition().X("e-");
    const scalarField& eeCells = ee.internalField();
    
    forAll(solvedVibEqSpecies(), speciei)
    {  
        /*if(speciei != thermo_.composition().vibTempAssociativity("e-"))
        {
            const volScalarField& pD = thermo_.composition().pD(speciei);
            const volScalarField& ev = thermo_.composition().hevel(speciei);
            const volScalarField& taueV = this->taueV_[speciei];
            volScalarField& QeV = this->QeV_[speciei];
            
            const scalarField& pDCells = pD.internalField();
            const scalarField& evCells = ev.internalField();
            const scalarField& taueVCells = taueV.internalField();
            scalarField& QeVCells = QeV.internalField();
            
            forAll(QeVCells, celli)
            {        
                Info << "taueVCells[celli] " << tab << taueVCells[celli] << "eeCells[celli] " << eeCells[celli] << endl;
                QeVCells[celli] = pDCells[celli]/taueVCells[celli]*(evCells[celli]-eeCells[celli]);
             // Info << "QeVCells[celli]" << tab << QeVCells[celli] << endl;
            }
            
            forAll(QeV.boundaryField(), patchi)
            {        
                const fvPatchScalarField& pee = ee.boundaryField()[patchi];
                const fvPatchScalarField& ppD = pD.boundaryField()[patchi];
                const fvPatchScalarField& pev = ev.boundaryField()[patchi];
                const fvPatchScalarField& ptaueV = taueV.boundaryField()[patchi];
                
                fvPatchScalarField& pQeV = QeV.boundaryField()[patchi];
                
                forAll(pQeV, facei)
                {
                    pQeV[facei] = ppD[facei]/ptaueV[facei]*(pev[facei] - pee[facei]);
                }
            }
        } */
    } 
} 


template<class ThermoType>
bool Foam::LeeLandauTellereV<ThermoType>::read()
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
