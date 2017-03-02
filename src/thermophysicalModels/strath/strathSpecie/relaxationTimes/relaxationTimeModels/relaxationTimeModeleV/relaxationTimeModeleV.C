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

#include "relaxationTimeModeleV.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(relaxationTimeModeleV, 0);
    defineRunTimeSelectionTable(relaxationTimeModeleV, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //  

  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relaxationTimeModeleV::relaxationTimeModeleV
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    IOdictionary
    (
        thermo.twoTemperatureDictionary()
    ),
    
    mesh_(thermo.Tt().mesh()), 
    thermo_(thermo),
    turbulence_(turbulence)
    
{  
    const word dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).substr
        (
            fileName(thermo.lookup("foamChemistryThermoFile")).find("constant/") + 9
        )
    );
    
    // Construct the relaxation time model
    taueViModel_.set
    (
        new eVModel
        (
            IOdictionary::name(),
            dictThermoPhy,
            species(), 
            thermo.composition().pP("e-"), 
            thermo.composition().Tv("e-")
        )
    ); 
    
    QeV_.setSize(solvedVibEqSpecies().size());
        
    forAll(solvedVibEqSpecies(), speciei)
    {
        QeV_.set
        (
            speciei, 
            new volScalarField
            (
                IOobject
                (
                    "QeV_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("QeV", dimensionSet(1,-1,-3,0,0), 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*Foam::tmp<Foam::volScalarField>
Foam::relaxationTimeModeleV::eVRelaxationSource()
{
    tmp<volScalarField> tQeV
    (
        new volScalarField
        (
            IOobject
            (
                "eVRelaxationSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("QeV", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    return tQeV;
}*/


bool Foam::relaxationTimeModeleV::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
