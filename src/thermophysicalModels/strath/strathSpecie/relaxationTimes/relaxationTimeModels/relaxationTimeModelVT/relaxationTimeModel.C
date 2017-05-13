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

#include "relaxationTimeModel.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(relaxationTimeModel, 0);
    defineRunTimeSelectionTable(relaxationTimeModel, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //  

  
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relaxationTimeModel::relaxationTimeModel
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
    const word dict2T(IOdictionary::name()), dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).substr
        (
            fileName(thermo.lookup("foamChemistryThermoFile")).find("constant/") + 9
        )
    );
    
    // Construct the relaxation time model
    tauVTijModel_.set
    (
        new VTModel
        (
            dict2T,
            dictThermoPhy,
            solvedVibEqSpecies(),
            species(), 
            thermo.p(), 
            thermo.Tt(), 
            thermo.composition().Tv(), 
            thermo.composition().nD()
        )
    );    
    
    QVT_.setSize(solvedVibEqSpecies().size()); //NEW VINCENT 05/08/2016
    //QVTmode_.setSize(species().size()); // TODO ONGOING WORK

    forAll(solvedVibEqSpecies(), speciei) //NEW VINCENT 05/08/2016
    {
        QVT_.set
        (
            speciei, 
            new volScalarField
            (
                IOobject
                (
                    "QVT_" + solvedVibEqSpecies()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("QVT", dimensionSet(1,-1,-3,0,0), 0.0)
            )
        );
    }
    
    /*forAll(QVTmode_, speciei) // TODO ONGOING WORK
    {
        QVTmode_.set
        (
            speciei,
            new PtrList<volScalarField>(thermo.composition().noVibrationalTemp(speciei))
        );
    }
    
    forAll(QVTmode_, speciei)
    {
      forAll(QVTmode_[speciei], vibMode)
      {
        QVTmode_[speciei].set
        (
            vibMode, 
            new volScalarField
            (
                IOobject
                (
                    "QVT_" + species()[speciei] + "." + word(vibMode+1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("QVT", dimensionSet(1,-1,-3,0,0), 0.0)
            )
        );
      }
    }*/
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*Foam::tmp<Foam::volScalarField>
Foam::relaxationTimeModel::VTRelaxationSource()
{
    tmp<volScalarField> tQVT
    (
        new volScalarField
        (
            IOobject
            (
                "VTRelaxationSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("QVT", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    return tQVT;
}*/


bool Foam::relaxationTimeModel::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
