/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "noMHD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(noMHD, 0);
        addToMhdRunTimeSelectionTables(noMHD);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noMHD::noMHD(const rho2ReactionThermo& thermo)
:
    mhdModel(thermo)
{}


noMHD::noMHD
(
    const dictionary& dict,
    const rho2ReactionThermo& thermo
)
:
    mhdModel(thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noMHD::~noMHD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool noMHD::read()
{
    return mhdModel::read();
}


void noMHD::update(const volVectorField& U)
{}


tmp<volScalarField> noMHD::jouleHeating(const volVectorField& U) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "jouleHeating",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "jouleHeating",
                dimensionSet(1, -1, -3, 0, 0, 0, 0),
                0.0
            )
        )
    );
}


tmp<volVectorField> noMHD::lorentzForce() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "lorentzForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "lorentzForce",
                dimensionSet(1, -2, -2, 0, 0, 0, 0),
                vector::zero
            )
        )
    );
}


volVectorField& noMHD::E()
{
    tmp<volVectorField> tE
    (
        new volVectorField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "E",
                dimensionSet(1, 1, -3, 0, 0, -1, 0),
                vector::zero
            )
        )
    );
    
    volVectorField& E = tE.ref();
    
    return E;
}


const volVectorField& noMHD::E() const
{
    tmp<volVectorField> tE
    (
        new volVectorField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "E",
                dimensionSet(1, 1, -3, 0, 0, -1, 0),
                vector::zero
            )
        )
    );
    
    const volVectorField& E = tE;
    
    return E;
}


volScalarField& noMHD::elecPot()
{
    tmp<volScalarField> telecPot
    (
        new volScalarField
        (
            IOobject
            (
                "elecPot",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "elecPot",
                dimensionSet(1, 2, -3, 0, 0, -1, 0),
                0.0
            )
        )
    );
    
    volScalarField& elecPot = telecPot.ref();
    
    return elecPot;
}


const volScalarField& noMHD::elecPot() const
{
    tmp<volScalarField> telecPot
    (
        new volScalarField
        (
            IOobject
            (
                "elecPot",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "elecPot",
                dimensionSet(1, 2, -3, 0, 0, -1, 0),
                0.0
            )
        )
    );
    
    const volScalarField& elecPot = telecPot;
    
    return elecPot;
}


volVectorField& noMHD::B()
{
    tmp<volVectorField> tB
    (
        new volVectorField
        (
            IOobject
            (
                "B",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "B",
                dimensionSet(1, 0, -2, 0, 0, -1, 0),
                vector::zero
            )
        )
    );
    
    volVectorField& B = tB.ref();
    
    return B;
}


const volVectorField& noMHD::B() const
{
    tmp<volVectorField> tB
    (
        new volVectorField
        (
            IOobject
            (
                "B",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "B",
                dimensionSet(1, 0, -2, 0, 0, -1, 0),
                vector::zero
            )
        )
    );
    
    const volVectorField& B = tB;
    
    return B;
}


const volTensorField& noMHD::sigma() const
{
    tmp<volTensorField> tsigma
    (
        new volTensorField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor
            (
                "sigma",
                dimensionSet(-1, -3, 3, 0, 0, 2, 0),
                tensor::zero
            )
        )
    );
    
    const volTensorField& sigma = tsigma;
    
    return sigma;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace FOAM

// ************************************************************************* //
