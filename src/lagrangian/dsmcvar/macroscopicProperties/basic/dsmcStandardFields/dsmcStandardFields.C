/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    dsmcStandardFields

Description

\*----------------------------------------------------------------------------*/

#include "dsmcStandardFields.H"
#include "dsmcCloud.H"
#include "zeroGradientFvPatchFields.H"


// namespace Foam
// {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and cloud and write (for dsmcInitialise)
Foam::dsmcStandardFields::dsmcStandardFields
(
    const fvMesh& mesh,
    dsmcCloud& cloud,
    bool write
)
:
    mesh_(mesh),
    cloud_(cloud),
    q_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    fD_
    (
        IOobject
        (
            "fD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), 0.0)
    ),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    linearKE_
    (
        IOobject
        (
            "linearKE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    rotationalE_
    (
        IOobject
        (
            "rotationalE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    rotationalDof_
    (
        IOobject
        (
            "rotationalDof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    vibrationalE_
    (
        IOobject
        (
            "vibrationalE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    vibrationalDof_
    (
        IOobject
        (
            "vibrationalDof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0)
    ),
    momentum_
    (
        IOobject
        (
            "momentum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -1, 0, 0),
            vector::zero
        )
    )
{
    volScalarField qCopy = q_;
    
    wordList qBFCopy = q_.boundaryField().types();
   
    forAll(mesh_.boundaryMesh(), i)
    {        
        if (isA<polyPatch>(mesh_.boundaryMesh()[i]))
        {
            if(mesh_.boundaryMesh()[i].type() == "patch")
            {  
                Info << "Remember to change the type entries in q, rhoN, etc to 'zeroGradient' for the '" <<
                mesh_.boundaryMesh()[i].name() << "' patch!" << endl << endl;
            }
        }
    }
}

// Construct from mesh and cloud (for dsmcFoam)
Foam::dsmcStandardFields::dsmcStandardFields
(
    const fvMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(mesh),
    cloud_(cloud),
    q_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    fD_
    (
        IOobject
        (
            "fD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    linearKE_
    (
        IOobject
        (
            "linearKE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rotationalE_
    (
        IOobject
        (
            "rotationalE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rotationalDof_
    (
        IOobject
        (
            "rotationalDof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    vibrationalE_
    (
        IOobject
        (
            "vibrationalE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    vibrationalDof_
    (
        IOobject
        (
            "vibrationalDof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    momentum_
    (
        IOobject
        (
            "momentum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dsmcStandardFields::~dsmcStandardFields()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcStandardFields::calculateFields()
{
    scalarField& rhoN = rhoN_.internalField();

    scalarField& rhoM = rhoM_.internalField();

    scalarField& dsmcRhoN = dsmcRhoN_.internalField();

    scalarField& linearKE = linearKE_.internalField();

    scalarField& rotationalE = rotationalE_.internalField();

    scalarField& rotationalDof = rotationalDof_.internalField();
    
    scalarField& vibrationalE = vibrationalE_.internalField();

    scalarField& vibrationalDof = vibrationalDof_.internalField();

    vectorField& momentum = momentum_.internalField();

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& p = iter();
        const label cellI = p.cell();

        rhoN[cellI]++;

        rhoM[cellI] += cloud_.constProps(p.typeId()).mass();

        dsmcRhoN[cellI]++;

        linearKE[cellI] += 0.5*cloud_.constProps(p.typeId()).mass()*(p.U() & p.U());

        rotationalE[cellI] += p.ERot();

        rotationalDof[cellI] += cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();

        vibrationalE[cellI] += p.vibLevel()*cloud_.constProps(p.typeId()).thetaV()*physicoChemical::k.value();

        vibrationalDof[cellI] += cloud_.constProps(p.typeId()).vibrationalDegreesOfFreedom();

        momentum[cellI] += cloud_.constProps(p.typeId()).mass()*p.U();
    }


    rhoN *= cloud_.nParticle()/mesh_.cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM *= cloud_.nParticle()/mesh_.cellVolumes();
    rhoM_.correctBoundaryConditions();

    dsmcRhoN_.correctBoundaryConditions();

    linearKE *= cloud_.nParticle()/mesh_.cellVolumes();
    linearKE_.correctBoundaryConditions();

    rotationalE *= cloud_.nParticle()/mesh_.cellVolumes();
    rotationalE_.correctBoundaryConditions();

    rotationalDof *= cloud_.nParticle()/mesh_.cellVolumes();
    rotationalDof_.correctBoundaryConditions();
    
    vibrationalE *= cloud_.nParticle()/mesh_.cellVolumes();
    vibrationalE_.correctBoundaryConditions();

    vibrationalDof *= cloud_.nParticle()/mesh_.cellVolumes();
    vibrationalDof_.correctBoundaryConditions();

    momentum *= cloud_.nParticle()/mesh_.cellVolumes();
    momentum_.correctBoundaryConditions();
}

void Foam::dsmcStandardFields::resetFields()
{
    q_ = dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0);

    fD_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -1, -2, 0, 0),
        vector::zero
    );

    rhoN_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL );

    rhoM_ =  dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), VSMALL);

    dsmcRhoN_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0);

    linearKE_ = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    rotationalE_ = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    rotationalDof_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL);
    
    vibrationalE_ = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    vibrationalDof_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL);

    momentum_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -2, -1, 0, 0),
        vector::zero
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// }  // End namespace Foam

// ************************************************************************* //
