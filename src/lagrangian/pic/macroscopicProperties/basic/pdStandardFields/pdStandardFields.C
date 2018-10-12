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
    pdStandardFields

Description

\*----------------------------------------------------------------------------*/

#include "pdStandardFields.H"
#include "pdCloud.H"
#include "zeroGradientFvPatchFields.H"


// namespace Foam
// {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and cloud and write (for pdInitialise)
Foam::pdStandardFields::pdStandardFields
(
    const fvMesh& mesh,
    pdCloud& cloud,
    bool write
)
:
    mesh_(mesh),
    cloud_(cloud),
    AveragePtr_(
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                "Average",
                mesh_.time().timeName(),
                mesh_
            ),
            cloud.particleProperties(),
            mesh_
        )
    ),
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
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        ),
        zeroGradientFvPatchField<vector>::typeName
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
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    pdRhoN_
    (
        IOobject
        (
            "pdRhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
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
        ),
        zeroGradientFvPatchField<vector>::typeName
    )/*,
    stdU_
    (
        IOobject
        (
            "stdU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero),
        zeroGradientFvPatchField<vector>::typeName
    ),
    transT_
    (
        IOobject
        (
            "transT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    )*/
{
    volScalarField qCopy = q_;

    wordList qBFCopy = q_.boundaryField().types();

    Info << "Remember to change the type entries in q, rhoN, etc to 'calculated' for the '" << endl;

    forAll(mesh_.boundaryMesh(), i)
    {
        if (isA<polyPatch>(mesh_.boundaryMesh()[i]))
        {
            /*if(mesh_.boundaryMesh()[i].type() == "patch")
            {
                Info << "   " << mesh_.boundaryMesh()[i].name() << "' patch!" << endl;
            }*/

            if(mesh_.boundaryMesh()[i].type() == "wall")
            {
                Info << "   " << mesh_.boundaryMesh()[i].name() << "' wall!" << endl;
            }
        }
    }
}

// Construct from mesh and cloud (for pdFoam)
Foam::pdStandardFields::pdStandardFields
(
    const fvMesh& mesh,
    pdCloud& cloud
)
:
    mesh_(mesh),
    cloud_(cloud),
    AveragePtr_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                "Average",
                mesh_.time().timeName(),
                mesh_
            ),
            cloud.particleProperties(),
            mesh_
        )
    ),
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
    pdRhoN_
    (
        IOobject
        (
            "pdRhoN",
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
    )/*,
    stdU_
    (
        IOobject
        (
            "stdU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    transT_
    (
        IOobject
        (
            "transT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )*/
{

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdStandardFields::~pdStandardFields()
{
    //Info << "pdStandardFields Destructor" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pdStandardFields::calculateFields()
{
    scalarField& rhoN = rhoN_.primitiveFieldRef();

    scalarField& rhoM = rhoM_.primitiveFieldRef();

    scalarField& pdRhoN = pdRhoN_.primitiveFieldRef();

    scalarField& linearKE = linearKE_.primitiveFieldRef();

    scalarField& rotationalE = rotationalE_.primitiveFieldRef();

    scalarField& rotationalDof = rotationalDof_.primitiveFieldRef();

    scalarField& vibrationalE = vibrationalE_.primitiveFieldRef();

    scalarField& vibrationalDof = vibrationalDof_.primitiveFieldRef();

    //scalarField& transT = transT_.internalField();

    vectorField& momentum = momentum_.primitiveFieldRef();

    //vectorField& stdU = stdU_.internalField();

    forAllConstIter(pdCloud, cloud_, iter)
    {
        const pdParcel& p = iter();
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh_);
        //const label cellI = p.cell();

        //- apply averaging method to properties
        List<List<scalar> > shape = AveragePtr_->stencil(p.position(), tetIs);
        List<scalar> shapeID = shape[0];
        List<scalar> shapeWeight = shape[1];

        forAll(shapeID,cellI)
        {
            rhoN[shapeID[cellI]]            += 1.0 * shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            rhoM[shapeID[cellI]]            += cloud_.constProps(p.typeId()).mass()* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            pdRhoN[shapeID[cellI]]          += 1.0* shapeWeight[cellI];

            linearKE[shapeID[cellI]]        += 0.5*cloud_.constProps(p.typeId()).mass()*(p.U() & p.U())* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            rotationalE[shapeID[cellI]]     += p.ERot()* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            rotationalDof[shapeID[cellI]]   += cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom()* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            vibrationalE[shapeID[cellI]]    += p.EVib()* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            vibrationalDof[shapeID[cellI]]  += cloud_.constProps(p.typeId()).vibrationalDegreesOfFreedom()* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];

            momentum[shapeID[cellI]]        += cloud_.constProps(p.typeId()).mass()*p.U()* shapeWeight[cellI] /mesh_.cellVolumes()[shapeID[cellI]];
        }
    }


    rhoN *= cloud_.nParticle();
    //rhoN *= cloud_.nParticle()/mesh_.cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM *= cloud_.nParticle();
    //rhoM *= cloud_.nParticle()/mesh_.cellVolumes();
    rhoM_.correctBoundaryConditions();

    pdRhoN_.correctBoundaryConditions();

    linearKE *= cloud_.nParticle();
    //linearKE *= cloud_.nParticle()/mesh_.cellVolumes();
    linearKE_.correctBoundaryConditions();

    rotationalE *= cloud_.nParticle();
    //rotationalE *= cloud_.nParticle()/mesh_.cellVolumes();
    rotationalE_.correctBoundaryConditions();

    rotationalDof *= cloud_.nParticle();
    //rotationalDof *= cloud_.nParticle()/mesh_.cellVolumes();
    rotationalDof_.correctBoundaryConditions();

    vibrationalE *= cloud_.nParticle();
    //vibrationalE *= cloud_.nParticle()/mesh_.cellVolumes();
    vibrationalE_.correctBoundaryConditions();

    vibrationalDof *= cloud_.nParticle();
    //vibrationalDof *= cloud_.nParticle()/mesh_.cellVolumes();
    vibrationalDof_.correctBoundaryConditions();

    momentum *= cloud_.nParticle();
    //momentum *= cloud_.nParticle()/mesh_.cellVolumes();
    momentum_.correctBoundaryConditions();

    //- added calculation of stdU and transT for use in rhoEF
    /*
    forAll(transT,cI)
    {
        scalar V = mesh_.cellVolumes()[cI];

        if(rhoN[cI] > VSMALL)
        {
            stdU[cI] = momentum[cI]/ (rhoM[cI]);

            transT[cI] = 2.0/(3.0*physicoChemical::k.value()*rhoN[cI])
                        *(linearKE[cI] - 0.5*rhoM[cI]*(stdU[cI] & stdU[cI]));
        }
        else
        {
            stdU[cI]       = vector::zero;
            transT[cI]      = 0.0;
        }

    }
    stdU_.correctBoundaryConditions();
    transT_.correctBoundaryConditions();*/

}

void Foam::pdStandardFields::resetFields()
{
    q_              = dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0);

    rhoN_           = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0);

    rhoM_           =  dimensionedScalar("zero",  dimensionSet(1, -3, 0, 0, 0), VSMALL);

    pdRhoN_         = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), 0.0);

    linearKE_       = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), VSMALL);

    rotationalE_    = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    rotationalDof_  = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL);

    vibrationalE_   = dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0);

    vibrationalDof_ = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0), VSMALL);

    //transT_         = dimensionedScalar("0.0", dimTemperature, VSMALL);

    fD_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -1, -2, 0, 0),
        vector::zero
    );

    momentum_ = dimensionedVector
    (
        "zero",
        dimensionSet(1, -2, -1, 0, 0),
        vector::zero
    );

    /*stdU_ = dimensionedVector
    (
        "zero",
        dimensionSet(0, 1, -1, 0, 0),
        vector::zero
    );*/
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// }  // End namespace Foam

// ************************************************************************* //
