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
    pdEmFields

Description

\*----------------------------------------------------------------------------*/

#include "pdEmFields.H"
#include "pdCloud.H"
//#include "zeroGradientFvPatchFields.H"

// namespace Foam
// {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and cloud and write (for pdInitialise)
Foam::pdEmFields::pdEmFields
(
    const fvMesh& mesh,
    pdCloud& cloud,
    bool write
)
:
    mesh_(mesh),
    cloud_(cloud),
    EAveragePtr_(NULL),
    BAveragePtr_(NULL),
    phiEAveragePtr_(NULL),
    //JpAveragePtr_(NULL),
    rhoAveragePtr_(NULL),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 1, -3, 0, 0, -1, 0),
            vector::zero
        )
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -2, 0, 0, -1, 0),
            vector::zero
        )
    ),
    phiE_
    (
        IOobject
        (
            "phiE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.0)
    ),
    Jp_
    (
        IOobject
        (
            "Jp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 0, -2, 0, 0, 1, 0),
            vector::zero
        ),
        zeroGradientFvPatchField<vector>::typeName
    ),
    rhoQ_
    (
        IOobject
        (
            "rhoQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0, 1, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    wallQ_
    (
        IOobject
        (
            "wallQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, 1, 0, 0, 1, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    rhoEF_
    (
        IOobject
        (
            "rhoEF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    /*lambdaD_
    (
        IOobject
        (
            "lambdaD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),*/
    debugCharge_(false),
    debugField_(false)
{

    volScalarField rhoQCopy = rhoQ_;

    wordList rhoQBFCopy = rhoQ_.boundaryField().types();

    Info << "Remember to change the type entries in phiE, E, B etc to the correct BC for:" << endl;

    forAll(mesh_.boundaryMesh(), i)
    {
        const polyPatch& patch = mesh.boundaryMesh()[i];

        if (isA<polyPatch>(patch))
        {
            if(mesh_.boundaryMesh()[i].type() == "patch")
            {
                Info << "   '" << mesh_.boundaryMesh()[i].name() << "' patch!" << endl;
            }

            if(mesh_.boundaryMesh()[i].type() == "wall")
            {
                Info << "   '" << mesh_.boundaryMesh()[i].name() << "' wall!" << endl;
            }
        }
    }

    Info << "Remember to change the type entries in Jp to calculated for:" << endl;

    forAll(mesh_.boundaryMesh(), i)
    {
        const polyPatch& patch = mesh.boundaryMesh()[i];

        if (isA<polyPatch>(patch))
        {
            if(mesh_.boundaryMesh()[i].type() == "wall")
            {
                Info << "   '" << mesh_.boundaryMesh()[i].name() << "' wall!" << endl;
            }
        }
    }
}

// Construct from mesh and cloud (for pdFoam)
Foam::pdEmFields::pdEmFields
(
    const fvMesh& mesh,
    pdCloud& cloud
)
:
    mesh_(mesh),
    cloud_(cloud),
    EAveragePtr_
    (
        AveragingMethod<vector>::New
        (
            IOobject
            (
                "EAverage",
                mesh_.time().timeName(),
                mesh_
            ),
            cloud_.particleProperties(),
            mesh_
        )
    ),
    BAveragePtr_
    (
        AveragingMethod<vector>::New
        (
            IOobject
            (
                "BAverage",
                mesh_.time().timeName(),
                mesh_
            ),
            cloud_.particleProperties(),
            mesh_
        )
    ),
    phiEAveragePtr_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                "phiEAverage",
                mesh_.time().timeName(),
                mesh_
            ),
            cloud_.particleProperties(),
            mesh_
        )
    ),
//    JpAveragePtr_
//    (
//        AveragingMethod<vector>::New
//        (
//            IOobject
//            (
//                "JpAverage",
//                mesh_.time().timeName(),
//                mesh_
//            ),
//            cloud_.particleProperties(),
//            mesh_
//        )
//    ),
    rhoAveragePtr_
    (
        AveragingMethod<scalar>::New
        (
            IOobject
            (
                "rhoAverage",
                mesh_.time().timeName(),
                mesh_
            ),
            cloud_.particleProperties(),
            mesh_
        )
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phiE_
    (
        IOobject
        (
            "phiE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Jp_
    (
        IOobject
        (
            "Jp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoQ_
    (
        IOobject
        (
            "rhoQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    wallQ_
    (
        IOobject
        (
            "wallQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rhoEF_
    (
        IOobject
        (
            "rhoEF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    /*lambdaD_
    (
        IOobject
        (
            "lambdaD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),*/
    debugCharge_(false),
    debugField_(false)
{
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pdEmFields::~pdEmFields()
{
    //Info << "pdEmFields Destructor" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//calculates fields based on particle distribution
void Foam::pdEmFields::calculateFields()
{
    /**scatter particle charge to mesh**/
    const scalar e = electromagnetic::e.value();   //- elementary charge

    //loop through particles and add charge to average mesh
    forAllConstIter(pdCloud, cloud_, part)
    {
        const pdParcel& p = part();
        const pdParcel::constantProperties& cP = cloud_.constProps(p.typeId());
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh_);

        const scalar q  = cP.Ze()*e;

        const vector U  =  p.U();

        //- get stencil to interpolate particle to mesh
        List<List<scalar> > shape = rhoAveragePtr_->stencil(p.position(), tetIs);
        List<scalar> shapeID = shape[0];
        List<scalar> shapeWeight = shape[1];

        forAll(shapeID,cellI)
        {
            rhoQ_[shapeID[cellI]]   += q * shapeWeight[cellI] /mesh_.V()[shapeID[cellI]];
            Jp_[shapeID[cellI]]     += q * U * shapeWeight[cellI] /mesh_.V()[shapeID[cellI]];
        }
    }

    rhoQ_   *= cloud_.nParticle();
    rhoQ_.correctBoundaryConditions();

    Jp_     *= cloud_.nParticle();
    Jp_.correctBoundaryConditions();

    if(debugCharge_)
    {
        scalar counter = 0;
        forAll(rhoQ_,cI)
        {
            counter += rhoQ_[cI]*mesh_.V()[cI];
        }
        Info << "Charge on mesh: " << counter << endl;

        counter = 0;
        rhoAveragePtr_->updateField(rhoQ_);
        forAllConstIter(pdCloud, cloud_, part)
        {
            const pdParcel& p = part();
            const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh_);

            counter += rhoAveragePtr_->interpolate(p.position(), tetIs)*mesh_.V()[p.cell()];
        }
        Info << "Charge reconstructed: " << counter << endl;
    }

    /** calculate potential **/
    cloud_.electronModel().calculatePotential();

    /** calculate field **/
    E_ = -fvc::grad(phiE_);
    E_.correctBoundaryConditions();

    //- add external electric field

    //- add external magnetic field

    /** setup fields for interpolation **/

    if(debugField_)
    {
        E_ = dimensionedVector("one",  dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::one);
    }

    EAveragePtr_->updateField(E_);
    phiEAveragePtr_->updateField(phiE_);
    //BAveragePtr_->updateField(B_);

}

//calculates forces on particles due to electromagnetic fields
void Foam::pdEmFields::calculateForces()
{
    //- intiailise constant variables
    const scalar e = electromagnetic::e.value();
    const scalar nParticles = cloud_.nParticle();

    //- loop through all particles and modify particle acceleration
    scalar counter_part = 0.0;

    forAllIter(pdCloud, cloud_, part)
    {
        pdParcel& p = part();                   //- reference to particle
        vector&   a = p.A();                    //- reference to particle acceleration

        //const vector& u = p.U();              //- reference to particle velocity

        const pdParcel::constantProperties& cP = cloud_.constProps(p.typeId());
        const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh_);

        const scalar q = nParticles*cP.Ze()*e; //-nParticle included for clarity, cancels out
        const scalar m = nParticles*cP.mass(); //-nParticle included for clarity, cancels out

        /** gather electric field from mesh and interpolate to particle **/
        const vector E = EAveragePtr_->interpolate(p.position(), tetIs);
        //const vector B = BAveragePtr_->interpolate(p.position(), tetIs);

        /** lorentz force calculation **/
        //a = q/m * (E + (U_ ^ B));
        a = q/m*E;

        if(debugField_)
        {
            counter_part++;
            Info << "p: " << counter_part << " E: " << E << endl;
        }

        /** record particle potential energy **/
        p.EPot() = phiEAveragePtr_->interpolate(p.position(), tetIs);
    }
}

//rolls back initialised particles to k-0.5
void Foam::pdEmFields::setupLeapfrog()
{

    forAllIter(pdCloud, cloud_, part)
    {
        pdParcel& p = part();

        if(p.newParcel() == 1)
        {
            vector& U_ = p.U();
            label& newParcel_ = p.newParcel();

            U_ -= p.A()*mesh_.time().deltaTValue();
            newParcel_ = 0;
        }
    }
}

//resets particle fields
void Foam::pdEmFields::resetFields()
{
    /** reset local fields **/
    //- electric field
    E_      = dimensionedVector("zero",  dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero);
    //- magnetic field
    //B_      = dimensionedVector("zero",  dimensionSet(1, 0, -2, 0, 0, -1, 0), vector::zero);
    //- total charge density (C/m^3)
    rhoQ_   = dimensionedScalar("zero",  dimensionSet(0, -3, 1, 0, 0, 1, 0), 0.0);
    //- electron cloud density (1/m^3)
    rhoEF_   = dimensionedScalar("zero",  dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0);
    //- debye length
    //lambdaD_   = dimensionedScalar("zero",  dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0);
    //- ion current density
    Jp_ = dimensionedVector("zero", dimensionSet(0, 0, -2, 0, 0, 1, 0), vector::zero);

    /** reset "cloud" fields **/
    //- cloud electric field
    EAveragePtr_->resetFields();
    //- cloud magnetic field
    BAveragePtr_->resetFields();
    //- cloud charge density
    rhoAveragePtr_->resetFields();
    //- cloud potential field
    phiEAveragePtr_->resetFields();
    //- cloud current density
    //JpAveragePtr_->resetFields();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// }  // End namespace Foam

// ************************************************************************* //
