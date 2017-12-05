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
    polyStandardFields

Description

\*----------------------------------------------------------------------------*/

#include "polyStandardFields.H"
#include "polyMoleculeCloud.H"
#include "zeroGradientFvPatchFields.H"


// namespace Foam
// {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and cloud and write (for mdInitialise)
Foam::polyStandardFields::polyStandardFields
(
    const fvMesh& mesh,
    polyMoleculeCloud& cloud,
    const word& fieldName,
    bool write
)
:
    mesh_(mesh),
    cloud_(cloud),
    fieldName_(fieldName),
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
    USAM_
    (
        IOobject
        (
            "U_SAM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 1, -1, 0, 0),
            vector::zero
        )
    )
    
    
{
//     volScalarField qCopy = q_;
    
//     wordList qBFCopy = q_.boundaryField().types();
   
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

// Construct from mesh and cloud (for mdFoam)
Foam::polyStandardFields::polyStandardFields
(
    const fvMesh& mesh,
    polyMoleculeCloud& cloud,
    const word& fieldName 
)
:
    mesh_(mesh),
    cloud_(cloud),
    fieldName_(fieldName),    
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
    USAM_
    (
        IOobject
        (
            "U_SAM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(0, 1, -1, 0, 0),
            vector::zero
        )
    )
{
    mass_.setSize();
    momentum_.setSize();
    
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyStandardFields::~polyStandardFields()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyStandardFields::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "volFieldsMethod_"+fieldName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("mass", mass_);
    dict.readIfPresent("momentum", momentum_);
    dict.readIfPresent("nTimeSteps", nTimeSteps_);
    
//     Info << "Some properties read in: "
//          << "mols = " << mols_[0] 
//          << ", mass = " << mass_[0]
//          << ", averagingCounter = " << averagingCounter_
//          << endl;
}

void polyStandardFields::writeOut()
{
    if (mesh_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFieldsMethod_"+fieldName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("mass", mass_);
        dict.add("momentum", momentum_);        
        dict.add("nTimeSteps", nTimeSteps_); 
        
        IOstream::streamFormat fmt = time_.time().writeFormat();
//         Pout << "fmt = " << fmt << endl;
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();
    
        dict.regIOobject::writeObject(fmt, ver, cmp);
        
//         Info<< "Some properties written out: "
//             << "mols = " << mols_[0]
//             << ", mass = " << mass_[0]
//             << ", averagingCounter = " << averagingCounter_
//             << endl;
    }
}

void Foam::polyStandardFields::calculateFields()
{
    scalarField& mass = mass_.internalField();
    
    scalarField& rhoN = rhoN_.internalField();

    scalarField& rhoM = rhoM_.internalField();

    vectorField& USAM = USAM_.internalField();

    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                label cellI mol().cell();
                
                const polyMolecule::constantProperties& constProp = molCloud_.constProps(mol().id());
                const scalar& massI = constProp.mass();
                
                mass[cellI] += massI;
                rhoN[cellI] += 1.0;                
                rhoM[cellI] += massI;
                
                const label& nMols = molCloud_.cellOccupancy[cellI];
                if(nMols > 0)
                {
                    USAM[cellI] += mol().v()/scalar(nMols);
                }
            }
        }
    }


    rhoN /= mesh_.cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM /= mesh_.cellVolumes();
    rhoM_.correctBoundaryConditions();

    USAM.correctBoundaryConditions();
}

void Foam::polyStandardFields::resetFields()
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
