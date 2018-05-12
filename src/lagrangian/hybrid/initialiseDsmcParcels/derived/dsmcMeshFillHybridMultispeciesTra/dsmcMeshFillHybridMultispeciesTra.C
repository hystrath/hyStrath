/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcMeshFillHybridMultispeciesTra.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMeshFillHybridMultispeciesTra, 0);

addToRunTimeSelectionTable(dsmcConfiguration, dsmcMeshFillHybridMultispeciesTra,
    dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMeshFillHybridMultispeciesTra::dsmcMeshFillHybridMultispeciesTra
(
    dsmcCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    dsmcConfiguration(cloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMeshFillHybridMultispeciesTra::~dsmcMeshFillHybridMultispeciesTra()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dsmcMeshFillHybridMultispeciesTra::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const dictionary& speciesDict
    (
        dsmcInitialiseDict_.subDict("species")
    );

    List<word> molecules(speciesDict.toc());

// Files to be read:------------------------------------------------------------
    volVectorField UInitial
    (
        IOobject
        (
            "UCFD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );
    PtrList<volScalarField> TtInitial(molecules.size());
    PtrList<volScalarField> TrInitial(molecules.size());
    PtrList<volScalarField> TvInitial(molecules.size());
    PtrList<volScalarField> pInitial(molecules.size());
    PtrList<volScalarField> kappatrInitial(molecules.size());
    PtrList<volScalarField> muInitial(molecules.size());
//------------------------------------------------------------------------------

    PtrList<volScalarField> numberDensitiesField(molecules.size());
    PtrList<volVectorField> qt(molecules.size());
    PtrList<volTensorField> tau(molecules.size());

    dimensionedScalar massDim("massDim", dimMass, 1.0);
    scalar maxTranslationalTemperature = 0.0;
    List<scalar> numberDensities(molecules.size());
    numberDensities = 0.0;

    forAll(molecules, moleculeI)
    {
        dimensionedScalar thetaV("thetaV", dimTemperature,
            cloud_.constProps(moleculeI).thetaV());
        TtInitial.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "TtCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
        maxTranslationalTemperature = max(maxTranslationalTemperature,
            max(TtInitial[moleculeI].internalField()));

        TrInitial.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "TrCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        TvInitial.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "TvCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        pInitial.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "pCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        numberDensitiesField.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "numberDensity_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                pInitial[moleculeI] / (physicoChemical::k
                    * TtInitial[moleculeI])
            )
        );

        numberDensitiesField[moleculeI] /= cloud_.nParticle();

        kappatrInitial.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "kappatrCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        muInitial.set
        (
            moleculeI,
            volScalarField
            (
                IOobject
                (
                    "mutrCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        scalar rotDOF = cloud_.constProps(moleculeI).rotationalDegreesOfFreedom();

        qt.set
        (
            moleculeI,
            volVectorField
            (
                IOobject
                (
                    "qtCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -( (15.0 / (15.0 + 2.0 * rotDOF))
                    * kappatrInitial[moleculeI] / pInitial[moleculeI])
                    * fvc::grad(TtInitial[moleculeI]) * sqrt(2.0
                    * cloud_.constProps(moleculeI).mass() * massDim
                    / (physicoChemical::k * TtInitial[moleculeI]))
            )
        );
        qt[moleculeI].write();

        tau.set
        (
            moleculeI,
            volTensorField
            (
                IOobject
                (
                    "tauCFD_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                muInitial[moleculeI] * (fvc::grad(UInitial) +
                    dev2(Foam::T(fvc::grad(UInitial))))
                    / pInitial[moleculeI]
            )
        );
        tau[moleculeI].write();
    }


    forAll(molecules, moleculeI)
    {
        // Average number of particles per specie
        numberDensities[moleculeI] =
            numberDensitiesField[moleculeI].weightedAverage(mesh_.V()).value();
    }

    forAll(mesh_.cells(), cellI)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            cellI
        );

        forAll(cellTets, tetI)
        {
            const tetIndices& cellTetIs = cellTets[tetI];

            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);

                label typeId(findIndex(cloud_.typeIdList(), moleculeName));

                if (typeId == -1)
                {
                    FatalErrorIn("Foam::dsmcCloud<dsmcParcel>::initialise")
                        << "typeId " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);

                scalar numberDensity = numberDensitiesField[i][cellI];
                scalar translationalTemperature = TtInitial[i][cellI];
                scalar rotationalTemperature = TrInitial[i][cellI];
                scalar vibrationalTemperature = TvInitial[i][cellI];
                vector velocity = UInitial[cellI];

                // Calculate the number of particles required
                scalar particlesRequired = numberDensity * tetVolume;

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.sample01<scalar>()
                )
                {
                    nParticlesToInsert++;
                }

                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = cloud_.chapmanEnskogVelocity
                    (
                        translationalTemperature,
                        cP.mass(),
                        qt[i][cellI],
                        tau[i][cellI]
                    );

                    scalar ERot = cloud_.equipartitionRotationalEnergy
                    (
                        rotationalTemperature,
                        cP.rotationalDegreesOfFreedom()
                    );
        
                    scalar EVib = cloud_.equipartitionVibrationalEnergy
                    (
                        vibrationalTemperature,
                        cP.vibrationalDegreesOfFreedom(),
                        i
                    );

                    U += velocity;
                    
                    label newParcel = -1;
                    
                    label classification = 0;

                    const scalar& RWF = cloud_.coordSystem().recalculateRWF(cellI);

                    cloud_.addNewParcel
                    (
                        p,
                        U,
                        RWF,
                        ERot,
                        ELevel,
                        cellI,
                        cellTetIs.face(),
                        cellTetIs.tetPt(),
                        typeId,
                        newParcel,
                        classification,
                        vibLevel
                    );
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const dsmcParcel::constantProperties& cP = cloud_.constProps
    (
        mostAbundantType
    );

    cloud_.sigmaTcRMax().primitiveFieldRef() = cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
    (
        maxTranslationalTemperature,
        cP.mass()
    );

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


} // End namespace Foam

// ************************************************************************* //
