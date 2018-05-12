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

#include "dsmcMeshFillFromCFDHybrid.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMeshFillFromCFDHybrid, 0);

addToRunTimeSelectionTable(dsmcConfiguration, dsmcMeshFillFromCFDHybrid,
    dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMeshFillFromCFDHybrid::dsmcMeshFillFromCFDHybrid
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

dsmcMeshFillFromCFDHybrid::~dsmcMeshFillFromCFDHybrid()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dsmcMeshFillFromCFDHybrid::setInitialConfiguration()
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
            "../meshCFD/U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );
    volScalarField TtInitial
    (
        IOobject
        (
            "../meshCFD/Tt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );
    volScalarField pTotInitial
    (
        IOobject
        (
            "../meshCFD/p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    /*volScalarField rhoTotInitial
    (
        IOobject
        (
            "../meshCFD/rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_
    );*/

    PtrList<volScalarField> TvInitial(molecules.size());
    PtrList<volScalarField> nDInitial(molecules.size());
    PtrList<volScalarField> kappatrInitial(molecules.size());
    PtrList<volScalarField> kappaveInitial(molecules.size());
    PtrList<volScalarField> muInitial(molecules.size());
    PtrList<volVectorField> diffFluxInitial(molecules.size());
//------------------------------------------------------------------------------

    PtrList<volScalarField> pInitial(molecules.size());
    PtrList<volVectorField> qt(molecules.size());
    PtrList<volVectorField> qr(molecules.size());
    PtrList<volVectorField> qv(molecules.size());
    PtrList<volVectorField> D(molecules.size());
    PtrList<volTensorField> tau(molecules.size());

    scalar maxTranslationalTemperature = gMax(TtInitial);
    List<scalar> numberDensities(molecules.size());
    numberDensities = 0.0;

    dimensionedScalar k = physicoChemical::k;

    forAll(molecules, moleculeI)
    {
        scalar thetaV = cloud_.constProps(moleculeI).thetaV();
        dimensionedScalar mass("mass", dimMass,
            cloud_.constProps(moleculeI).mass());

        TvInitial.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "../meshCFD/Tv_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimTemperature, 0.0)
            )
        );

        nDInitial.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "../meshCFD/nD_" + molecules[moleculeI],
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
            new volScalarField
            (
                IOobject
                (
                    "../meshCFD/p_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                nDInitial[moleculeI] * k * TtInitial
            )
        );

        kappatrInitial.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "../meshCFD/kappatr_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        kappaveInitial.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "../meshCFD/kappave_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", kappatrInitial[moleculeI].dimensions(), 0.0)
            )
        );

        muInitial.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "../meshCFD/mu_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        diffFluxInitial.set
        (
            moleculeI,
            new volVectorField
            (
                IOobject
                (
                    "../meshCFD/J_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", dimMass / dimArea / dimTime, vector::zero)
            )
        );

        scalar rotDoF =
            cloud_.constProps()[moleculeI].rotationalDegreesOfFreedom();
        scalar vibDoF =
            cloud_.constProps()[moleculeI].vibrationalDegreesOfFreedom();
        scalar omega =
            cloud_.constProps()[moleculeI].omega();

        qt.set
        (
            moleculeI,
            new volVectorField
            (
                IOobject
                (
                    "qt_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", dimless, vector::zero)
            )
        );

        qr.set
        (
            moleculeI,
            new volVectorField
            (
                IOobject
                (
                    "qr_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", dimless, vector::zero)
            )
        );

        qv.set
        (
            moleculeI,
            new volVectorField
            (
                IOobject
                (
                    "qv_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", dimless, vector::zero)
            )
        );

        D.set
        (
            moleculeI,
            volVectorField
            (
                IOobject
                (
                    "D_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", dimless, vector::zero)
            )
        );

        tau.set
        (
            moleculeI,
            volTensorField
            (
                IOobject
                (
                    "tau_" + molecules[moleculeI],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedTensor("zero", dimless, tensor::zero)
            )
        );

        volVectorField dummyGradT
        (
            IOobject
            (
                "dummygradT" + molecules[moleculeI],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(TtInitial)
        );
        volVectorField dummyGradTv
        (
            IOobject
            (
                "dummygradTv" + molecules[moleculeI],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(TvInitial[moleculeI])
        );

        volTensorField dummyTau
        (
            IOobject
            (
                "dummyTau" + molecules[moleculeI],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            muInitial[moleculeI] * (fvc::grad(UInitial)
                + dev2(Foam::T(fvc::grad(UInitial))))
        );

        dummyGradT.write();
        dummyGradTv.write();
        dummyTau.write();

        forAll(mesh_.cells(), cellI)
        {
            if(nDInitial[moleculeI][cellI] > VSMALL)
            {
                scalar twoBeta = sqrt(2.0 * mass.value() / k.value() / TtInitial[cellI]);

                qt[moleculeI][cellI] = -(75.0 / (75.0 + 2.0 * (7.0 - 2.0
                    * omega) * rotDoF)) * (kappatrInitial[moleculeI][cellI]
//                    / pInitial[moleculeI][cellI])
                    / pTotInitial[cellI])
                    * dummyGradT[cellI] * twoBeta;

                if(rotDoF > SMALL)
                {
                    qr[moleculeI][cellI] = -(2.0 * (7.0 - 2.0 * omega) * rotDoF
                        / (75.0 + 2.0 * (7.0 - 2.0 * omega) * rotDoF))
                        * (kappatrInitial[moleculeI][cellI] * 2.0 / rotDoF
//                        / pInitial[moleculeI][cellI])
                        / pTotInitial[cellI])
                        * dummyGradT[cellI] * twoBeta;
                }

                if(vibDoF > SMALL)
                {
                    qv[moleculeI][cellI] = -(kappaveInitial[moleculeI][cellI]
//                        / pInitial[moleculeI][cellI])
                        / pTotInitial[cellI])
                        * dummyGradTv[cellI]
                        * twoBeta * (TtInitial[cellI]
                        * TvInitial[moleculeI][cellI] / pow(thetaV, 2))
                        * (exp(thetaV / TvInitial[moleculeI][cellI]) - 1.0)
                        * (exp(thetaV / TvInitial[moleculeI][cellI]) - 1.0)
                        / exp(thetaV / TvInitial[moleculeI][cellI]);
                }

                D[moleculeI][cellI] = diffFluxInitial[moleculeI][cellI]
                    * 0.5 * twoBeta / mass.value()
                    / nDInitial[moleculeI][cellI];

                tau[moleculeI][cellI] = dummyTau[cellI]
//                    / pInitial[moleculeI][cellI];
                    / pTotInitial[cellI];
            }
        }

        qt[moleculeI].write();
        qr[moleculeI].write();
        qv[moleculeI].write();
        D[moleculeI].write();
        tau[moleculeI].write();
    }


    forAll(molecules, moleculeI)
    {
        // Average number of particles per specie
        numberDensities[moleculeI] =
            nDInitial[moleculeI].weightedAverage(mesh_.V()).value();
        nDInitial[moleculeI] /= cloud_.nParticle();
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

                const dsmcParcel::constantProperties& cP
                    = cloud_.constProps(typeId);

                scalar numberDensity = nDInitial[i][cellI];
                scalar translationalTemperature = TtInitial[cellI];
                scalar rotationalTemperature = TtInitial[cellI];
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

                    vector U;

                    scalar ERot;
                    scalar EVib;

                    cloud_.generalisedChapmanEnskog
                    (
                        i,
                        translationalTemperature,
                        rotationalTemperature,
                        vibrationalTemperature,
                        cP.mass(),
                        D[i][cellI],
                        qt[i][cellI],
                        qr[i][cellI],
                        qv[i][cellI],
                        tau[i][cellI],
                        ERot,
                        EVib,
                        U
                    );

                    U += velocity;
                    
                    label newParcel = -1;
                    
                    label classification = 0;

                    label stuckToWall = 0;
                    
                    scalarField wallTemperature(4, 0.0);
                    
                    vectorField wallVectors(4, vector::zero);
                    
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
                        stuckToWall,
                        wallTemperature,
                        wallVectors,
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
