/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
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

Description

Measures overall temperature, including vibrational temperature, for a single species gas
or a gas mixture and writes the results to a volume scalar field that can be viewed in Paraview.

Translational, rotatational and vibrational temperature field will also be written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements class and are also written.

\*---------------------------------------------------------------------------*/

#include "dsmcVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVolFields, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVolFields, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void dsmcVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];

        if (isA<wallPolyPatch>(pPatch))
        {
            const vectorField& fC = pPatch.faceCentres();

            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei]
                    /mag(pPatch.faceAreas()[facei]);

                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei] - mesh_.points()[mesh_.faces()[pPatch.start() + facei][0]];
                t1_[patchi][facei] /= mag(t1_[patchi][facei]);

                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei] /= mag(t2_[patchi][facei]);
            }

            //Info << "n_[patchi]: " << n_[patchi] << endl;
            //Info << "t1_[patchi]: " << t1_[patchi] << endl;
            //Info << "t2_[patchi]: " << t2_[patchi] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVolFields::dsmcVolFields
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    sampleInterval_(1),
    sampleCounter_(0),
    mfpReferenceTemperature_(273.0),
    fieldName_(propsDict_.lookup("fieldName")),
    dsmcN_
    (
        IOobject
        (
            "dsmcN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    dsmcNMean_
    (
        IOobject
        (
            "dsmcNMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimPressure, 0.0)
    ),
    translationalT_
    (
        IOobject
        (
            "translationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    rotationalT_
    (
        IOobject
        (
            "rotationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    vibrationalT_
    (
        IOobject
        (
            "vibrationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    electronicT_
    (
        IOobject
        (
            "electronicT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    overallT_
    (
        IOobject
        (
            "overallT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    q_
    (
        IOobject
        (
            "surfaceHeatTransfer_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    tau_
    (
        IOobject
        (
            "surfaceShearStress_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimPressure, 0.0)
    ),
    meanFreePath_
    (
        IOobject
        (
            "variableHardSphereMeanFreePath_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    mfpCellRatio_
    (
        IOobject
        (
            "mfpCellRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    cellMfpRatio_
    (
        IOobject
        (
            "cellMfpRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    meanCollisionRate_
    (
        IOobject
        (
            "meanCollisionRate_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    meanCollisionTime_
    (
        IOobject
        (
            "meanCollisionTime_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, 1, 0, 0), 0.0)
    ),
    meanCollisionTimeTimeStepRatio_
    (
        IOobject
        (
            "mctTimeStepRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimless, 0.0)
    ),
    cr_
    (
        IOobject
        (
            "measuredCollisionRate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    meanCollisionSeparation_
    (
        IOobject
        (
            "meanCollisionSeparation_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    SOF_
    (
        IOobject
        (
            "SOF_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    Ma_
    (
        IOobject
        (
            "Ma_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIDistribution_
    (
        IOobject
        (
            "classIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIIDistribution_
    (
        IOobject
        (
            "classIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIIIDistribution_
    (
        IOobject
        (
            "classIIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    densityError_
    (
        IOobject
        (
            "densityError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    velocityError_
    (
        IOobject
        (
            "velocityError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    temperatureError_
    (
        IOobject
        (
            "temperatureError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureError_
    (
        IOobject
        (
            "pressureError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    UMean_
    (
        IOobject
        (
            "UMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    fD_
    (
        IOobject
        (
            "fD_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -3, 0, 0),
            vector::zero
        )
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    nTimeSteps_(0.0),
    typeIds_(),
    dsmcNCum_(mesh_.nCells(), 0.0),
    dsmcNInstantaneous_(mesh_.nCells(), 0.0),
    nCum_(mesh_.nCells(), 0.0),
    dsmcNWithRotDofCum_(mesh_.nCells(), 0.0),
    dsmcNEleCum_(mesh_.nCells(), 0.0),
    dsmcMCum_(mesh_.nCells(), 0.0),
    mCum_(mesh_.nCells(), 0.0),
    dsmcLinKinEnCum_(mesh_.nCells(), 0.0),
    linKinEnCum_(mesh_.nCells(), 0.0),
    dsmcRotEnCum_(mesh_.nCells(), 0.0),
    dsmcRotDofCum_(mesh_.nCells(), 0.0),
    dsmcMuuCum_(mesh_.nCells(), 0.0),
    dsmcMuvCum_(mesh_.nCells(), 0.0),
    dsmcMuwCum_(mesh_.nCells(), 0.0),
    dsmcMvvCum_(mesh_.nCells(), 0.0),
    dsmcMvwCum_(mesh_.nCells(), 0.0),
    dsmcMwwCum_(mesh_.nCells(), 0.0),
    dsmcMccCum_(mesh_.nCells(), 0.0),
    dsmcMccuCum_(mesh_.nCells(), 0.0),
    dsmcMccvCum_(mesh_.nCells(), 0.0),
    dsmcMccwCum_(mesh_.nCells(), 0.0),
    dsmcEuCum_(mesh_.nCells(), 0.0),
    dsmcEvCum_(mesh_.nCells(), 0.0),
    dsmcEwCum_(mesh_.nCells(), 0.0),
    dsmcECum_(mesh_.nCells(), 0.0),
    totalvDof_(mesh_.nCells(), 0.0),
    dsmcNClassICum_(mesh_.nCells(), 0.0),
    dsmcNClassIICum_(mesh_.nCells(), 0.0),
    dsmcNClassIIICum_(mesh_.nCells(), 0.0),
    collisionSeparation_(mesh_.nCells(), 0.0),
    dsmcNCollsCum_(mesh_.nCells(), 0.0),
    dsmcMomentumCum_(mesh.nCells(), vector::zero),
    momentumCum_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    dsmcVibEnSpeciesModeCum_(),
    dsmcEleEnSpeciesCum_(),
    dsmcNSpeciesCum_(),
    nSpeciesCum_(),
    dsmcMccSpeciesCum_(),
    vibT_(),
    dsmcNGroundEleLvlSpeciesCum_(),
    dsmcNFirstEleLvlSpeciesCum_(),
    vDof_(),
    mfp_(),
    mcr_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    qBF_(),
    vibTxvDofBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    speciesRhoNElecBF_(),
    momentumBF_(),
    fDBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    speciesRhoNBF_(),
    mccSpeciesBF_(),
    vibTBF_(),
    vDofBF_(),
    evmsBF_(),
    n_(),
    t1_(),
    t2_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVolFields::~dsmcVolFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVolFields::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "volFieldsMethod_"+fieldName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("nTimeSteps", nTimeSteps_);

    // DSMC parcel related cumulative values
    dict.readIfPresent("dsmcNCum", dsmcNCum_);
    dict.readIfPresent("dsmcMCum", dsmcMCum_);

    dict.readIfPresent("dsmcLinKinEnCum", dsmcLinKinEnCum_);
    dict.readIfPresent("dsmcMomentumCum", dsmcMomentumCum_);
    dict.readIfPresent("dsmcRotEnCum", dsmcRotEnCum_);
    dict.readIfPresent("dsmcRotDofCum", dsmcRotDofCum_);
    dict.readIfPresent("dsmcEleEnSpeciesCum", dsmcEleEnSpeciesCum_);
    dict.readIfPresent("dsmcNSpeciesCum", dsmcNSpeciesCum_);
    dict.readIfPresent("dsmcMccSpeciesCum", dsmcMccSpeciesCum_);
    dict.readIfPresent("dsmcMuuCum", dsmcMuuCum_);
    dict.readIfPresent("dsmcMuvCum", dsmcMuvCum_);
    dict.readIfPresent("dsmcMuwCum", dsmcMuwCum_);
    dict.readIfPresent("dsmcMvvCum", dsmcMvvCum_);
    dict.readIfPresent("dsmcMvwCum", dsmcMvwCum_);
    dict.readIfPresent("dsmcMwwCum", dsmcMwwCum_);
    dict.readIfPresent("dsmcMccCum", dsmcMccCum_);
    dict.readIfPresent("dsmcMccuCum", dsmcMccuCum_);
    dict.readIfPresent("dsmcMccvCum", dsmcMccvCum_);
    dict.readIfPresent("dsmcMccwCum", dsmcMccwCum_);
    dict.readIfPresent("dsmcEuCum", dsmcEuCum_);
    dict.readIfPresent("dsmcEvCum", dsmcEvCum_);
    dict.readIfPresent("dsmcEwCum", dsmcEwCum_);
    dict.readIfPresent("dsmcECum", dsmcECum_);
    dict.readIfPresent("dsmcNWithRotDofCum", dsmcNWithRotDofCum_);
    dict.readIfPresent("dsmcNEleCum", dsmcNEleCum_);
    dict.readIfPresent("dsmcNGroundEleLvlSpeciesCum", dsmcNGroundEleLvlSpeciesCum_);
    dict.readIfPresent("dsmcNFirstEleLvlSpeciesCum", dsmcNFirstEleLvlSpeciesCum_);
    dict.readIfPresent("dsmcNClassICum", dsmcNClassICum_);
    dict.readIfPresent("dsmcNClassIICum", dsmcNClassIICum_);
    dict.readIfPresent("dsmcNClassIIICum", dsmcNClassIIICum_);
    dict.readIfPresent("dsmcVibEnSpeciesModeCum", dsmcVibEnSpeciesModeCum_);
    dict.readIfPresent("dsmcNCollsCum", dsmcNCollsCum_);

    // cumulative values
    dict.readIfPresent("nCum", nCum_);
    dict.readIfPresent("mCum", mCum_);
    dict.readIfPresent("nSpeciesCum", nSpeciesCum_);
    dict.readIfPresent("momentumCum", momentumCum_);
    dict.readIfPresent("linKinEnCum", linKinEnCum_);
    dict.readIfPresent("totalvDof", totalvDof_);
}


void dsmcVolFields::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFieldsMethod_"+fieldName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("nTimeSteps", nTimeSteps_);

        // DSMC parcel related cumulative values
        dict.add("dsmcNCum", dsmcNCum_);
        dict.add("dsmcMCum", dsmcMCum_);

        dict.add("dsmcLinKinEnCum", dsmcLinKinEnCum_);
        dict.add("dsmcMomentumCum", dsmcMomentumCum_);
        dict.add("dsmcRotEnCum", dsmcRotEnCum_);
        dict.add("dsmcRotDofCum", dsmcRotDofCum_);
        dict.add("dsmcEleEnSpeciesCum", dsmcEleEnSpeciesCum_);
        dict.add("dsmcNSpeciesCum", dsmcNSpeciesCum_);
        dict.add("dsmcMccSpeciesCum", dsmcMccSpeciesCum_);
        dict.add("dsmcMuuCum", dsmcMuuCum_);
        dict.add("dsmcMuvCum", dsmcMuvCum_);
        dict.add("dsmcMuwCum", dsmcMuwCum_);
        dict.add("dsmcMvvCum", dsmcMvvCum_);
        dict.add("dsmcMvwCum", dsmcMvwCum_);
        dict.add("dsmcMwwCum", dsmcMwwCum_);
        dict.add("dsmcMccCum", dsmcMccCum_);
        dict.add("dsmcMccuCum", dsmcMccuCum_);
        dict.add("dsmcMccvCum", dsmcMccvCum_);
        dict.add("dsmcMccwCum", dsmcMccwCum_);
        dict.add("dsmcEuCum", dsmcEuCum_);
        dict.add("dsmcEvCum", dsmcEvCum_);
        dict.add("dsmcEwCum", dsmcEwCum_);
        dict.add("dsmcECum", dsmcECum_);
        dict.add("dsmcNWithRotDofCum", dsmcNWithRotDofCum_);
        dict.add("dsmcNEleCum", dsmcNEleCum_);
        dict.add("dsmcNGroundEleLvlSpeciesCum", dsmcNGroundEleLvlSpeciesCum_);
        dict.add("dsmcNFirstEleLvlSpeciesCum", dsmcNFirstEleLvlSpeciesCum_);
        dict.add("dsmcNClassICum", dsmcNClassICum_);
        dict.add("dsmcNClassIICum", dsmcNClassIICum_);
        dict.add("dsmcNClassIIICum", dsmcNClassIIICum_);
        dict.add("dsmcVibEnSpeciesModeCum", dsmcVibEnSpeciesModeCum_);
        dict.add("dsmcNCollsCum", dsmcNCollsCum_);

        // cumulative values
        dict.add("nCum", nCum_);
        dict.add("mCum", mCum_);
        dict.add("nSpeciesCum", nSpeciesCum_);
        dict.add("momentumCum", momentumCum_);
        dict.add("linKinEnCum", linKinEnCum_);
        dict.add("totalvDof", totalvDof_);

        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();

        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}


//- Initial configuration
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;

    const List<word>& molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if (findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName = moleculesReduced[i];

        const label typeId = findIndex(cloud_.typeIdList(), moleculeName);

        if (typeId == -1)
        {
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    vibT_.setSize(typeIds_.size());
    dsmcNGroundEleLvlSpeciesCum_.setSize(typeIds_.size());
    dsmcNFirstEleLvlSpeciesCum_.setSize(typeIds_.size());
    vDof_.setSize(typeIds_.size());
    dsmcVibEnSpeciesModeCum_.setSize(typeIds_.size());
    dsmcEleEnSpeciesCum_.setSize(typeIds_.size());
    dsmcNSpeciesCum_.setSize(typeIds_.size());
    nSpeciesCum_.setSize(typeIds_.size());
    dsmcMccSpeciesCum_.setSize(typeIds_.size());
    mfp_.setSize(typeIds_.size());
    mcr_.setSize(typeIds_.size());

    boundaryCells_.setSize(mesh_.boundaryMesh().size());

    forAll(typeIds_, i)
    {
        vibT_[i].setSize(mesh_.nCells());
        dsmcNGroundEleLvlSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcNFirstEleLvlSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        vDof_[i].setSize(mesh_.nCells(), 0.0);
        dsmcEleEnSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcNSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        nSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcMccSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        mfp_[i].setSize(mesh_.nCells(), 0.0);
        mcr_[i].setSize(mesh_.nCells(), 0.0);

        dsmcVibEnSpeciesModeCum_[i].setSize
        (
            cloud_.constProps(typeIds_[i]).nVibrationalModes()
        );

        forAll(dsmcVibEnSpeciesModeCum_[i], j)
        {
            dsmcVibEnSpeciesModeCum_[i][j].setSize(mesh_.nCells(), 0.0);
        }
    }

    forAll(boundaryCells_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        boundaryCells_[j].setSize(patch.size());

        forAll(boundaryCells_[j], k)
        {
            boundaryCells_[j][k] = patch.faceCells()[k];
        }
    }

    //- boundary fields initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());
    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());

    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        fDBF_[j].setSize(patch.size(), vector::zero);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);

        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }

    calculateWallUnitVectors();

    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());
    evmsBF_.setSize(typeIds_.size());

    forAll(typeIds_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());
        evmsBF_[i].setSize(cloud_.constProps(typeIds_[i]).thetaV().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }

        forAll(evmsBF_[i], mode)
        {
            evmsBF_[i][mode].setSize(mesh_.boundaryMesh().size());
            forAll(evmsBF_[i][mode], j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                evmsBF_[i][mode][j].setSize(patch.size(), 0.0);
            }
        }
    }

    sampleInterval_ = propsDict_.lookupOrDefault("sampleInterval", 1);

    measureClassifications_ =
        propsDict_.lookupOrDefault<bool>("measureClassifications", false);

    measureErrors_ = propsDict_.lookupOrDefault<bool>("measureErrors", false);

    densityOnly_ = propsDict_.lookupOrDefault<bool>("densityOnly", false);

    measureHeatFluxShearStress_ =
        propsDict_.lookupOrDefault<bool>("measureHeatFluxShearStress", false);

    measureMeanFreePath_ =
        propsDict_.lookupOrDefault<bool>("measureMeanFreePath", false);

    mfpReferenceTemperature_ =
        propsDict_.lookupOrDefault<scalar>("mfpReferenceTemperature", 273.0);

    averagingAcrossManyRuns_ =
        propsDict_.lookupOrDefault<bool>("averagingAcrossManyRuns", false);

    //- read in stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        Info<< "Averaging across many runs enabled." << nl
            << endl;

        readIn();
    }
}


void dsmcVolFields::calculateField()
{
    sampleCounter_++;

    dsmcNInstantaneous_ = 0.0;

    const scalar kB = physicoChemical::k.value();
    const scalar NAvo = physicoChemical::NA.value();

    if (sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;
        const scalar nAvTimeSteps = nTimeSteps_;

        if (densityOnly_)
        {
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label iD = findIndex(typeIds_, p.typeId());

                if (iD != -1 && p.isFree())
                {
                    const label cell = p.cell();
                    const scalar nParticles = cloud_.nParticles(cell);
                    const scalar mass = cloud_.constProps(p.typeId()).mass();

                    // cumulative number of DSMC parcels
                    dsmcNCum_[cell] += 1.0;
                    // instantaneous number of DSMC parcels in this time step
                    dsmcNInstantaneous_[cell] += 1.0;

                    // cumulative number of simulated real particles
                    nCum_[cell] += nParticles;
                    // cumulative mass
                    mCum_[cell] += mass*nParticles;
                }
            }
        }
        else
        {
            //scalar timer = mesh_.time().elapsedCpuTime(); // TODO VINCENT

            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label iD = findIndex(typeIds_, p.typeId());

                if (iD != -1 && p.isFree())
                {
                    const dsmcParcel::constantProperties& cP = cloud_.constProps(p.typeId());

                    const label cell = p.cell();
                    const scalar nParticles = cloud_.nParticles(cell);
                    const scalar mass = cP.mass();
                    const scalar massBySqMagU = mass*(p.U() & p.U());
                    const scalar rotationalDof = cP.rotationalDegreesOfFreedom();
                    const scalarList& electronicEnergies = cP.electronicEnergyList();

                    scalar vibEn = 0.0;
                    forAll(cP.thetaV(), mode)
                    {
                        const scalar eVib_m = cP.eVib_m(mode, p.vibLevel()[mode]);
                        dsmcVibEnSpeciesModeCum_[iD][mode][cell] += eVib_m;
                        vibEn += eVib_m;
                    }

                    // cumulative number of DSMC parcels
                    dsmcNCum_[cell] += 1.0;
                    // instantaneous number of DSMC parcels in this time step
                    dsmcNInstantaneous_[cell] += 1.0;
                    // cumulative mass of the DSMC parcels
                    dsmcMCum_[cell] += mass;
                    // cumulative linear kinetic energy of the DSMC parcels
                    dsmcLinKinEnCum_[cell] += massBySqMagU;
                    // cumulative momentum of the DSMC parcels
                    dsmcMomentumCum_[cell] += mass*p.U();
                    // cumulative rotational energy of the DSMC parcels
                    dsmcRotEnCum_[cell] += p.ERot();
                    // cumulative rotational degrees of freedom of the DSMC
                    // parcels
                    dsmcRotDofCum_[cell] += rotationalDof;
                    // cumulative electronic energy of the DSMC parcels of
                    // each species
                    dsmcEleEnSpeciesCum_[iD][cell] +=
                        electronicEnergies[p.ELevel()];
                    // cumulative number of DSMC parcels of each species
                    dsmcNSpeciesCum_[iD][cell] += 1.0;
                    dsmcMccSpeciesCum_[iD][cell] += massBySqMagU;
                    // cumulative number of simulated real particles of each
                    // species
                    nSpeciesCum_[iD][cell] += nParticles;
                    // cumulative number of simulated real particles
                    nCum_[cell] += nParticles;
                    // cumulative mass
                    mCum_[cell] += mass*nParticles;
                    // cumulative momentum of the simulated real particles
                    momentumCum_[cell] += mass*(p.U())*nParticles;
                    // cumulative linear kinetic energy of the simulated
                    // particles
                    linKinEnCum_[cell] += massBySqMagU*nParticles;

                    // cumulative counters for mass times squared velocity
                    // components of the DSMC parcels
                    // notation: p.U = c = (u v w)
                    dsmcMuuCum_[cell] += mass*sqr(p.U().x());
                    dsmcMuvCum_[cell] += mass*p.U().x()*p.U().y();
                    dsmcMuwCum_[cell] += mass*p.U().x()*p.U().z();
                    dsmcMvvCum_[cell] += mass*sqr(p.U().y());
                    dsmcMvwCum_[cell] += mass*p.U().y()*p.U().z();
                    dsmcMwwCum_[cell] += mass*sqr(p.U().z());

                    dsmcMccCum_[cell] += massBySqMagU;
                    dsmcMccuCum_[cell] += massBySqMagU*(p.U().x());
                    dsmcMccvCum_[cell] += massBySqMagU*(p.U().y());
                    dsmcMccwCum_[cell] += massBySqMagU*(p.U().z());

                    // cumulative energy times velocity components of the DSMC
                    // parcels
                    dsmcEuCum_[cell] += (p.ERot() + vibEn)*p.U().x(); // TODO missing Eel
                    dsmcEvCum_[cell] += (p.ERot() + vibEn)*p.U().y(); // TODO missing Eel
                    dsmcEwCum_[cell] += (p.ERot() + vibEn)*p.U().z(); // TODO missing Eel
                    dsmcECum_[cell] += (p.ERot() + vibEn); // TODO missing Eel

                    if (rotationalDof > 0)
                    {
                        dsmcNWithRotDofCum_[cell] += 1.0;
                    }

                    const label& nElecLevels = cP.nElectronicLevels();

                    if (nElecLevels > 1)
                    {
                        dsmcNEleCum_[cell] += 1.0;

                        if (p.ELevel() == 0)
                        {
                            dsmcNGroundEleLvlSpeciesCum_[iD][cell]++;
                        }
                        if (p.ELevel() == 1)
                        {
                            dsmcNFirstEleLvlSpeciesCum_[iD][cell]++;
                        }
                    }

                    if (measureClassifications_)
                    {
                        const label& classification = p.classification();

                        if (classification == 0)
                        {
                            dsmcNClassICum_[cell] += 1.0;
                        }

                        if (classification == 1)
                        {
                            dsmcNClassIICum_[cell] += 1.0;
                        }

                        if (classification == 2)
                        {
                            dsmcNClassIIICum_[cell] += 1.0;
                        }
                    }
                }
            }

            //Info<< "fields myCal" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl; // TODO VINCENT

            //- Obtain collision quality measurements and mixture translational
            //  temperature
            forAll(cloud_.cellPropMeasurements().collisionSeparation(), cell)
            {
                collisionSeparation_[cell] +=
                    cloud_.cellPropMeasurements().collisionSeparation()[cell];
                dsmcNCollsCum_[cell] += cloud_.cellPropMeasurements().nColls()[cell];

                if (dsmcNCum_[cell] > 1e-3)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];

                    dsmcNMean_[cell] = dsmcNCum_[cell]/nAvTimeSteps;

                    const scalar rhoNMean = nCum_[cell]
                        /(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean = mCum_[cell]
                        /(nAvTimeSteps*cellVolume);

                    rhoN_[cell] = rhoNMean;
                    rhoM_[cell] = rhoMMean;

                    UMean_[cell] = momentumCum_[cell]
                        /mCum_[cell];

                    const scalar linearKEMean = 0.5
                        *linKinEnCum_[cell]
                        /(cellVolume*nAvTimeSteps);

                    //- Translational temperature
                    translationalT_[cell] =
                        2.0/(3.0*kB*rhoNMean)
                       *(
                            linearKEMean - 0.5*rhoMMean
                           *(
                                UMean_[cell] & UMean_[cell]
                            )
                        );

                    p_[cell] = rhoNMean*kB*translationalT_[cell];
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_[cell] = 0.001;
                    //dsmcN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;
                    p_[cell] = 0.0;
                }
            }

            //- Obtain boundary measurements
            forAll(typeIds_, i)
            {
                const label typeId = typeIds_[i];

                forAll(mesh_.boundaryMesh(), j)
                {
                    forAll(mesh_.boundaryMesh()[j], k)
                    {
                        rhoNBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[typeId][j][k];
                        rhoMBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoMBF()[typeId][j][k];
                        linearKEBF_[j][k] += cloud_.boundaryFluxMeasurements().linearKEBF()[typeId][j][k];
                        momentumBF_[j][k] += cloud_.boundaryFluxMeasurements().momentumBF()[typeId][j][k];
                        rotationalEBF_[j][k] += cloud_.boundaryFluxMeasurements().rotationalEBF()[typeId][j][k];
                        rotationalDofBF_[j][k] += cloud_.boundaryFluxMeasurements().rotationalDofBF()[typeId][j][k];
                        qBF_[j][k] += cloud_.boundaryFluxMeasurements().qBF()[typeId][j][k];
                        fDBF_[j][k] += cloud_.boundaryFluxMeasurements().fDBF()[typeId][j][k];
                        speciesRhoNBF_[i][j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[typeId][j][k];
                        vibrationalEBF_[i][j][k] += cloud_.boundaryFluxMeasurements().vibrationalEBF()[typeId][j][k];
                        electronicEBF_[i][j][k] += cloud_.boundaryFluxMeasurements().electronicEBF()[typeId][j][k];
                        mccSpeciesBF_[i][j][k] += cloud_.boundaryFluxMeasurements().mccSpeciesBF()[typeId][j][k];
                        speciesRhoNIntBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNIntBF()[typeId][j][k];
                        speciesRhoNElecBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNElecBF()[typeId][j][k];
                    }
                }

                forAll(evmsBF_[i], mode)
                {
                    forAll(mesh_.boundaryMesh(), j)
                    {
                        forAll(mesh_.boundaryMesh()[j], k)
                        {
                            evmsBF_[i][mode][j][k] +=
                                cloud_.boundaryFluxMeasurements()
                                    .evmsBF()[typeId][mode][j][k];
                        }
                    }
                }
            }
        }

        sampleCounter_ = 0;
    }

    if (time_.time().outputTime())
    {
        const scalar nAvTimeSteps = nTimeSteps_;

        if (densityOnly_)
        {
            forAll(dsmcNCum_, cell)
            {
                if (dsmcNCum_[cell] > SMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];

                    dsmcNMean_[cell] = dsmcNCum_[cell]/nAvTimeSteps;

                    rhoN_[cell] = nCum_[cell]
                        /(nAvTimeSteps*cellVolume);

                    rhoM_[cell] = mCum_[cell]
                        /(nAvTimeSteps*cellVolume);
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_[cell] = 0.001;
                    //dsmcN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                }

                if (dsmcNInstantaneous_[cell] > SMALL)
                {
                    dsmcN_[cell] = dsmcNInstantaneous_[cell];
                }
                else
                {
                    dsmcN_[cell] = 0.001;
                }
            }
        }
        else
        {
            const label nCells = mesh_.nCells();

            vibrationalT_.primitiveFieldRef() = 0.0;
            scalarField vibTForOverallT(nCells, 0.0);

            //- Heat capacity at constant volume/(0.5*kB), trans-rotational
            scalarField molarCv_transrot(nCells, 0.0);
            //- Heat capacity at constant pressure/(0.5*kB), trans-rotational
            scalarField molarCp_transrot(nCells, 0.0);
            scalarField molecularMass(nCells, 0.0);
            scalarField particleCv(nCells, 0.0);
            scalarField totalvDofOverall(nCells, 0.0);

            forAll(dsmcNCum_, cell)
            {
                /*if (dsmcNCum_[cell] > 1e-3) // TODO moved up already
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];

                    dsmcNMean_[cell] = dsmcNCum_[cell]/nAvTimeSteps;

                    //dsmcN_[cell] = dsmcNInstantaneous_[cell];

                    const scalar rhoNMean = nCum_[cell]
                        /(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean = mCum_[cell]
                        /(nAvTimeSteps*cellVolume);

                    rhoN_[cell] = rhoNMean;
                    rhoM_[cell] = rhoMMean;

                    UMean_[cell] = momentumCum_[cell]
                        /mCum_[cell];

                    const scalar linearKEMean = 0.5
                        *linKinEnCum_[cell]
                        /(cellVolume*nAvTimeSteps);

                    translationalT_[cell] =
                        2.0/(3.0*kB*rhoNMean)
                       *(
                            linearKEMean - 0.5*rhoMMean
                           *(
                                UMean_[cell] & UMean_[cell]
                            )
                        );

                    p_[cell] = rhoNMean*kB*translationalT_[cell];
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_[cell] = 0.001;
                    //dsmcN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;
                    p_[cell] = 0.0;
                }*/

                if (dsmcNInstantaneous_[cell] > SMALL)
                {
                    dsmcN_[cell] = dsmcNInstantaneous_[cell];
                }
                else
                {
                    dsmcN_[cell] = 0.001;
                }

                //- Rotational energy mode
                const scalar totalrDof
                (
                    dsmcNCum_[cell] > SMALL
                  ? dsmcRotDofCum_[cell]/dsmcNCum_[cell]
                  : 0.0
                );

                rotationalT_[cell] =
                (
                    dsmcRotDofCum_[cell] > SMALL
                  ? 2.0*dsmcRotEnCum_[cell]/(kB*dsmcRotDofCum_[cell])
                  : 0.0
                );

                //- Vibrational energy mode
                scalarList degreesOfFreedomSpecies(typeIds_.size(), 0.0);
                scalarList vibTID(typeIds_.size(), 0.0);
                List<scalarList> degreesOfFreedomMode(typeIds_.size());
                List<scalarList> vibTMode(typeIds_.size());

                forAll(degreesOfFreedomMode, i)
                {
                    degreesOfFreedomMode[i].setSize
                    (
                        cloud_.constProps(typeIds_[i])
                            .nVibrationalModes(), 0.0
                    );
                    vibTMode[i].setSize
                    (
                        cloud_.constProps(typeIds_[i])
                            .nVibrationalModes(), 0.0
                    );
                }

                forAll(typeIds_, i)
                {
                    forAll(dsmcVibEnSpeciesModeCum_[i], mode)
                    {
                        if (dsmcVibEnSpeciesModeCum_[i][mode][cell] > VSMALL
                            && dsmcNSpeciesCum_[i][cell] > SMALL
                            && degreesOfFreedomMode.size() > SMALL)
                        {
                            const scalar thetaV =
                                cloud_.constProps(typeIds_[i]).thetaV_m(mode);

                            const scalar vibrationalEMean =
                                dsmcVibEnSpeciesModeCum_[i][mode][cell]
                               /dsmcNSpeciesCum_[i][cell];

                            const scalar iMean = vibrationalEMean/(kB*thetaV);

                            vibTMode[i][mode] = thetaV/log(1.0 + 1.0/iMean);

                            degreesOfFreedomMode[i][mode] =
                                2.0*thetaV/vibTMode[i][mode]
                              /(exp(thetaV/vibTMode[i][mode]) - 1.0);

                            degreesOfFreedomSpecies[i] += degreesOfFreedomMode[i][mode];
                        }
                    }

                    forAll(degreesOfFreedomMode[i], mode)
                    {
                        if (degreesOfFreedomSpecies[i] > SMALL)
                        {
                            vibTID[i] += vibTMode[i][mode]
                                *degreesOfFreedomMode[i][mode]
                                /degreesOfFreedomSpecies[i];
                        }
                    }

                    totalvDof_[cell] += degreesOfFreedomSpecies[i];

                    if
                    (
                         dsmcNWithRotDofCum_[cell] > VSMALL
                      && dsmcNCum_[cell] > VSMALL
                      && dsmcNSpeciesCum_[i][cell] > SMALL
                    )
                    {
                        const scalar fraction = dsmcNSpeciesCum_[i][cell]
                            /dsmcNWithRotDofCum_[cell];

                        const scalar fractionOverall = dsmcNSpeciesCum_[i][cell]
                            /dsmcNCum_[cell];

                        totalvDofOverall[cell] +=
                            totalvDof_[cell]*fractionOverall/fraction;

                        //- TODO
                        vibrationalT_[cell] += vibTID[i]*fraction;
                    }
                }

                //- Electronic energy mode
                scalar totalEDof = 0.0;
                scalar elecT = 0.0;

                forAll(dsmcNSpeciesCum_, i) // TODO
                {
                    /*const scalarList& electronicEnergies =
                        cloud_.constProps(typeIds_[i]).electronicEnergyList();
                    const labelList& degeneracies =
                        cloud_.constProps(typeIds_[i]).electronicDegeneracyList();

                    if
                    (    dsmcNGroundEleLvlSpeciesCum_[i][cell] > SMALL
                      && dsmcNFirstEleLvlSpeciesCum_[i][cell] > SMALL
                      && dsmcNFirstEleLvlSpeciesCum_[i][cell]*degeneracies[0]
                            != dsmcNGroundEleLvlSpeciesCum_[i][cell]*degeneracies[1]
                    )
                    {
                        const scalar elecTID =
                            (electronicEnergies[1] - electronicEnergies[0])
                           /(
                                kB*log
                                (
                                    dsmcNGroundEleLvlSpeciesCum_[i][cell]*degeneracies[1]
                                   /(dsmcNFirstEleLvlSpeciesCum_[i][cell]*degeneracies[0])
                                )
                            );

                        const scalar fraction = dsmcNSpeciesCum_[i][cell]
                            /dsmcNEleCum_[cell];

                        if (elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }

                        const scalar eDof = 2.0*dsmcEleEnSpeciesCum_[i][cell]
                            /dsmcNSpeciesCum_[i][cell]/(kB*elecTID);

                        totalEDof += fraction*eDof;
                    }*/

                     label nElectronicLevels = cloud_.constProps(typeIds_[i]).nElectronicLevels();

                     if (nElectronicLevels > 1 && dsmcNSpeciesCum_[i][cell] > SMALL && dsmcNEleCum_[cell] > SMALL)
                     {
                         const scalarList& electronicEnergies = cloud_.constProps(typeIds_[i]).electronicEnergyList();
                         const labelList& degeneracies = cloud_.constProps(typeIds_[i]).electronicDegeneracyList();

                         const scalar translationalTSpecies =
                            1.0/(3.0*kB)
                           *(
                               dsmcMccSpeciesCum_[i][cell]/dsmcNSpeciesCum_[i][cell]
                             - (
                                   cloud_.constProps(typeIds_[i]).mass()
                                  *mag(UMean_[cell])*mag(UMean_[cell])
                               )
                           );

                         const scalar fraction = dsmcNSpeciesCum_[i][cell]/dsmcNEleCum_[cell];

                         //Info << "translationalTSpecies" << tab << translationalTSpecies << endl;
                         //Info << "fraction" << tab << fraction << endl;

                         if (translationalTSpecies > SMALL && dsmcEleEnSpeciesCum_[i][cell] > VSMALL)
                         {
                             scalar sum1 = 0.0;
                             scalar sum2 = 0.0;

                             forAll(electronicEnergies, ii)
                             {
                                 sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(kB*translationalTSpecies));
                                 sum2 += degeneracies[ii]*electronicEnergies[ii]/(kB*translationalTSpecies)
                                             *exp(-electronicEnergies[ii]/(kB*translationalTSpecies));
                             }

                             if (sum1 > VSMALL && sum2 > VSMALL)
                             {
                                 const scalar electronicTSpecies =
                                    dsmcEleEnSpeciesCum_[i][cell]
                                   /(kB*dsmcNSpeciesCum_[i][cell])*sum1/sum2;

                                 if (electronicTSpecies > SMALL && electronicTSpecies < GREAT)
                                 {
                                     elecT += fraction*electronicTSpecies;

                                     const scalar eDof =
                                        2.0*dsmcEleEnSpeciesCum_[i][cell]/dsmcNSpeciesCum_[i][cell]
                                       /(kB*translationalTSpecies);

                                     totalEDof += fraction*eDof;
                                 }
                             }
                         }
                     }
                }

                electronicT_[cell] = elecT;
//Info << "totalEDof" << tab << totalEDof << tab << "electronicT_[cell]" << tab << electronicT_[cell] << endl;


                /*scalar totalEDof = 0.0;
                scalar elecT = 0.0;

                forAll(dsmcNSpeciesCum_, iD)
                {
                    const scalarList& electronicEnergies =
                        cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const labelList& degeneracies =
                        cloud_.constProps(typeIds_[iD]).electronicDegeneracyList();

                    if
                    (
                        dsmcNGroundEleLvlSpeciesCum_[iD][cell] > VSMALL
                     && dsmcNFirstEleLvlSpeciesCum_[iD][cell] > VSMALL
                     && dsmcNFirstEleLvlSpeciesCum_[iD][cell]*degeneracies[0]
                            != dsmcNGroundEleLvlSpeciesCum_[iD][cell]*degeneracies[1]
                    )
                    {
                        const scalar elecTID =
                            (electronicEnergies[1] - electronicEnergies[0])/
                            (
                                kB*
                                log
                                (
                                    (dsmcNGroundEleLvlSpeciesCum_[iD][cell]*degeneracies[1])/
                                    (dsmcNFirstEleLvlSpeciesCum_[iD][cell]*degeneracies[0])
                                )
                            );


                        const scalar fraction = dsmcNSpeciesCum_[iD][cell]/dsmcNEleCum_[cell];

                        if(elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }

                        const scalar eDof =
                            (
                                2.0*(dsmcEleEnSpeciesCum_[iD][cell]/
                             dsmcNSpeciesCum_[iD][cell])
                            )
                            /
                            (
                                kB*elecTID
                            );

                        totalEDof += fraction*eDof;
                    }
                }

                electronicT_[cell] = elecT;*/

                //- Overall temperature
                overallT_[cell] =
                    (
                        3.0*translationalT_[cell]
                        + totalrDof*rotationalT_[cell]
                        + totalvDof_[cell]*vibrationalT_[cell]
                        + totalEDof*electronicT_[cell]
                    ) /
                    (3.0 + totalrDof + totalvDof_[cell] + totalEDof);


                if (measureHeatFluxShearStress_)
                {
                    if (dsmcNCum_[cell] > SMALL)
                    {
                        pressureTensor_[cell].xx() =
                            rhoN_[cell]*
                            (
                                dsmcMuuCum_[cell]/dsmcNCum_[cell]
                              - dsmcMCum_[cell]/dsmcNCum_[cell]
                                  *sqr(UMean_[cell].x())
                            );
                        pressureTensor_[cell].xy() =
                            rhoN_[cell]*
                            (
                                dsmcMuvCum_[cell]/dsmcNCum_[cell]
                              - dsmcMCum_[cell]/dsmcNCum_[cell]
                                    *UMean_[cell].x()*UMean_[cell].y()
                            );
                        pressureTensor_[cell].xz() =
                            rhoN_[cell]*
                            (
                                dsmcMuwCum_[cell]/dsmcNCum_[cell]
                              - dsmcMCum_[cell]/dsmcNCum_[cell]
                                  *UMean_[cell].x()*UMean_[cell].z()
                            );

                        pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                        pressureTensor_[cell].yy() =
                            rhoN_[cell]*
                            (
                                dsmcMvvCum_[cell]/(dsmcNCum_[cell])
                              - dsmcMCum_[cell]/dsmcNCum_[cell]
                                    *sqr(UMean_[cell].y())
                            );
                        pressureTensor_[cell].yz() =
                            rhoN_[cell]*
                            (
                                dsmcMvwCum_[cell]/dsmcNCum_[cell]
                              - dsmcMCum_[cell]/dsmcNCum_[cell]
                                    *UMean_[cell].y()*UMean_[cell].z()
                            );

                        pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                        pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                        pressureTensor_[cell].zz() =
                            rhoN_[cell]*
                            (
                                dsmcMwwCum_[cell]/dsmcNCum_[cell]
                              - dsmcMCum_[cell]/dsmcNCum_[cell]
                                    *sqr(UMean_[cell].z())
                            );

                        const scalar scalarPressure =
                            1.0/3.0*
                            (
                                pressureTensor_[cell].xx()
                              + pressureTensor_[cell].yy()
                              + pressureTensor_[cell].zz()
                            );

                        shearStressTensor_[cell] = -pressureTensor_[cell];
                        shearStressTensor_[cell].xx() += scalarPressure;
                        shearStressTensor_[cell].yy() += scalarPressure;
                        shearStressTensor_[cell].zz() += scalarPressure;

                        //- terms involving pressure tensor should not be
                        // multiplied by the number density
                        // (see Bird corrigendum)

                        heatFluxVector_[cell].x() =
                            rhoN_[cell]*
                            (
                                0.5*dsmcMccuCum_[cell]/dsmcNCum_[cell]
                              - 0.5*dsmcMccCum_[cell]/dsmcNCum_[cell]*UMean_[cell].x()
                              + dsmcEuCum_[cell]/dsmcNCum_[cell]
                              - dsmcECum_[cell]/dsmcNCum_[cell]*UMean_[cell].x()
                            )
                          - pressureTensor_[cell].xx()*UMean_[cell].x()
                          - pressureTensor_[cell].xy()*UMean_[cell].y()
                          - pressureTensor_[cell].xz()*UMean_[cell].z();

                        heatFluxVector_[cell].y() =
                            rhoN_[cell]*
                            (
                                0.5*dsmcMccvCum_[cell]/dsmcNCum_[cell]
                              - 0.5*dsmcMccCum_[cell]/dsmcNCum_[cell]*UMean_[cell].y()
                              + dsmcEvCum_[cell]/dsmcNCum_[cell]
                              - dsmcECum_[cell]/dsmcNCum_[cell]*UMean_[cell].y()
                            )
                          - pressureTensor_[cell].yx()*UMean_[cell].x()
                          - pressureTensor_[cell].yy()*UMean_[cell].y()
                          - pressureTensor_[cell].yz()*UMean_[cell].z();

                        heatFluxVector_[cell].z() =
                            rhoN_[cell]*
                            (
                                0.5*dsmcMccwCum_[cell]/dsmcNCum_[cell]
                              - 0.5*dsmcMccCum_[cell]/dsmcNCum_[cell]*UMean_[cell].z()
                              + dsmcEwCum_[cell]/dsmcNCum_[cell]
                              - dsmcECum_[cell]/dsmcNCum_[cell]*UMean_[cell].z()
                            )
                          - pressureTensor_[cell].zx()*UMean_[cell].x()
                          - pressureTensor_[cell].zy()*UMean_[cell].y()
                          - pressureTensor_[cell].zz()*UMean_[cell].z();
                    }
                    else
                    {
                        pressureTensor_[cell] = tensor::zero;
                        shearStressTensor_[cell] = tensor::zero;
                        heatFluxVector_[cell] = vector::zero;
                    }
                }

                forAll(dsmcNSpeciesCum_, i)
                {
                    const label iD = typeIds_[i];

                    if (dsmcNCum_[cell] > SMALL)
                    {
                        const scalar molarFrac = dsmcNSpeciesCum_[i][cell]
                            /dsmcNCum_[cell];

                        molecularMass[cell] += molarFrac
                            *cloud_.constProps(iD).mass();
                        molarCv_transrot[cell] += molarFrac
                           *(
                                3.0
                              + cloud_.constProps(iD)
                                  .rotationalDegreesOfFreedom()
                            );
                        molarCp_transrot[cell] += molarFrac
                           *(
                                5.0
                              + cloud_.constProps(iD)
                                  .rotationalDegreesOfFreedom()
                            );
                    }
                }

                particleCv[cell] = molarCv_transrot[cell]/NAvo;

                const scalar gasConstant =
                (
                    molecularMass[cell] > 0.0 && translationalT_[cell] > SMALL
                  ? kB/molecularMass[cell]
                  : 0.0
                );

                const scalar gamma =
                (
                    molarCv_transrot[cell] > SMALL
                  ? molarCp_transrot[cell]/molarCv_transrot[cell]
                  : 0.0
                );

                const scalar speedOfSound =
                    sqrt
                    (
                        gamma*gasConstant*translationalT_[cell]
                    );

                Ma_[cell] =
                (
                    speedOfSound > SMALL
                  ? mag(UMean_[cell])/speedOfSound
                  : 0.0
                );

                if (measureMeanFreePath_)
                {
                    forAll(typeIds_, i)
                    {
                        mfp_[i][cell] = 0.0;
                        mcr_[i][cell] = 0.0;

                        forAll(typeIds_, qspec)
                        {
                            const scalar dPQ =
                                0.5*
                                (
                                    cloud_.constProps(typeIds_[i]).d()
                                  + cloud_.constProps(typeIds_[qspec]).d()
                                );
                            const scalar omegaPQ =
                                0.5*
                                (
                                    cloud_.constProps(typeIds_[i]).omega()
                                  + cloud_.constProps(typeIds_[qspec]).omega()
                                );
                            const scalar massRatio =
                                cloud_.constProps(typeIds_[i]).mass()
                               /cloud_.constProps(typeIds_[qspec]).mass();

                            if
                            (
                                dsmcNSpeciesCum_[qspec][cell] > SMALL
                             && translationalT_[cell] > SMALL
                            )
                            {
                                const scalar nDensQ =
                                    nSpeciesCum_[qspec][cell]
                                   /(mesh_.cellVolumes()[cell]*nAvTimeSteps);
                                const scalar reducedMass =
                                    cloud_.constProps(typeIds_[i]).mass()
                                   *cloud_.constProps(typeIds_[qspec]).mass()
                                   /
                                    (
                                       cloud_.constProps(typeIds_[i]).mass()
                                     + cloud_.constProps(typeIds_[qspec]).mass()
                                    );

                                mfp_[i][cell] += pi*sqr(dPQ)*nDensQ
                                    *pow
                                    (
                                        mfpReferenceTemperature_
                                       /translationalT_[cell],
                                        omegaPQ - 0.5
                                    )*sqrt(1.0+massRatio); //Bird, eq (4.76)

                                mcr_[i][cell] += 2.0*sqrt(pi)*sqr(dPQ)
                                    *nDensQ
                                    *pow
                                    (
                                        translationalT_[cell]
                                       /mfpReferenceTemperature_,
                                        1.0-omegaPQ
                                    )
                                    *sqrt
                                    (
                                        2.0*kB*mfpReferenceTemperature_
                                       /reducedMass
                                    ); // Bird, eq (4.74)
                            }
                        }

                        if (mfp_[i][cell] > SMALL)
                        {
                            mfp_[i][cell] = 1.0/mfp_[i][cell];
                        }
                    }

                    meanCollisionSeparation_[cell] =
                    (
                        dsmcNCollsCum_[cell] > SMALL
                      ? collisionSeparation_[cell]/dsmcNCollsCum_[cell]
                      : GREAT
                    );

                    const scalar deltaT = cloud_.deltaTValue(cell);

                    if (rhoN_[cell] > SMALL)
                    {
                        //const scalar symmFactor = 2.0; // TODO (i == qspec ? 1.0 : 2.0);

                        cr_[cell] = dsmcNCollsCum_[cell]*cloud_.nParticles(cell)
                            /(rhoN_[cell]*mesh_.cellVolumes()[cell]
                                *nAvTimeSteps*deltaT);
                    }

                    meanFreePath_[cell] = 0.0;
                    meanCollisionRate_[cell] = 0.0;

                    forAll(typeIds_, i)
                    {
                        if (rhoN_[cell] > SMALL)
                        {
                            const scalar nDensP = nSpeciesCum_[i][cell]
                                /(mesh_.cellVolumes()[cell]*nAvTimeSteps);

                            meanFreePath_[cell] += mfp_[i][cell]*nDensP
                                /rhoN_[cell]; //Bird, eq (4.77)

                            meanCollisionRate_[cell] += mcr_[i][cell]*nDensP
                                /rhoN_[cell]; //Bird, eq (1.38)
                        }
                    }

                    if (meanFreePath_[cell] < SMALL)
                    {
                        meanFreePath_[cell] = GREAT;
                    }

                    if (meanCollisionRate_[cell] > SMALL)
                    {
                        meanCollisionTime_[cell] = 1.0/meanCollisionRate_[cell];
                        meanCollisionTimeTimeStepRatio_[cell] =
                            meanCollisionTime_[cell]/deltaT;
                    }
                    else
                    {
                        meanCollisionTime_[cell] = GREAT;
                        meanCollisionTimeTimeStepRatio_[cell] = GREAT;
                    }

                    if (meanFreePath_[cell] != GREAT)
                    {
                        scalar largestCellDimension = 0.0;

                        const labelList& pLabels
                        (
                            mesh_.cells()[cell].labels(mesh_.faces())
                        );
                        pointField pLocal(pLabels.size(), vector::zero);

                        forAll (pLabels, pointi)
                        {
                            pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                        }

                        scalarField dimension(2, 0.0);

                        dimension[0] = Foam::max(pLocal & vector(1,0,0))
                            - Foam::min(pLocal & vector(1,0,0));
                        dimension[1] = Foam::max(pLocal & vector(0,1,0))
                            - Foam::min(pLocal & vector(0,1,0));

                        largestCellDimension = dimension[0];

                        forAll(dimension, dim)
                        {
                            if (dimension[dim] > largestCellDimension)
                            {
                                largestCellDimension = dimension[dim];
                            }
                        }

                        mfpCellRatio_[cell] = meanFreePath_[cell]
                            /largestCellDimension;

                        SOF_[cell] =
                        (
                            meanFreePath_[cell] > SMALL
                          ? meanCollisionSeparation_[cell]/meanFreePath_[cell]
                          : 0.0
                        );
                    }
                    else
                    {
                        mfpCellRatio_[cell] = GREAT;
                        SOF_[cell] = GREAT;
                    }

                    // when no particle in cell, refines when it should not
                    // the condition should eliminates this
                    if (dsmcN_[cell] >= 4.0)
                    {
                        cellMfpRatio_[cell] = 1.0/mfpCellRatio_[cell];
                    }
                }

                if (measureClassifications_)
                {
                    if (dsmcNCum_[cell] > SMALL)
                    {
                        classIDistribution_[cell] = dsmcNClassICum_[cell]
                            /dsmcNCum_[cell];
                        classIIDistribution_[cell] = dsmcNClassIICum_[cell]
                            /dsmcNCum_[cell];
                        classIIIDistribution_[cell] = dsmcNClassIIICum_[cell]
                            /dsmcNCum_[cell];
                    }
                }

                if (measureErrors_)
                {
                    // TODO
                    if (
                           dsmcNMean_[cell] > SMALL && Ma_[cell] > SMALL
                        && gamma > SMALL && particleCv[cell] > SMALL
                       )
                    {
                        densityError_[cell] = 1.0
                            /sqrt(dsmcNMean_[cell]*nAvTimeSteps);
                        velocityError_[cell] = (1.0/sqrt(dsmcNMean_[cell]*nAvTimeSteps))*(1.0/(Ma_[cell]*sqrt(gamma)));
                        temperatureError_[cell] =
                            (1.0/sqrt(dsmcNMean_[cell]*nAvTimeSteps))
                           *sqrt(kB/particleCv[cell]);
                        pressureError_[cell] = sqrt(gamma)
                            /sqrt(dsmcNMean_[cell]*nAvTimeSteps);
                    }

                }
            }

            const label nPatches = mesh_.boundaryMesh().size();

            List<scalarField> molecularMassBF(nPatches);
            //- Heat capacity at constant volume/(0.5*kB), trans-rotational
            List<scalarField> molarCvBF_transrot(nPatches);
            //- Heat capacity at constant pressure/(0.5*kB), trans-rotational
            List<scalarField> molarCpBF_transrot(nPatches);
            List<scalarField> particleCvBF(nPatches);

            //- computing boundary measurements
            forAll(rhoNBF_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                molecularMassBF[j].setSize(patch.size(), 0.0);
                molarCvBF_transrot[j].setSize(patch.size(), 0.0);
                molarCpBF_transrot[j].setSize(patch.size(), 0.0);
                particleCvBF[j].setSize(patch.size(), 0.0);

                if (isA<wallPolyPatch>(patch))
                {
                    forAll(rhoN_.boundaryField()[j], k)
                    {
                        // Note: do not use the nParticles value that includes
                        // the RWF here. This is wrong because boundary
                        // measurements are performend during move steps. Hence
                        // radial weighting (if simulation is axisymmetric or
                        // spherical) is performed after the measurement. That
                        // is why the parcel RWFs are already included during
                        // the measurement step and we only need the FNUM value
                        // here.
                        const scalar nParticles = cloud_.coordSystem().dtModel()
                            .nParticles(j, k);

                        const scalar rhoNMean =
                            rhoNBF_[j][k]*nParticles/nAvTimeSteps;
                        const scalar rhoMMean =
                            rhoMBF_[j][k]*nParticles/nAvTimeSteps;
                        const scalar linearKEMean =
                            linearKEBF_[j][k]*nParticles/nAvTimeSteps;

                        rhoN_.boundaryFieldRef()[j][k] = rhoNMean;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMMean;

                        //- Translational energy mode and velocity
                        if (rhoMMean > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]
                                /rhoMBF_[j][k];

                            translationalT_.boundaryFieldRef()[j][k] =
                                2.0/(3.0*kB*rhoNMean)*
                                (
                                    linearKEMean - 0.5*rhoMMean*
                                    (
                                        UMean_.boundaryField()[j][k]
                                      & UMean_.boundaryField()[j][k]
                                    )
                                );
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                            translationalT_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Rotational energy mode
                        const scalar totalrDof =
                        (
                            rhoNBF_[j][k] > SMALL
                          ? rotationalDofBF_[j][k]/rhoNBF_[j][k]
                          : 0.0
                        );

                        rotationalT_.boundaryFieldRef()[j][k] =
                        (
                            rotationalDofBF_[j][k] > SMALL
                          ? 2.*rotationalEBF_[j][k]/(kB*rotationalDofBF_[j][k])
                          : 0.0
                        );

                        //- Vibrational energy mode
                        vibrationalT_.boundaryFieldRef()[j][k] = 0.0;
                        totalvDofBF_[j][k] = 0.0;

                        scalar molarFracByvDof = 0.0;

                        forAll(typeIds_, i)
                        {
                            scalarList vDofms(evmsBF_[i].size(), 0.0);
                            scalarList Tvms(vDofms.size(), 0.0);

                            vDofBF_[i][j][k] = 0.0;
                            vibTBF_[i][j][k] = 0.0;

                            if (speciesRhoNBF_[i][j][k] > 0)
                            {
                                forAll(vDofms, mode)
                                {
                                    const scalar thetaV = cloud_
                                        .constProps(typeIds_[i]).thetaV()[mode];

                                    const scalar iMean = evmsBF_[i][mode][j][k]
                                        /(kB*thetaV*speciesRhoNBF_[i][j][k]);

                                    Tvms[mode] =
                                    (
                                        iMean > 0
                                      ? thetaV/log(1.0 + 1.0/iMean)
                                      : 0.0
                                    );

                                    vDofms[mode] =
                                    (
                                      Tvms[mode] > 0.0
                                      ? 2.0*thetaV/Tvms[mode]
                                          /(exp(thetaV/Tvms[mode]) - 1.0)
                                      : 0.0
                                    );

                                    vDofBF_[i][j][k] += vDofms[mode];
                                    vibTBF_[i][j][k] += vDofms[mode]*Tvms[mode];
                                }
                            }

                            if (vDofBF_[i][j][k] > SMALL)
                            {
                                vibTBF_[i][j][k] /= vDofBF_[i][j][k];
                            }
                            else
                            {
                                vibTBF_[i][j][k] = 0.0;
                            }

                            const scalar molarFrac =
                            (
                                rhoNMean > SMALL
                              ? speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]
                              : 0.0
                            );

                            vibrationalT_.boundaryFieldRef()[j][k] +=
                                molarFrac*vDofBF_[i][j][k]*vibTBF_[i][j][k];

                            molarFracByvDof += molarFrac*vDofBF_[i][j][k];

                            totalvDofBF_[j][k] += vDofBF_[i][j][k];
                        }

                        if (totalvDofBF_[j][k] > VSMALL)
                        {
                            vibrationalT_.boundaryFieldRef()[j][k] /=
                                molarFracByvDof;
                        }
                        else
                        {
                            vibrationalT_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Electronic energy mode
                        scalar totalelDof = 0.0; // TODO
                        scalar elecT = 0.0; // TODO

                        forAll(typeIds_, i)
                        {
                            /*const label nElectronicLevels = cloud_.constProps(typeIds_[i]).nElectronicLevels();

                            if (nElectronicLevels > 1 && speciesRhoNBF_[i][j][k] > VSMALL && speciesRhoNElecBF_[j][k] > VSMALL)
                            {
                                const scalarList& electronicEnergies = cloud_.constProps(typeIds_[i]).electronicEnergyList();
                                const labelList& degeneracies = cloud_.constProps(typeIds_[i]).electronicDegeneracyList();

                                scalar speciesTransT =
                                    1.0/(3.0*kB)
                                   *(
                                        mccSpeciesBF_[i][j][k]/speciesRhoNBF_[i][j][k]
                                     - (
                                         cloud_.constProps(typeIds_[i]).mass()
                                        *sqr(mag(UMean_.boundaryField()[j][k]))
                                       )
                                   );

                                scalar fraction = speciesRhoNBF_[i][j][k]/speciesRhoNElecBF_[j][k];

                                if (speciesTransT > SMALL)
                                {
                                   scalar sum1 = 0.0;
                                   scalar sum2 = 0.0;

                                    forAll(electronicEnergies, ii)
                                    {
                                       sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(kB*speciesTransT));
                                       sum2 += degeneracies[ii]*(electronicEnergies[ii]/(kB*speciesTransT))
                                                   *exp(-electronicEnergies[ii]/(kB*speciesTransT));
                                    }

                                    if (sum2 > VSMALL && sum1> VSMALL)
                                    {
                                       scalar elecTID = electronicEBF_[i][j][k]
                                          /(kB*speciesRhoNBF_[i][j][k])*(sum1/sum2);

                                       if (elecTID > SMALL && elecTID < GREAT)
                                       {
                                           elecT += fraction*elecTID;

                                           scalar eDof = 2.0*electronicEBF_[i][j][k]
                                              /speciesRhoNBF_[i][j][k]/(kB*speciesTransT);

                                           totalelDof += fraction*eDof;
                                       }
                                    }
                                }
                            }*/
                        }

                        electronicT_.boundaryFieldRef()[j][k] = elecT; // TODO

                        overallT_.boundaryFieldRef()[j][k] =
                            (
                                (3.0*translationalT_.boundaryField()[j][k])
                              + (totalrDof*rotationalT_.boundaryField()[j][k])
                              + (
                                    totalvDofBF_[j][k]
                                   *vibrationalT_.boundaryField()[j][k]
                                )
                              + (totalelDof*elecT)
                            )
                           /(
                                3.0 + totalrDof + totalvDofBF_[j][k]
                              + totalelDof
                            );

                        //- Mach number
                        forAll(typeIds_, i)
                        {
                            const label iD = typeIds_[i];

                            if (rhoNBF_[j][k] > SMALL)
                            {
                                const scalar molarFrac = speciesRhoNBF_[i][j][k]
                                    /rhoNBF_[j][k];

                                molecularMassBF[j][k] += molarFrac
                                    *cloud_.constProps(iD).mass();

                                molarCvBF_transrot[j][k] += molarFrac
                                   *(
                                        3.0
                                      + cloud_.constProps(iD)
                                          .rotationalDegreesOfFreedom()
                                    );

                                molarCpBF_transrot[j][k] += molarFrac
                                   *(
                                        5.0
                                      + cloud_.constProps(iD)
                                          .rotationalDegreesOfFreedom()
                                    );
                            }
                        }

                        particleCvBF[j][k] = molarCvBF_transrot[j][k]/NAvo;

                        const scalar gasConstant =
                        (
                              molecularMassBF[j][k] > 0.0
                            ? kB/molecularMassBF[j][k]
                            : 0.0
                        );

                        const scalar gamma =
                        (
                            molarCvBF_transrot[j][k] > SMALL
                          ? molarCpBF_transrot[j][k]/molarCvBF_transrot[j][k]
                          : 0.0
                        );

                        const scalar speedOfSound =
                            sqrt
                            (
                                gamma*gasConstant
                               *translationalT_.boundaryField()[j][k]
                            );

                        Ma_.boundaryFieldRef()[j][k] =
                        (
                            speedOfSound > SMALL
                          ? mag(UMean_.boundaryField()[j][k])/speedOfSound
                          : 0.0
                        );

                        //- Heat flux and force density
                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;

                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;
                    }

                    forAll(p_.boundaryField()[j], facei)
                    {
                        //- Surface pressure and wall shear stress
                        p_.boundaryFieldRef()[j][facei] =
                            fD_.boundaryField()[j][facei] & n_[j][facei];

                        tau_.boundaryFieldRef()[j][facei] =
                          sqrt
                          (
                              sqr(fD_.boundaryField()[j][facei] & t1_[j][facei])
                            + sqr(fD_.boundaryField()[j][facei] & t2_[j][facei])
                          );
                    }
                }
            }

            forAll(boundaryCells_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                const bool isNonEmptyNonCyclicPatch = isA<polyPatch>(patch)
                    && !isA<emptyPolyPatch>(patch)
                    && !isA<cyclicPolyPatch>(patch);

                const bool isNonEmptyNonCyclicNonWallPatch =
                    isNonEmptyNonCyclicPatch && !isA<wallPolyPatch>(patch);

                if (isNonEmptyNonCyclicPatch)
                {
                    if (measureMeanFreePath_)
                    {
                        forAll(boundaryCells_[j], k)
                        {
                            const label celli = boundaryCells_[j][k];

                            meanFreePath_.boundaryFieldRef()[j][k] =
                                meanFreePath_[celli];
                            SOF_.boundaryFieldRef()[j][k] = SOF_[celli];
                            mfpCellRatio_.boundaryFieldRef()[j][k] =
                                mfpCellRatio_[celli];
                            meanCollisionRate_.boundaryFieldRef()[j][k] =
                                meanCollisionRate_[celli];
                            meanCollisionTime_.boundaryFieldRef()[j][k] =
                                meanCollisionTime_[celli];
                            meanCollisionTimeTimeStepRatio_
                              .boundaryFieldRef()[j][k] =
                                meanCollisionTimeTimeStepRatio_[celli];
                        }
                    }

                    if (measureHeatFluxShearStress_)
                    {
                        forAll(boundaryCells_[j], k)
                        {
                            const label celli = boundaryCells_[j][k];

                            shearStressTensor_.boundaryFieldRef()[j][k] =
                                shearStressTensor_[celli];
                            heatFluxVector_.boundaryFieldRef()[j][k] =
                                heatFluxVector_[celli];
                            pressureTensor_.boundaryFieldRef()[j][k] =
                                pressureTensor_[celli];
                        }
                    }

                    if (measureClassifications_)
                    {
                        forAll(boundaryCells_[j], k)
                        {
                            const label celli = boundaryCells_[j][k];

                            classIDistribution_.boundaryFieldRef()[j][k] =
                                classIDistribution_[celli];
                            classIIDistribution_.boundaryFieldRef()[j][k] =
                                classIIDistribution_[celli];
                            classIIIDistribution_.boundaryFieldRef()[j][k] =
                                classIIIDistribution_[celli];
                        }
                    }
                }

                if (isNonEmptyNonCyclicNonWallPatch)
                {
                    forAll(boundaryCells_[j], k)
                    {
                        const label celli = boundaryCells_[j][k];

                        translationalT_.boundaryFieldRef()[j][k] =
                            translationalT_[celli];
                        rotationalT_.boundaryFieldRef()[j][k] =
                            rotationalT_[celli];
                        vibrationalT_.boundaryFieldRef()[j][k] =
                            vibrationalT_[celli];
                        overallT_.boundaryFieldRef()[j][k] = overallT_[celli];
                        dsmcNMean_.boundaryFieldRef()[j][k] =
                            dsmcNMean_[celli];
                        rhoN_.boundaryFieldRef()[j][k] = rhoN_[celli];
                        rhoM_.boundaryFieldRef()[j][k] = rhoM_[celli];
                        p_.boundaryFieldRef()[j][k] = p_[celli];
                        Ma_.boundaryFieldRef()[j][k] = Ma_[celli];
                        UMean_.boundaryFieldRef()[j][k] = UMean_[celli];
                    }
                }
            } // end loop boundary cells for generalPatches

            if (measureMeanFreePath_)
            {
                meanFreePath_.write();
                mfpCellRatio_.write();
                //meanCollisionRate_.write();
                meanCollisionTime_.write();
                meanCollisionTimeTimeStepRatio_.write();
                //meanCollisionSeparation_.write();
                //cr_.write();
                //SOF_.write();
            }

            if (measureClassifications_)
            {
                classIDistribution_.write();
                classIIDistribution_.write();
                classIIIDistribution_.write();
            }

            if (measureErrors_)
            {
                densityError_.write();
                velocityError_.write();
                temperatureError_.write();
                pressureError_.write();
            }

            if (measureHeatFluxShearStress_)
            {
                heatFluxVector_.write();
                pressureTensor_.write();
                shearStressTensor_.write();
            }

            UMean_.write();
            p_.write();
            translationalT_.write();
            rotationalT_.write();
            vibrationalT_.write();
            electronicT_.write();
            overallT_.write();
            Ma_.write();
            q_.write();
            fD_.write();
            tau_.write();
        }

        //- reset
        if (time_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0.0;

            forAll(dsmcNCum_, celli)
            {
                dsmcNCum_[celli] = 0.0;
                dsmcNInstantaneous_[celli] = 0.0;
                dsmcMCum_[celli] = 0.0;
                dsmcLinKinEnCum_[celli] = 0.0;
                dsmcMomentumCum_[celli] = vector::zero;
                dsmcRotEnCum_[celli] = 0.0;
                dsmcRotDofCum_[celli] = 0.0;
                dsmcNWithRotDofCum_[celli] = 0.0;
                dsmcNEleCum_[celli] = 0.0,
                dsmcNClassICum_[celli] = 0.0;
                dsmcNClassIICum_[celli] = 0.0;
                dsmcNClassIIICum_[celli] = 0.0;
                collisionSeparation_[celli] = 0.0;
                dsmcNCollsCum_[celli] = 0.0;
                cr_[celli] = 0.0;
                dsmcMuuCum_[celli] = 0.0;
                dsmcMuvCum_[celli] = 0.0;
                dsmcMuwCum_[celli] = 0.0;
                dsmcMvvCum_[celli] = 0.0;
                dsmcMvwCum_[celli] = 0.0;
                dsmcMwwCum_[celli] = 0.0;
                dsmcMccCum_[celli] = 0.0;
                dsmcMccuCum_[celli] = 0.0;
                dsmcMccvCum_[celli] = 0.0;
                dsmcMccwCum_[celli] = 0.0;
                dsmcEuCum_[celli] = 0.0;
                dsmcEvCum_[celli] = 0.0;
                dsmcEwCum_[celli] = 0.0;
                dsmcECum_[celli] = 0.0;
                totalvDof_[celli] = 0.0;
                nCum_[celli] = 0.0;
                mCum_[celli] = 0.0;
                momentumCum_[celli] = vector::zero;
                linKinEnCum_[celli] = 0.0;
            }

            forAll(dsmcVibEnSpeciesModeCum_, i)
            {
                forAll(dsmcVibEnSpeciesModeCum_[i], cell)
                {
                    vibT_[i][cell] = 0.0;
                    vDof_[i][cell] = 0.0;
                }

                forAll(dsmcVibEnSpeciesModeCum_[i], mode)
                {
                    forAll(dsmcVibEnSpeciesModeCum_[i][mode], cell)
                    {
                        dsmcVibEnSpeciesModeCum_[i][mode][cell] = 0.0;
                    }
                }
            }

            forAll(dsmcNSpeciesCum_, i)
            {
                forAll(dsmcNSpeciesCum_[i], cell)
                {
                    dsmcEleEnSpeciesCum_[i][cell] = 0.0;
                    dsmcMccSpeciesCum_[i][cell] = 0.0;
                    dsmcNSpeciesCum_[i][cell] = 0.0;
                    dsmcNGroundEleLvlSpeciesCum_[i][cell] = 0.0;
                    dsmcNFirstEleLvlSpeciesCum_[i][cell] = 0.0;
                    nSpeciesCum_[i][cell] = 0.0;

                    mfp_[i][cell] = 0.0;
                    mcr_[i][cell] = 0.0;
                }
            }


            //- reset boundary information
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                speciesRhoNIntBF_[j] = 0.0;
                speciesRhoNElecBF_[j] = 0.0;
                rotationalEBF_[j] = 0.0;
                rotationalDofBF_[j] = 0.0;
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }

            forAll(vibTxvDofBF_, j)
            {
                vibTxvDofBF_[j] = 0.0;
                totalvDofBF_[j] = 0.0;
            }

            forAll(speciesRhoNBF_, i)
            {
                forAll(speciesRhoNBF_[i], j)
                {
                    speciesRhoNBF_[i][j] = 0.0;
                    vibrationalEBF_[i][j] = 0.0;
                    electronicEBF_[i][j] = 0.0;
                    mccSpeciesBF_[i][j] = 0.0;
                }
            }

            forAll(vibTBF_, i)
            {
                forAll(vibTBF_[i], j)
                {
                    vibTBF_[i][j] = 0.0;
                    vDofBF_[i][j] = 0.0;
                }

                forAll(evmsBF_[i], mode)
                {
                    forAll(evmsBF_[i][mode], j)
                    {
                        evmsBF_[i][mode][j] = 0.0;
                    }
                }
            }
        }

        if (averagingAcrossManyRuns_)
        {
            writeOut();
        }

    }
}


//- reset fields when mesh is edited
void dsmcVolFields::resetField()
{
    nTimeSteps_ = 0.0;

    dsmcNCum_.clear();
    dsmcNInstantaneous_.clear();
    dsmcMCum_.clear();
    dsmcLinKinEnCum_.clear();
    dsmcMomentumCum_.clear();
    dsmcRotEnCum_.clear();
    dsmcRotDofCum_.clear();
    dsmcNWithRotDofCum_.clear();
    dsmcNEleCum_.clear();
    dsmcNClassICum_.clear();
    dsmcNClassIICum_.clear();
    dsmcNClassIIICum_.clear();
    collisionSeparation_.clear();
    dsmcNCollsCum_.clear();
    dsmcMuuCum_.clear();
    dsmcMuvCum_.clear();
    dsmcMuwCum_.clear();
    dsmcMvvCum_.clear();
    dsmcMvwCum_.clear();
    dsmcMwwCum_.clear();
    dsmcMccCum_.clear();
    dsmcMccuCum_.clear();
    dsmcMccvCum_.clear();
    dsmcMccwCum_.clear();
    dsmcEuCum_.clear();
    dsmcEvCum_.clear();
    dsmcEwCum_.clear();
    dsmcECum_.clear();
    totalvDof_.clear();
    nCum_.clear();
    mCum_.clear();
    momentumCum_.clear();
    linKinEnCum_.clear();


    dsmcNCum_.setSize(mesh_.nCells(), 0.0);
    dsmcNInstantaneous_.setSize(mesh_.nCells(), 0.0);
    dsmcMCum_.setSize(mesh_.nCells(), 0.0);
    dsmcLinKinEnCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMomentumCum_.setSize(mesh_.nCells(), vector::zero);
    dsmcRotEnCum_.setSize(mesh_.nCells(), 0.0);
    dsmcRotDofCum_.setSize(mesh_.nCells(), 0.0);
    dsmcNWithRotDofCum_.setSize(mesh_.nCells(), 0.0);
    dsmcNEleCum_.setSize(mesh_.nCells(), 0.0);
    dsmcNClassICum_.setSize(mesh_.nCells(), 0.0);
    dsmcNClassIICum_.setSize(mesh_.nCells(), 0.0);
    dsmcNClassIIICum_.setSize(mesh_.nCells(), 0.0);
    collisionSeparation_.setSize(mesh_.nCells(), 0.0);
    dsmcNCollsCum_.setSize(mesh_.nCells(), 0.0);
    cr_.setSize(mesh_.nCells(), 0.0);
    dsmcMuuCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMuvCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMuwCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMvvCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMvwCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMwwCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMccCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMccuCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMccvCum_.setSize(mesh_.nCells(), 0.0);
    dsmcMccwCum_.setSize(mesh_.nCells(), 0.0);
    dsmcEuCum_.setSize(mesh_.nCells(), 0.0);
    dsmcEvCum_.setSize(mesh_.nCells(), 0.0);
    dsmcEwCum_.setSize(mesh_.nCells(), 0.0);
    dsmcECum_.setSize(mesh_.nCells(), 0.0);
    totalvDof_.setSize(mesh_.nCells(), 0.0);
    nCum_.setSize(mesh_.nCells(), 0.0);
    mCum_.setSize(mesh_.nCells(), 0.0);
    momentumCum_.setSize(mesh_.nCells(), vector::zero);
    linKinEnCum_.setSize(mesh_.nCells(), 0.0);

    forAll(dsmcNSpeciesCum_, i)
    {
        dsmcEleEnSpeciesCum_[i].clear();
        dsmcMccSpeciesCum_[i].clear();
        dsmcNSpeciesCum_[i].clear();
        dsmcNGroundEleLvlSpeciesCum_[i].clear();
        dsmcNFirstEleLvlSpeciesCum_[i].clear();
        nSpeciesCum_[i].clear();

        mfp_[i].clear();
        mcr_[i].clear();

        dsmcEleEnSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcMccSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcNSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcNGroundEleLvlSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        dsmcNFirstEleLvlSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);
        nSpeciesCum_[i].setSize(mesh_.nCells(), 0.0);

        mfp_[i].setSize(mesh_.nCells(), 0.0);
        mcr_[i].setSize(mesh_.nCells(), 0.0);
    }


    forAll(dsmcVibEnSpeciesModeCum_, i)
    {
        forAll(dsmcVibEnSpeciesModeCum_[i], mode)
        {
           dsmcVibEnSpeciesModeCum_[i][mode].clear();

           dsmcVibEnSpeciesModeCum_[i][mode].setSize(mesh_.nCells(), 0.0);
        }
    }

    forAll(vDof_, i)
    {
        vibT_[i].clear();
        vDof_[i].clear();

        vibT_[i].setSize(mesh_.nCells(), 0.0);
        vDof_[i].setSize(mesh_.nCells(), 0.0);
    }

    // reset boundary information
    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        rhoNBF_[j].clear();
        rhoMBF_[j].clear();
        linearKEBF_[j].clear();
        momentumBF_[j].clear();
        rotationalEBF_[j].clear();
        rotationalDofBF_[j].clear();
        qBF_[j].clear();
        fDBF_[j].clear();
        vibTxvDofBF_[j].clear();
        totalvDofBF_[j].clear();
        speciesRhoNIntBF_[j].clear();
        speciesRhoNElecBF_[j].clear();

        n_[j].clear();
        t1_[j].clear();
        t2_[j].clear();

        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        fDBF_[j].setSize(patch.size(), vector::zero);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);

        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }

    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].clear();
        electronicEBF_[i].clear();
        speciesRhoNBF_[i].clear();
        mccSpeciesBF_[i].clear();
        vibTBF_[i].clear();
        vDofBF_[i].clear();
        evmsBF_[i].clear();

        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            vibrationalEBF_[i][j].clear();
            electronicEBF_[i][j].clear();
            speciesRhoNBF_[i][j].clear();
            mccSpeciesBF_[i][j].clear();
            vibTBF_[i][j].clear();
            vDofBF_[i][j].clear();

            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }

        forAll(evmsBF_[i], mode)
        {
            evmsBF_[i][mode].setSize(mesh_.boundaryMesh().size());

            forAll(evmsBF_[i][mode], j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                evmsBF_[i][mode][j].setSize(patch.size(), 0.0);
            }
        }
    }

    forAll(boundaryCells_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        boundaryCells_[j].clear();
        boundaryCells_[j].setSize(patch.size());

        forAll(boundaryCells_[j], k)
        {
            boundaryCells_[j][k] = patch.faceCells()[k];
        }
    }
}


void dsmcVolFields::writeField()
{}


void dsmcVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);
}

} // End namespace Foam

// ************************************************************************** //

