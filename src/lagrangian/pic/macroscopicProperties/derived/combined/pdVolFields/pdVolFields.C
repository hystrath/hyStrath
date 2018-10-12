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

Measures overall temperature, including vibrational temperature, for a single species gas
or a gas mixture and writes the results to a volume scalar field that can be viewed in Paraview.

Translational, rotatational and vibrational temperature field will also be written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements class and are also written.

\*---------------------------------------------------------------------------*/

#include "pdVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdVolFields, 0);

addToRunTimeSelectionTable(pdField, pdVolFields, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdVolFields::pdVolFields
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    mfpReferenceTemperature_(273.0),
    fieldName_(propsDict_.lookup("fieldName")),
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
    pdRhoN_
    (
        IOobject
        (
            "pdRhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
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
    rhoQ_
    (
        IOobject
        (
            "rhoQ_"+ fieldName_,
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
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
    Ma_
    (
        IOobject
        (
            "Ma_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
    phiE_
    (
        IOobject
        (
            "phiE_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.0)
    ),
    phiEMean_
    (
        IOobject
        (
            "phiEMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.0)
    ),
    E_
    (
        IOobject
        (
            "E_"+ fieldName_,
            time_.time().timeName(),
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
    EMean_
    (
        IOobject
        (
            "EMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 1, -3, 0, 0, -1, 0),
            vector::zero
        )
    ),
    U_
    (
        IOobject
        (
            "U_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
    fEM_
    (
        IOobject
        (
            "fEM_"+ fieldName_,
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
    Jp_
    (
        IOobject
        (
            "Jp_"+ fieldName_,
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -3, 0, 0),
            vector::zero
        )
    ),
    maxwellTensor_
    (
        IOobject
        (
            "maxwellTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            tensor::zero
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoQMean_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    rotationalEMean_(mesh_.nCells(), 0.0),
    rotationalDofMean_(mesh_.nCells(), 0.0),
    //EMean_(mesh_.nCells(), vector::zero),
    JpMean_(mesh_.nCells(), vector::zero),
    muu_(mesh_.nCells(), 0.0),
    muv_(mesh_.nCells(), 0.0),
    muw_(mesh_.nCells(), 0.0),
    mvv_(mesh_.nCells(), 0.0),
    mvw_(mesh_.nCells(), 0.0),
    mww_(mesh_.nCells(), 0.0),
    mcc_(mesh_.nCells(), 0.0),
    mccu_(mesh_.nCells(), 0.0),
    mccv_(mesh_.nCells(), 0.0),
    mccw_(mesh_.nCells(), 0.0),
    eu_(mesh_.nCells(), 0.0),
    ev_(mesh_.nCells(), 0.0),
    ew_(mesh_.nCells(), 0.0),
    e_(mesh_.nCells(), 0.0),
    totalvDof_(mesh_.nCells(), 0.0),
    nClassI_(mesh_.nCells(), 0.0),
    nClassII_(mesh_.nCells(), 0.0),
    nClassIII_(mesh_.nCells(), 0.0),
    momentumMean_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    vibrationalETotal_(),
    nParcels_(),
    vibT_(),
    vDof_(),
    mfp_(),
    mcr_(),
    rhoNBF_(),
    rhoQBF_(),
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    qBF_(),
    vibTxvDofBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    momentumBF_(),
    EBF_(),
    fEMBF_(),
    fDBF_(),
    JpBF_(),
    vibrationalEBF_(),
    speciesRhoNBF_(),
    vibTBF_(),
    vDofBF_(),
    //EBF_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false)
{

    // standard to reading typeIds ------------
    const List<word> molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("pdVolFields::pdVolFields()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    /*******************************************************************/

    //- outer list is typeIds, inner list is number of cells on the mesh

    //- set size of outer lists
    vibT_.setSize(typeIds_.size());
    vDof_.setSize(typeIds_.size());
    vibrationalETotal_.setSize(typeIds_.size());
    nParcels_.setSize(typeIds_.size());
    mfp_.setSize(typeIds_.size());
    mcr_.setSize(typeIds_.size());

    //- set inner list size
    forAll(mfp_, i)
    {
        vibT_[i].setSize(mesh_.nCells());
        vDof_[i].setSize(mesh_.nCells());
        vibrationalETotal_[i].setSize(mesh_.nCells());
        nParcels_[i].setSize(mesh_.nCells());
        mfp_[i].setSize(mesh_.nCells());
        mcr_[i].setSize(mesh_.nCells());
    }
    /*******************************************************************/

    boundaryCells_.setSize(mesh.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];

        boundaryCells_[p].setSize(patch.size());

        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }

    //- outer list is the patch, inner list is the face

    //- initialise outer boundary lists
    //phiEBF_.setSize(mesh_.boundaryMesh().size());
    //EBF_.setSize(mesh_.boundaryMesh().size());
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoQBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    JpBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    EBF_.setSize(mesh_.boundaryMesh().size());
    fEMBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());
    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());

    //- initialise inner face list
    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        //phiEBF_[j].setSize(patch.size(), 0.0);
        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoQBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        //EBF_[j].setSize(patch.size(), vector::zero);
        JpBF_[j].setSize(patch.size(), vector::zero);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        EBF_[j].setSize(patch.size(), vector::zero);
        fEMBF_[j].setSize(patch.size(), vector::zero);
        fDBF_[j].setSize(patch.size(), vector::zero);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
    }

    //- outer list is type iD, inner list is boundary iD and inner inner list is face iD
    vibrationalEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());

    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }
    }

    if (propsDict_.found("measureClassifications"))
    {
        measureClassifications_ = Switch(propsDict_.lookup("measureClassifications"));
    }

    if (propsDict_.found("measureErrors"))
    {
        measureErrors_ = Switch(propsDict_.lookup("measureErrors"));
    }

    if(propsDict_.found("measureMeanFreePath"))
    {
        measureMeanFreePath_ = Switch(propsDict_.lookup("measureMeanFreePath"));

        if(measureMeanFreePath_)
        {
            mfpReferenceTemperature_ = readScalar(propsDict_.lookup("mfpReferenceTemperature"));
        }
    }

    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = Switch(propsDict_.lookup("averagingAcrossManyRuns"));

        // read in stored data from dictionary
        if(averagingAcrossManyRuns_)
        {
            Info << nl << "Averaging across many runs initiated." << nl << endl;

            readIn();
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdVolFields::~pdVolFields()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pdVolFields::readIn()
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

    dict.readIfPresent("rhoNMean", rhoNMean_);
    dict.readIfPresent("rhoQMean", rhoQMean_);
    dict.readIfPresent("rhoMMean", rhoMMean_);
    dict.readIfPresent("linearKEMean", linearKEMean_);
    dict.readIfPresent("momentumMean", momentumMean_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);
    dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
    dict.readIfPresent("nParcels", nParcels_);
    dict.readIfPresent("rhoNMeanInt", rhoNMeanInt_);
    dict.readIfPresent("JpMean", JpMean_);
    //dict.readIfPresent("phiEMean", phiEMean_);

    dict.readIfPresent("nTimeSteps", nTimeSteps_);
}

void pdVolFields::writeOut()
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

        dict.add("rhoNMean", rhoNMean_);
        dict.add("rhoQMean", rhoQMean_);
        dict.add("rhoMMean", rhoMMean_);
        dict.add("linearKEMean", linearKEMean_);
        dict.add("momentumMean", momentumMean_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);
        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("nParcels", nParcels_);
        dict.add("rhoNMeanInt", rhoNMeanInt_);
        dict.add("JpMean", JpMean_);

        //dict.add("phiEMean", phiEMean_);

        dict.add("nTimeSteps", nTimeSteps_);

        IOstream::streamFormat fmt = time_.time().writeFormat();
//         Pout << "fmt = " << fmt << endl;
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();

        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

//- initial conditions
void pdVolFields::createField()
{
    Info << "Initialising pdVolFields field" << endl;
}


void pdVolFields::calculateField()
{
    nTimeSteps_ += 1.0;
    const scalar eps0   = electromagnetic::epsilon0.value();    //- permitivity of space

    forAllConstIter(pdCloud, cloud_, iter)
    {
        const pdParcel& p = iter();
        label iD = findIndex(typeIds_, p.typeId());

        if(iD != -1)
        {
            //const label& cell = p.cell();
            const tetIndices tetIs(p.cell(), p.tetFace(), p.tetPt(), mesh_);

            const scalar& mass = cloud_.constProps(p.typeId()).mass();
            const scalar& rotationalDof = cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();

            const scalar& Ze = cloud_.constProps(p.typeId()).Ze();
            const scalar e = electromagnetic::e.value();   //- elementary charge

            List<List<scalar> > shape = AveragePtr_->stencil(p.position(), tetIs);

            List<scalar> shapeID = shape[0];
            List<scalar> shapeWeight = shape[1];

            forAll(shapeID,sID)
            {
                rhoNMean_[shapeID[sID]]             += 1.0                  * shapeWeight[sID];
                rhoQMean_[shapeID[sID]]             += Ze*e                 * shapeWeight[sID];
                rhoMMean_[shapeID[sID]]             += mass                 * shapeWeight[sID];
                JpMean_[shapeID[sID]]               += Ze*e*p.U()           * shapeWeight[sID];
                linearKEMean_[shapeID[sID]]         += mass*(p.U() & p.U()) * shapeWeight[sID];
                momentumMean_[shapeID[sID]]         += mass*p.U()           * shapeWeight[sID];
                rotationalEMean_[shapeID[sID]]      += p.ERot()             * shapeWeight[sID];
                rotationalDofMean_[shapeID[sID]]    += rotationalDof        * shapeWeight[sID];

                muu_[shapeID[sID]]  += mass*sqr(p.U().x())              * shapeWeight[sID];
                muv_[shapeID[sID]]  += mass*( (p.U().x()) * (p.U().y()))* shapeWeight[sID];
                muw_[shapeID[sID]]  += mass*( (p.U().x()) * (p.U().z()))* shapeWeight[sID];
                mvv_[shapeID[sID]]  += mass*sqr(p.U().y())              * shapeWeight[sID];
                mvw_[shapeID[sID]]  += mass*( (p.U().y()) * (p.U().z()))* shapeWeight[sID];
                mww_[shapeID[sID]]  += mass*sqr(p.U().z())              * shapeWeight[sID];

                mcc_[shapeID[sID]]  += mass*mag(p.U())*mag(p.U())            * shapeWeight[sID];
                mccu_[shapeID[sID]] += mass*mag(p.U())*mag(p.U())*(p.U().x())* shapeWeight[sID];
                mccv_[shapeID[sID]] += mass*mag(p.U())*mag(p.U())*(p.U().y())* shapeWeight[sID];
                mccw_[shapeID[sID]] += mass*mag(p.U())*mag(p.U())*(p.U().z())* shapeWeight[sID];

                eu_[shapeID[sID]]   += ( p.ERot() + p.EVib() )*(p.U().x())* shapeWeight[sID];
                ev_[shapeID[sID]]   += ( p.ERot() + p.EVib() )*(p.U().y())* shapeWeight[sID];
                ew_[shapeID[sID]]   += ( p.ERot() + p.EVib() )*(p.U().z())* shapeWeight[sID];
                e_[shapeID[sID]]    += ( p.ERot() + p.EVib() )            * shapeWeight[sID];

                vibrationalETotal_[iD][shapeID[sID]]    += p.EVib() *shapeWeight[sID];
                nParcels_[iD][shapeID[sID]]             += 1.0      *shapeWeight[sID];

                if(rotationalDof > VSMALL)
                {
                    rhoNMeanInt_[shapeID[sID]] += 1.0*shapeWeight[sID];
                }

                if(measureClassifications_)
                {
                    label classification = p.classification();

                    if(classification == 0)
                    {
                        nClassI_[shapeID[sID]] += 1.0*shapeWeight[sID];
                    }

                    if(classification == 1)
                    {
                        nClassII_[shapeID[sID]] += 1.0*shapeWeight[sID];
                    }

                    if(classification == 2)
                    {
                        nClassIII_[shapeID[sID]] += 1.0*shapeWeight[sID];
                    }
                }
            }
        }
    }


    phiEMean_.ref()   += cloud_.emFields().phiE();
    EMean_.ref()      += cloud_.emFields().E();

    // obtain boundary measurements
    forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
    {
        label iD = findIndex(typeIds_, i);

        if(iD != -1)
        {
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
            {
                forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
                {
                    rhoNBF_[j][k]           += cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                    rhoQBF_[j][k]           += cloud_.boundaryFluxMeasurements().rhoQBF()[i][j][k];
                    rhoMBF_[j][k]           += cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                    JpBF_[j][k]             += cloud_.boundaryFluxMeasurements().JpBF()[i][j][k];
                    linearKEBF_[j][k]       += cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                    momentumBF_[j][k]       += cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                    rotationalEBF_[j][k]    += cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                    rotationalDofBF_[j][k]  += cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                    qBF_[j][k]              += cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                    fDBF_[j][k]             += cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
                }
            }
        }
    }

    phiEMean_.boundaryFieldRef()        += cloud_.emFields().phiE_.boundaryField().boundaryInternalField();
    EMean_.boundaryFieldRef()           += cloud_.emFields().E_.boundaryField().boundaryInternalField();

    forAll(speciesRhoNBF_, i)
    {
        label iD = findIndex(typeIds_, i);

        if(iD != -1)
        {
            forAll(speciesRhoNBF_[i], j)
            {
                forAll(speciesRhoNBF_[i][j], k)
                {
                    speciesRhoNBF_[i][j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[typeIds_[i]][j][k];
                    vibrationalEBF_[i][j][k] += cloud_.boundaryFluxMeasurements().vibrationalEBF()[typeIds_[i]][j][k];
                }
            }
        }
    }

    forAll(speciesRhoNIntBF_, j)
    {
        forAll(speciesRhoNIntBF_[j], k)
        {
            speciesRhoNIntBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNIntBF()[j][k];
        }
    }

    /** If it's time to output averages **/
    if(time_.time().outputTime())
    {
        scalar nAvTimeSteps = nTimeSteps_;

        forAll(rhoNMean_, cell)
        {

            //- Average electrodynamic fields
            phiE_[cell] = phiEMean_.internalField()[cell]/nAvTimeSteps;
            E_[cell]    = EMean_.internalField()[cell]/nAvTimeSteps;

            //- Average fields if there are particles (div(0) catch)
            if(rhoNMean_[cell] > VSMALL)
            {
                scalar V = mesh_.cellVolumes()[cell];

                //- Average kinetic fields
                pdRhoN_[cell] = rhoNMean_[cell]/(nAvTimeSteps);

                rhoN_[cell] = (rhoNMean_[cell]*cloud_.nParticle())/(nAvTimeSteps*V);

                rhoQ_[cell] = (rhoQMean_[cell]*cloud_.nParticle())/(nAvTimeSteps*V);

                Jp_[cell] = (JpMean_[cell]*cloud_.nParticle())/(nAvTimeSteps*V);

                scalar rhoMMean = rhoMMean_[cell]*cloud_.nParticle()/(nAvTimeSteps*V);

                U_[cell] = momentumMean_[cell]*cloud_.nParticle() / (rhoMMean*V*nAvTimeSteps);
                scalar linearKEMean = 0.5*linearKEMean_[cell]*cloud_.nParticle()
                                        / (V*nAvTimeSteps);
                scalar rhoNMean = rhoNMean_[cell]*cloud_.nParticle()/(V*nAvTimeSteps);


                translationalT_[cell] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                *(linearKEMean - 0.5*rhoMMean*(U_[cell] & U_[cell]));

                p_[cell] = rhoN_[cell]*physicoChemical::k.value()*translationalT_[cell];
            }
            else
            {
                pdRhoN_[cell] = 0.0;
                rhoN_[cell] = 0.0;
                rhoQ_[cell] = 0.0;
                rhoM_[cell] = 0.0;
                translationalT_[cell] = 0.0;
                p_[cell] = 0.0;

                Jp_[cell] = vector::zero;
                U_[cell] = vector::zero;
            }

            if(rotationalDofMean_[cell] > VSMALL && nAvTimeSteps > VSMALL)
            {
                scalar rotationalEMean = rotationalEMean_[cell] / nAvTimeSteps;
                scalar rotationalDofMean = rotationalDofMean_[cell] / nAvTimeSteps;

                rotationalT_[cell] = (2.0/physicoChemical::k.value())*(rotationalEMean/rotationalDofMean);
            }
            else
            {
                rotationalT_[cell] = 0.0;
            }

            scalarField vibT(mesh_.nCells(), scalar(0.0));
            scalarField vibTForOverallT(mesh_.nCells(), scalar(0.0));

            forAll(vibrationalETotal_, iD)
            {
                if(vibrationalETotal_[iD][cell] > VSMALL && nParcels_[iD][cell] > VSMALL)
                {
                    scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV();

                    scalar vibrationalEMean = (vibrationalETotal_[iD][cell]/nParcels_[iD][cell]);

                    scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);

                    vibT_[iD][cell] = thetaV / log(1.0 + (1.0/iMean));

                    scalar fraction = nParcels_[iD][cell]/rhoNMeanInt_[cell];

                    vibT[cell] += vibT_[iD][cell]*fraction;

                    vDof_[iD][cell] = fraction*(2.0*thetaV/vibT_[iD][cell]) / (exp(thetaV/vibT_[iD][cell]) - 1.0);

                    totalvDof_[cell] += vDof_[iD][cell];
                }
            }

            vibrationalT_[cell] = vibT[cell];

            scalar nRotDof = 0.0;

            if(rhoNMean_[cell] > VSMALL)
            {
                nRotDof = rotationalDofMean_[cell] / rhoNMean_[cell];
            }

            overallT_[cell] = (
                                    (3.0*translationalT_[cell])
                                    + (nRotDof*rotationalT_[cell])
                                    + (totalvDof_[cell]*vibrationalT_[cell])
                                ) /
                                (3.0 + nRotDof + totalvDof_[cell]);

            if(rhoNMean_[cell] > VSMALL)
            {
                pressureTensor_[cell].xx() = rhoN_[cell]*( muu_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*U_[cell].x()*U_[cell].x()) );
                pressureTensor_[cell].xy() = rhoN_[cell]*( muv_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell])))*U_[cell].x()*U_[cell].y() );
                pressureTensor_[cell].xz() = rhoN_[cell]*( muw_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*U_[cell].x()*U_[cell].z()) );
                pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                pressureTensor_[cell].yy() = rhoN_[cell]*( mvv_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell])))*U_[cell].y()*U_[cell].y() );
                pressureTensor_[cell].yz() = rhoN_[cell]*( mvw_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*U_[cell].y()*U_[cell].z()) );
                pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                pressureTensor_[cell].zz() = rhoN_[cell]*(mww_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*U_[cell].z()*U_[cell].z()));

                scalar scalarPressure = (1.0/3.0)*(pressureTensor_[cell].xx() + pressureTensor_[cell].yy() + pressureTensor_[cell].zz());

                shearStressTensor_[cell] = -pressureTensor_[cell];
                shearStressTensor_[cell].xx() += scalarPressure;
                shearStressTensor_[cell].yy() += scalarPressure;
                shearStressTensor_[cell].zz() += scalarPressure;

                //- Calculate maxwell stress tensor field
                scalar Emag = E_[cell] & E_[cell];
                scalar Ex   = E_[cell].x();
                scalar Ey   = E_[cell].y();
                scalar Ez   = E_[cell].z();

                maxwellTensor_[cell].xx() = Ex*Ex - 0.5*Emag;
                maxwellTensor_[cell].xy() = Ex*Ey;
                maxwellTensor_[cell].xz() = Ex*Ez;
                maxwellTensor_[cell].yx() = maxwellTensor_[cell].xy();
                maxwellTensor_[cell].yy() = Ey*Ey - 0.5*Emag;
                maxwellTensor_[cell].yz() = Ey*Ez;
                maxwellTensor_[cell].zx() = maxwellTensor_[cell].xz();
                maxwellTensor_[cell].zy() = maxwellTensor_[cell].yz();
                maxwellTensor_[cell].zz() = Ez*Ez - 0.5*Emag;

                maxwellTensor_[cell] *= eps0;

                heatFluxVector_[cell].x() = rhoN_[cell]*(
                                        0.5*(mccu_[cell]/(rhoNMean_[cell]))
                                        - 0.5*(mcc_[cell]/(rhoNMean_[cell]))*U_[cell].x()
                                        + eu_[cell]/(rhoNMean_[cell])
                                        - (e_[cell]/(rhoNMean_[cell]))*U_[cell].x()
                                )
                                        - pressureTensor_[cell].xx()*U_[cell].x()
                                        - pressureTensor_[cell].xy()*U_[cell].y()
                                        - pressureTensor_[cell].xz()*U_[cell].z();

                //terms involving pressure tensor should not be multiplied by the number density (see Bird corrigendum)

                heatFluxVector_[cell].y() = rhoN_[cell]*(
                                        0.5*(mccv_[cell]/(rhoNMean_[cell]))
                                        - 0.5*(mcc_[cell]/(rhoNMean_[cell]))*U_[cell].y()
                                        + ev_[cell]/(rhoNMean_[cell])
                                        - (e_[cell]/(rhoNMean_[cell]))*U_[cell].y()
                                )
                                        - pressureTensor_[cell].yx()*U_[cell].x()
                                        - pressureTensor_[cell].yy()*U_[cell].y()
                                        - pressureTensor_[cell].yz()*U_[cell].z();

                heatFluxVector_[cell].z() = rhoN_[cell]*(
                                        0.5*(mccw_[cell]/(rhoNMean_[cell]))
                                        - 0.5*(mcc_[cell]/(rhoNMean_[cell]))*U_[cell].z()
                                        + ew_[cell]/(rhoNMean_[cell])
                                        - (e_[cell]/(rhoNMean_[cell]))*U_[cell].z()
                                )
                                        - pressureTensor_[cell].zx()*U_[cell].x()
                                        - pressureTensor_[cell].zy()*U_[cell].y()
                                        - pressureTensor_[cell].zz()*U_[cell].z();
            }
            else
            {
                pressureTensor_[cell] = tensor::zero;
                shearStressTensor_[cell] = tensor::zero;
                heatFluxVector_[cell] = vector::zero;
            }

            totalvDof_ = scalar(0.0);

            scalarField molarconstantPressureSpecificHeat(mesh_.nCells(), scalar(0.0));
            scalarField molarconstantVolumeSpecificHeat(mesh_.nCells(), scalar(0.0));
            scalarField molecularMass(mesh_.nCells(), scalar(0.0));
            scalarField particleConstantVolumeSpecificHeat(mesh_.nCells(), scalar(0.0));

            forAll(nParcels_, iD)
            {
                const label& typeId = typeIds_[iD];

                if(rhoNMean_[cell] > VSMALL)
                {
                    molecularMass[cell] += cloud_.constProps(typeId).mass()*(nParcels_[iD][cell]/rhoNMean_[cell]);
                    molarconstantPressureSpecificHeat[cell] += (5.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*(nParcels_[iD][cell]/rhoNMean_[cell]);
                    molarconstantVolumeSpecificHeat[cell] += (3.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*(nParcels_[iD][cell]/rhoNMean_[cell]);
                }
            }

            particleConstantVolumeSpecificHeat[cell] = molarconstantVolumeSpecificHeat[cell]/6.02214129e23;

            scalar gasConstant = 0.0;
            scalar gamma = 0.0;
            scalar speedOfSound = 0.0;

            if(molecularMass[cell] > VSMALL)
            {
                gasConstant = physicoChemical::k.value()/molecularMass[cell]; // R = k/m
            }

            if(molarconstantVolumeSpecificHeat[cell] > VSMALL)
            {
                gamma = molarconstantPressureSpecificHeat[cell]/molarconstantVolumeSpecificHeat[cell]; // gamma = cP/cV
            }

            if(gamma > VSMALL && gasConstant > VSMALL && translationalT_[cell] > VSMALL)
            {
                speedOfSound = sqrt(gamma*gasConstant*translationalT_[cell]);
            }

            if(speedOfSound > VSMALL)
            {
                Ma_[cell] = mag(U_[cell])/speedOfSound;
            }
            else
            {
                Ma_[cell] = 0.0;
            }

            if(measureMeanFreePath_)
            {
                forAll(mfp_, iD)
                {
                    label qspec = 0;

                    for (qspec=0; qspec<typeIds_.size(); qspec++)
                    {
                        scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());
                        scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());
                        scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();

                        if(nParcels_[qspec][cell] > VSMALL && translationalT_[cell] > VSMALL)
                        {
                            scalar nDensQ = (cloud_.nParticle()*nParcels_[qspec][cell])/(mesh_.cellVolumes()[cell]*nTimeSteps_);
                            scalar reducedMass = (cloud_.constProps(typeIds_[iD]).mass()*cloud_.constProps(typeIds_[qspec]).mass())
                                                / (cloud_.constProps(typeIds_[iD]).mass()+cloud_.constProps(typeIds_[qspec]).mass());

                            mfp_[iD][cell] += (pi*dPQ*dPQ*nDensQ*pow(mfpReferenceTemperature_/translationalT_[cell],omegaPQ-0.5)*sqrt(1.0+massRatio)); //Bird, eq (4.76)

                            mcr_[iD][cell] += (2.0*sqrt(pi)*dPQ*dPQ*nDensQ*pow(translationalT_[cell]/mfpReferenceTemperature_,1.0-omegaPQ)
                                                *sqrt(2.0*physicoChemical::k.value()*mfpReferenceTemperature_/reducedMass)); // Bird, eq (4.74)
                        }
                    }

                    if(mfp_[iD][cell] > VSMALL)
                    {
                        mfp_[iD][cell] = 1.0/mfp_[iD][cell];
                    }
                }

                meanFreePath_[cell] = 0.0;
                mfpCellRatio_[cell] = 0.0;
                meanCollisionRate_[cell] = 0.0;
                meanCollisionTime_[cell] = 0.0;
                meanCollisionTimeTimeStepRatio_[cell] = 0.0;

                forAll(mfp_, iD)
                {
                    if(rhoN_[cell] > VSMALL)
                    {
                        scalar nDensP = (cloud_.nParticle()*nParcels_[iD][cell])/(mesh_.cellVolumes()[cell]*nTimeSteps_);

                        meanFreePath_[cell] += mfp_[iD][cell]*nDensP/rhoN_[cell]; //Bird, eq (4.77)

                        meanCollisionRate_[cell] += mcr_[iD][cell]*nDensP/rhoN_[cell]; //Bird, eq (1.38)
                    }
                }

                if(meanFreePath_[cell] < VSMALL)
                {
                    meanFreePath_[cell] = GREAT;
                }

                const scalar deltaT = mesh_.time().deltaTValue();

                if(meanCollisionRate_[cell] > VSMALL)
                {
                    meanCollisionTime_[cell] = 1.0/meanCollisionRate_[cell];
                    meanCollisionTimeTimeStepRatio_[cell] = meanCollisionTime_[cell]/deltaT;
                }
                else
                {
                    meanCollisionTime_[cell] = GREAT;
                    meanCollisionTimeTimeStepRatio_[cell] = GREAT;
                }

                forAll(mfp_, iD)
                {
                    mfp_[iD] = scalar(0.0);
                    mcr_[iD] = scalar(0.0);
                }

                scalar largestCellDimension = 0.0;

                labelList pLabels(mesh_.cells()[cell].labels(mesh_.faces()));
                pointField pLocal(pLabels.size(), vector::zero);

                forAll (pLabels, pointi)
                {
                    pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                }

                scalarField dimension;

                dimension.setSize(3, 0.0);

                dimension[0] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
                dimension[1] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
                dimension[2] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));

                largestCellDimension = dimension[0];

                label dim = 0;

                for (dim=0; dim<dimension.size(); dim++)
                {
                    if(dimension[dim] > largestCellDimension)
                    {
                        largestCellDimension = dimension[dim];
                    }
                }

                mfpCellRatio_[cell] = meanFreePath_[cell]/largestCellDimension;
            }

            if(measureClassifications_)
            {
                if(rhoNMean_[cell] > VSMALL)
                {
                    classIDistribution_[cell] = nClassI_[cell]/rhoNMean_[cell];
                    classIIDistribution_[cell] = nClassII_[cell]/rhoNMean_[cell];
                    classIIIDistribution_[cell] = nClassIII_[cell]/rhoNMean_[cell];
                }
            }

            if(measureErrors_)
            {
                if(pdRhoN_[cell] > VSMALL && Ma_[cell] > VSMALL && gamma > VSMALL && particleConstantVolumeSpecificHeat[cell] > VSMALL)
                {
                    densityError_[cell] = 1.0/sqrt(pdRhoN_[cell]*nTimeSteps_);
                    velocityError_[cell] = (1.0/sqrt(pdRhoN_[cell]*nTimeSteps_))*(1.0/(Ma_[cell]*sqrt(gamma)));
                    temperatureError_[cell] = (1.0/sqrt(pdRhoN_[cell]*nTimeSteps_))
                        *sqrt(physicoChemical::k.value()/particleConstantVolumeSpecificHeat[cell]);
                    pressureError_[cell] = sqrt(gamma)/sqrt(pdRhoN_[cell]*nTimeSteps_);
                }

            }
        }

        forAll(boundaryCells_, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            if(isA<polyPatch>(patch))
            {
                if(!isA<emptyPolyPatch>(patch))
                {
                    if(!isA<cyclicPolyPatch>(patch))
                    {
                        forAll(boundaryCells_[j], k)
                        {
                            translationalT_.boundaryFieldRef()[j][k]   = translationalT_[boundaryCells_[j][k]];
                            rotationalT_.boundaryFieldRef()[j][k]      = rotationalT_[boundaryCells_[j][k]];
                            vibrationalT_.boundaryFieldRef()[j][k]     = vibrationalT_[boundaryCells_[j][k]];
                            overallT_.boundaryFieldRef()[j][k]         = overallT_[boundaryCells_[j][k]];
                            pdRhoN_.boundaryFieldRef()[j][k]           = pdRhoN_[boundaryCells_[j][k]];
                            rhoN_.boundaryFieldRef()[j][k]             = rhoN_[boundaryCells_[j][k]];
                            rhoQ_.boundaryFieldRef()[j][k]             = rhoQ_[boundaryCells_[j][k]];
                            rhoM_.boundaryFieldRef()[j][k]             = rhoM_[boundaryCells_[j][k]];
                            p_.boundaryFieldRef()[j][k]                = p_[boundaryCells_[j][k]];
                            Ma_.boundaryFieldRef()[j][k]               = Ma_[boundaryCells_[j][k]];
                            U_.boundaryFieldRef()[j][k]                = U_[boundaryCells_[j][k]];
                            Jp_.boundaryFieldRef()[j][k]               = Jp_[boundaryCells_[j][k]];
                            E_.boundaryFieldRef()[j][k]                = E_[boundaryCells_[j][k]];
                            phiE_.boundaryFieldRef()[j][k]             = phiE_[boundaryCells_[j][k]];
                        }
                    }
                }
            }
            if(measureMeanFreePath_)
            {
                if(!isA<emptyPolyPatch>(patch))
                {
                    if(!isA<cyclicPolyPatch>(patch))
                    {
                        forAll(boundaryCells_[j], k)
                        {
                            meanFreePath_.boundaryFieldRef()[j][k] = meanFreePath_[boundaryCells_[j][k]];
                            mfpCellRatio_.boundaryFieldRef()[j][k] = mfpCellRatio_[boundaryCells_[j][k]];
                            meanCollisionRate_.boundaryFieldRef()[j][k] = meanCollisionRate_[boundaryCells_[j][k]];
                            meanCollisionTime_.boundaryFieldRef()[j][k] = meanCollisionTime_[boundaryCells_[j][k]];
                            meanCollisionTimeTimeStepRatio_.boundaryFieldRef()[j][k] = meanCollisionTimeTimeStepRatio_[boundaryCells_[j][k]];
                        }
                    }
                }
            }
            if(measureClassifications_)
            {
                if(isA<polyPatch>(patch))
                {
                    if(!isA<emptyPolyPatch>(patch))
                    {
                        if(!isA<cyclicPolyPatch>(patch))
                        {
                            forAll(boundaryCells_[j], k)
                            {
                                classIDistribution_.boundaryFieldRef()[j][k] = classIDistribution_[boundaryCells_[j][k]];
                                classIIDistribution_.boundaryFieldRef()[j][k] = classIIDistribution_[boundaryCells_[j][k]];
                                classIIIDistribution_.boundaryFieldRef()[j][k] = classIIIDistribution_[boundaryCells_[j][k]];
                            }
                        }
                    }
                }
            }
        }

        if(measureMeanFreePath_)
        {
            meanFreePath_.write();
            mfpCellRatio_.write();
            meanCollisionRate_.write();
            meanCollisionTime_.write();
            meanCollisionTimeTimeStepRatio_.write();
        }

        if(measureClassifications_)
        {
            classIDistribution_.write();
            classIIDistribution_.write();
            classIIIDistribution_.write();
        }

        if(measureErrors_)
        {
            densityError_.write();
            velocityError_.write();
            temperatureError_.write();
            pressureError_.write();
        }

        pdRhoN_.boundaryFieldRef() = pdRhoN_.boundaryField().boundaryInternalField();

        phiE_.boundaryFieldRef()   = phiEMean_.boundaryField().boundaryInternalField()/nAvTimeSteps;
        E_.boundaryFieldRef()      = EMean_.boundaryField().boundaryInternalField()/nAvTimeSteps;

        /** Computing boundary measurements **/

        forAll(mesh_.boundaryMesh(), j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            //if(isA<polyPatch>(patch))
            //{
                if(!isA<emptyPolyPatch>(patch))
                {
                    forAll(mesh_.boundaryMesh()[j], k)
                    {
                        rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                        rhoQ_.boundaryFieldRef()[j][k] = rhoQBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                        Jp_.boundaryFieldRef()[j][k]   = JpBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;

                        scalar Emag = E_.boundaryField()[j][k] & E_.boundaryField()[j][k];
                        scalar Ex   = E_.boundaryField()[j][k].x();
                        scalar Ey   = E_.boundaryField()[j][k].y();
                        scalar Ez   = E_.boundaryField()[j][k].z();

                        maxwellTensor_.boundaryFieldRef()[j][k].xx() = Ex*Ex - 0.5*Emag;
                        maxwellTensor_.boundaryFieldRef()[j][k].xy() = Ex*Ey;
                        maxwellTensor_.boundaryFieldRef()[j][k].xz() = Ex*Ez;
                        maxwellTensor_.boundaryFieldRef()[j][k].yx() = maxwellTensor_.boundaryField()[j][k].xy();
                        maxwellTensor_.boundaryFieldRef()[j][k].yy() = Ey*Ey - 0.5*Emag;
                        maxwellTensor_.boundaryFieldRef()[j][k].yz() = Ey*Ez;
                        maxwellTensor_.boundaryFieldRef()[j][k].zx() = maxwellTensor_.boundaryField()[j][k].xz();
                        maxwellTensor_.boundaryFieldRef()[j][k].zy() = maxwellTensor_.boundaryField()[j][k].yz();
                        maxwellTensor_.boundaryFieldRef()[j][k].zz() = Ez*Ez - 0.5*Emag;

                        maxwellTensor_.boundaryFieldRef()[j][k] *= eps0;
                    }
                }
            //}
        }

        rhoN_.correctBoundaryConditions();
        rhoQ_.correctBoundaryConditions();
        rhoM_.correctBoundaryConditions();
        Jp_.correctBoundaryConditions();
        //phiE_.correctBoundaryConditions();
        //E_.correctBoundaryConditions();
        //maxwellTensor_.correctBoundaryConditions();

        List<scalarField> vibTBF(mesh_.boundaryMesh().size());
        List<scalarField> molecularMass(mesh_.boundaryMesh().size());
        List<scalarField> molarconstantPressureSpecificHeat(mesh_.boundaryMesh().size());
        List<scalarField> molarconstantVolumeSpecificHeat(mesh_.boundaryMesh().size());
        List<scalarField> particleConstantVolumeSpecificHeat(mesh_.boundaryMesh().size());

        forAll(vibTBF, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            vibTBF[j].setSize(patch.size(), 0.0);
            molecularMass[j].setSize(patch.size(), 0.0);
            molarconstantPressureSpecificHeat[j].setSize(patch.size(), 0.0);
            molarconstantVolumeSpecificHeat[j].setSize(patch.size(), 0.0);
            particleConstantVolumeSpecificHeat[j].setSize(patch.size(), 0.0);
        }

        // computing boundary measurements
        forAll(mesh_.boundaryMesh(), j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            if(isA<polyPatch>(patch))
            {
                if(!isA<emptyPolyPatch>(patch))
                {
                    forAll(mesh_.boundaryMesh()[j], k)
                    {
                        if(rhoM_.boundaryField()[j][k] > VSMALL)
                        {
                            U_.boundaryFieldRef()[j][k] = momentumBF_[j][k]*cloud_.nParticle()/(rhoM_.boundaryField()[j][k]*nAvTimeSteps);
                        }
                        else
                        {
                            U_.boundaryFieldRef()[j][k] = vector::zero;
                        }

                        scalar rhoMMean = rhoMBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;

                        scalar rhoQMean = rhoQBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;

                        scalar linearKEMean = linearKEBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                        scalar rhoNMean     = rhoNBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;

                        if(rhoNMean > VSMALL)
                        {
                            translationalT_.boundaryFieldRef()[j][k] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                            *(linearKEMean - 0.5*rhoMMean*(U_.boundaryField()[j][k] & U_.boundaryField()[j][k]));

                        }
                        else
                        {
                            translationalT_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        if(rotationalDofBF_[j][k] > VSMALL)
                        {
                            rotationalT_.boundaryFieldRef()[j][k] = (2.0/physicoChemical::k.value())*(rotationalEBF_[j][k]/rotationalDofBF_[j][k]);
                        }
                        else
                        {
                            rotationalT_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        forAll(vibrationalEBF_, i)
                        {
                            if(rhoNBF_[j][k] > VSMALL)
                            {
                                molecularMass[j][k] +=  cloud_.constProps(typeIds_[i]).mass()
                                            *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);

                                molarconstantPressureSpecificHeat[j][k] += (5.0 + cloud_.constProps(typeIds_[i]).rotationalDegreesOfFreedom())
                                            *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);

                                molarconstantVolumeSpecificHeat[j][k] += (3.0 + cloud_.constProps(typeIds_[i]).rotationalDegreesOfFreedom())
                                            *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
                            }

                            if(vibrationalEBF_[i][j][k] > VSMALL && speciesRhoNBF_[i][j][k] > VSMALL)
                            {
                                const scalar& thetaV = cloud_.constProps(typeIds_[i]).thetaV();

                                scalar vibrationalEMean = (vibrationalEBF_[i][j][k]/speciesRhoNBF_[i][j][k]);

                                scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);

                                vibTBF_[i][j][k] = thetaV / log(1.0 + (1.0/iMean));

                                scalar fraction = speciesRhoNBF_[i][j][k]/speciesRhoNIntBF_[j][k];

                                vDofBF_[i][j][k] = fraction*(2.0*thetaV/vibTBF_[i][j][k]) / (exp(thetaV/vibTBF_[i][j][k]) - 1.0);

                                vibTBF[j][k] += fraction*vibTBF_[i][j][k];

                                totalvDofBF_[j][k] += vDofBF_[i][j][k];
                            }

                        }

                        if(totalvDofBF_[j][k] > VSMALL)
                        {
                            vibrationalT_.boundaryFieldRef()[j][k] = vibTBF[j][k];
                        }
                        else
                        {
                            vibrationalT_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        scalar nRotDof = 0.0;

                        if(rhoNBF_[j][k] > VSMALL)
                        {
                            nRotDof = rotationalDofBF_[j][k] / rhoNBF_[j][k];
                        }

                        overallT_.boundaryFieldRef()[j][k] = (
                                                (3.0*translationalT_.boundaryField()[j][k])
                                                + (nRotDof*rotationalT_.boundaryField()[j][k])
                                                + (totalvDofBF_[j][k]*vibrationalT_.boundaryField()[j][k])
                                            ) /
                                            (3.0 + nRotDof + totalvDofBF_[j][k]);

                        totalvDofBF_[j] = scalar(0.0);

                        particleConstantVolumeSpecificHeat[j][k] = molarconstantVolumeSpecificHeat[j][k]/6.02214129e23;

                        scalar gasConstant = 0.0;
                        scalar gamma = 0.0;
                        scalar speedOfSound = 0.0;

                        if(molecularMass[j][k] > VSMALL)
                        {
                            gasConstant = physicoChemical::k.value()/molecularMass[j][k]; // R = k/m
                        }

                        if(molarconstantVolumeSpecificHeat[j][k] > VSMALL)
                        {
                            gamma = molarconstantPressureSpecificHeat[j][k]/molarconstantVolumeSpecificHeat[j][k]; // gamma = cP/cV
                        }

                        if(gamma > VSMALL && gasConstant > VSMALL && translationalT_.boundaryField()[j][k] > VSMALL)
                        {
                            speedOfSound = sqrt(gamma*gasConstant*translationalT_.boundaryField()[j][k]);
                        }

                        if(speedOfSound > VSMALL)
                        {
                            Ma_.boundaryFieldRef()[j][k] = mag(U_.boundaryField()[j][k])/speedOfSound;
                        }
                        else
                        {
                            Ma_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;

                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;

                        Jp_.boundaryFieldRef()[j][k] = JpBF_[j][k]/nAvTimeSteps;
                    }

                    fEM_.boundaryFieldRef()[j] = maxwellTensor_.boundaryField()[j] & (patch.faceAreas()/mag(patch.faceAreas()));

                    p_.boundaryFieldRef()[j] =
                        fD_.boundaryField()[j]
                    & (patch.faceAreas()/mag(patch.faceAreas()));

                    // Wall tangential unit vector. Use the direction between the
                    // face centre and the first vertex in the list
                    vectorField t1(patch.size(), vector::zero);

                    labelField faces(patch.size(),0);

                    //- loop through all faces and set the boundary faces

                    for(label f = 0; f < patch.size(); f++)
                    {
                        label globalFaceI = patch.start() + f;

                        faces[f] = globalFaceI;
                    }

                    forAll(t1, f)
                    {
                        t1[f] = patch.faceCentres()[f] - mesh_.points()[mesh_.faces()[faces[f]][0]];
                    }

                    t1 /= mag(t1);

                    // Other tangential unit vector.  Rescaling in case face is not
                    // flat and n and t1 aren't perfectly orthogonal
                    vectorField t2 = (patch.faceAreas()/mag(patch.faceAreas()))^t1;
                    t2 /= mag(t2);

                    tau_.boundaryFieldRef()[j] = sqrt( sqr(fD_.boundaryField()[j] & t1) + sqr(fD_.boundaryField()[j] & t2));
                }
            }
        }

//         U_.correctBoundaryConditions();
//         translationalT_.correctBoundaryConditions();
//         rotationalT_.correctBoundaryConditions();
//         vibrationalT_.correctBoundaryConditions();
//         overallT_.correctBoundaryConditions();
//         Ma_.correctBoundaryConditions();
//         p_.correctBoundaryConditions();
//         q_.correctBoundaryConditions();
//         fD_.correctBoundaryConditions();
//         tau_.correctBoundaryConditions();

        //- reset
        if(time_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0.0;

            phiEMean_ = dimensionedScalar("zero",  dimensionSet(1, 2, -3, 0, 0, -1, 0), 0.0);

            forAll(rhoNMean_, c)
            {
                rhoNMean_[c] = scalar(0.0);
                rhoQMean_[c] = scalar(0.0);
                rhoMMean_[c] = scalar(0.0);
                JpMean_[c] = vector::zero;
                linearKEMean_[c] = scalar(0.0);
                momentumMean_[c] = vector::zero;
                rotationalEMean_[c] = scalar(0.0);
                rotationalDofMean_[c] = scalar(0.0);
                rhoNMeanInt_[c] = scalar(0.0);
                nClassI_[c] = scalar(0.0);
                nClassII_[c] = scalar(0.0);
                nClassIII_[c] = scalar(0.0);
                muu_[c] = scalar(0.0);
                muv_[c] = scalar(0.0);
                muw_[c] = scalar(0.0);
                mvv_[c] = scalar(0.0);
                mvw_[c] = scalar(0.0);
                mww_[c] = scalar(0.0);
                mcc_[c] = scalar(0.0);
                mccu_[c] = scalar(0.0);
                mccv_[c] = scalar(0.0);
                mccw_[c] = scalar(0.0);
                eu_[c] = scalar(0.0);
                ev_[c] = scalar(0.0);
                ew_[c] = scalar(0.0);
                e_[c] = scalar(0.0);

                //phiEMean_[c]   = scalar(0.0);
                EMean_[c]      = vector::zero;
            }

            forAll(vibrationalETotal_, iD)
            {
                forAll(vibrationalETotal_[iD], cell)
                {
                    vibrationalETotal_[iD][cell] = scalar(0.0);
                    nParcels_[iD][cell] = 0.0;
                }
            }

            // reset boundary information

            forAll(rhoNBF_, j)
            {
                //phiEBF_[j] = 0.0;
                rhoNBF_[j] = 0.0;
                rhoQBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                JpBF_[j] = vector::zero;
                linearKEBF_[j] = 0.0;
                speciesRhoNIntBF_[j] = 0.0;
                rotationalEBF_[j] = 0.0;
                rotationalDofBF_[j] = 0.0;
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                fEMBF_[j] = vector::zero;
                //EBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }

            forAll(speciesRhoNBF_, i)
            {
                forAll(speciesRhoNBF_[i], j)
                {
                    speciesRhoNBF_[i][j] = 0.0;
                    vibrationalEBF_[i][j] = 0.0;
                }
            }
        }

        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
    }
}


//- write field
void pdVolFields::writeField()
{}

void pdVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //

