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
#include "fvCFD.H"
#include "OFstream.H"
#include "hybridSubrelax.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(hybridSubrelax, 0);

addToRunTimeSelectionTable(dsmcField, hybridSubrelax, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
hybridSubrelax::hybridSubrelax
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
/*
    readScalar()
*/
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    mfpReferenceTemperature_(propsDict_.lookupOrDefault<scalar>(
        "mfpReferenceTemperature", 273.0)),
    fieldName_(propsDict_.lookup("fieldName")),
    itBeforeRelaxing_(propsDict_.lookupOrDefault<label>(
        "itBeforeRelaxing", 1000)),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN_"+ fieldName_,
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
    UMean_
    (
        IOobject
        (
            "UMean_"+ fieldName_,
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
    nTimeSteps_(0),
    typeIds_(),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    rotationalEMean_(mesh_.nCells(), 0.0),
    rotationalDofMean_(mesh_.nCells(), 0.0),
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
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    qBF_(),
    vibTxvDofBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    momentumBF_(),
    fDBF_(),
    vibrationalEBF_(),
    speciesRhoNBF_(),
    vibTBF_(),
    vDofBF_(),
    //starting fields
    rhoNMean0_(mesh_.nCells(), 0.0),
    rhoNMeanInt0_(mesh_.nCells(), 0.0),
    rhoMMean0_(mesh_.nCells(), 0.0),
    linearKEMean0_(mesh_.nCells(), 0.0),
    rotationalEMean0_(mesh_.nCells(), 0.0),
    rotationalDofMean0_(mesh_.nCells(), 0.0),
    muu0_(mesh_.nCells(), 0.0),
    muv0_(mesh_.nCells(), 0.0),
    muw0_(mesh_.nCells(), 0.0),
    mvv0_(mesh_.nCells(), 0.0),
    mvw0_(mesh_.nCells(), 0.0),
    mww0_(mesh_.nCells(), 0.0),
    mcc0_(mesh_.nCells(), 0.0),
    mccu0_(mesh_.nCells(), 0.0),
    mccv0_(mesh_.nCells(), 0.0),
    mccw0_(mesh_.nCells(), 0.0),
    eu0_(mesh_.nCells(), 0.0),
    ev0_(mesh_.nCells(), 0.0),
    ew0_(mesh_.nCells(), 0.0),
    e0_(mesh_.nCells(), 0.0),
    nClassI0_(mesh_.nCells(), 0.0),
    nClassII0_(mesh_.nCells(), 0.0),
    nClassIII0_(mesh_.nCells(), 0.0),
    momentumMean0_(mesh.nCells(), vector::zero),
    vibrationalETotal0_(),
    nParcels0_(),
    //starting boundary measures
    rhoNBF0_(),
    rhoMBF0_(),
    linearKEBF0_(),
    rotationalEBF0_(),
    rotationalDofBF0_(),
    qBF0_(),
    speciesRhoNIntBF0_(),
    momentumBF0_(),
    fDBF0_(),
    vibrationalEBF0_(),
    speciesRhoNBF0_(),

    pName_(propsDict_.lookup("pName")),
    UName_(propsDict_.lookup("UName")),
    TtName_(propsDict_.lookup("TtName")),
    TrName_(propsDict_.lookup("TrName")),
    mutrName_(propsDict_.lookup("mutrName")),
    kappatrName_(propsDict_.lookup("kappatrName")),
    rhoNNames_(propsDict_.lookup("rhoNames")),
    TvNames_(propsDict_.lookup("TvNames")),
    kappavNames_(propsDict_.lookup("kappavNames")),
    nEq_(readScalar(propsDict_.lookup("nEq"))),
    thetaV_(readList<scalar>(propsDict_.lookup("thetaVs"))),
    mass_(readList<scalar>(propsDict_.lookup("masses"))),
    rotDOF_(readList<scalar>(propsDict_.lookup("rotDOFs"))),

    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    stopSamplingReset_(false)
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
            FatalErrorIn("hybridSubrelax::hybridSubrelax()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    // ---------------------------------------------------
    
        // Note; outer list is typeIds, inner list is number of cells on the mesh
    
    vibT_.setSize(typeIds_.size());
        
    forAll(vibT_, i)
    {
        vibT_[i].setSize(mesh_.nCells());
    }

    vDof_.setSize(typeIds_.size());
    
    forAll(vDof_, i)
    {
        vDof_[i].setSize(mesh_.nCells());
    }
    
    vibrationalETotal_.setSize(typeIds_.size());
    vibrationalETotal0_.setSize(typeIds_.size());
    
    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize(mesh_.nCells());
        vibrationalETotal0_[i].setSize(mesh_.nCells());
    }
    
    nParcels_.setSize(typeIds_.size());
    nParcels0_.setSize(typeIds_.size());
    
    forAll(nParcels_, i)
    {
        nParcels_[i].setSize(mesh_.nCells());
        nParcels0_[i].setSize(mesh_.nCells());
    }
    
    mfp_.setSize(typeIds_.size());
    
    forAll(mfp_, i)
    {
        mfp_[i].setSize(mesh_.nCells());
    }
    
    mcr_.setSize(typeIds_.size());
    
    forAll(mcr_, i)
    {
        mcr_[i].setSize(mesh_.nCells());
    }
    
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
        
    // initialisation
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
    //initial boundaries
    rhoNBF0_.setSize(mesh_.boundaryMesh().size());
    rhoMBF0_.setSize(mesh_.boundaryMesh().size());
    linearKEBF0_.setSize(mesh_.boundaryMesh().size());
    momentumBF0_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF0_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF0_.setSize(mesh_.boundaryMesh().size());
    qBF0_.setSize(mesh_.boundaryMesh().size());
    fDBF0_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF0_.setSize(mesh_.boundaryMesh().size());
        
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
        // initial state
        rhoNBF0_[j].setSize(patch.size(), 0.0);
        rhoMBF0_[j].setSize(patch.size(), 0.0);
        linearKEBF0_[j].setSize(patch.size(), 0.0);
        momentumBF0_[j].setSize(patch.size(), vector::zero);
        rotationalEBF0_[j].setSize(patch.size(), 0.0);
        rotationalDofBF0_[j].setSize(patch.size(), 0.0);
        qBF0_[j].setSize(patch.size(), 0.0);
        fDBF0_[j].setSize(patch.size(), vector::zero);
        speciesRhoNIntBF0_[j].setSize(patch.size(), 0.0);
    }
    
    vibrationalEBF_.setSize(typeIds_.size());
    vibrationalEBF0_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    speciesRhoNBF0_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());
    
    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        vibrationalEBF0_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF0_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());
        
        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            vibrationalEBF0_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF0_[i][j].setSize(patch.size(), 0.0);
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
    }
    
    if(measureMeanFreePath_)
    {
        mfpReferenceTemperature_ = readScalar(propsDict_.lookup("mfpReferenceTemperature"));
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

    setInitialFields(molecules);
    writeFiles();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hybridSubrelax::~hybridSubrelax()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hybridSubrelax::readIn()
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
    dict.readIfPresent("rhoMMean", rhoMMean_);
    dict.readIfPresent("linearKEMean", linearKEMean_);
    dict.readIfPresent("momentumMean", momentumMean_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);
    dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
    dict.readIfPresent("nParcels", nParcels_);     
    dict.readIfPresent("rhoNMeanInt", rhoNMeanInt_);     
    
    dict.readIfPresent("nTimeSteps", nTimeSteps_);
    
//     Info << "Some properties read in: "
//          << "mols = " << mols_[0] 
//          << ", mass = " << mass_[0]
//          << ", averagingCounter = " << averagingCounter_
//          << endl;
}

void hybridSubrelax::writeOut()
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
        dict.add("rhoMMean", rhoMMean_);
        dict.add("linearKEMean", linearKEMean_);
        dict.add("momentumMean", momentumMean_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);
        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("nParcels", nParcels_);     
        dict.add("rhoNMeanInt", rhoNMeanInt_);     
        
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

void hybridSubrelax::writeFiles()
{
    fileName outputFile1("fields.dat");
    fileName outputFile2("boundaries.dat");
    OFstream os2(mesh_.time().timeName()/outputFile2);
    OFstream os(mesh_.time().timeName()/outputFile1);

    os << "# rhoN rhoNInt rhoM linearKE rotE rotDOF muu muv muw mvv mvw mww "
        << "mcc mccu mccv mccw eu ev ew e pX pY pZ";
    forAll(typeIds_, id)
    {
        os << " vibE" << id;
    }
    forAll(typeIds_, id)
    {
        os << " nParcels" << id;
    }
    os << nl;
    forAll(mesh_.cells(), cellI)
    {
        os << rhoNMean_[cellI] << " ";
        os << rhoNMeanInt_[cellI] << " ";
        os << rhoMMean_[cellI] << " ";
        os << linearKEMean_[cellI] << " ";
        os << rotationalEMean_[cellI] << " ";
        os << rotationalDofMean_[cellI] << " ";
        os << muu_[cellI] << " ";
        os << muv_[cellI] << " ";
        os << muw_[cellI] << " ";
        os << mvv_[cellI] << " ";
        os << mvw_[cellI] << " ";
        os << mww_[cellI] << " ";
        os << mcc_[cellI] << " ";
        os << mccu_[cellI] << " ";
        os << mccv_[cellI] << " ";
        os << mccw_[cellI] << " ";
        os << eu_[cellI] << " ";
        os << ev_[cellI] << " ";
        os << ew_[cellI] << " ";
        os << e_[cellI] << " ";
        os << momentumMean_[cellI][0] << " ";
        os << momentumMean_[cellI][1] << " ";
        os << momentumMean_[cellI][2];
        forAll(typeIds_, id)
        {
            os << " " << vibrationalETotal_[id][cellI];
        }
        forAll(typeIds_, id)
        {
            os << " " << nParcels_[id][cellI];
        }
        os << nl;
    }

    os2 << "# rhoN rhoM linearKE rotE rotDOF q sRhoNInt pX pY pZ fX fY fZ";
    forAll(typeIds_, id)
    {
        os2 << " vibE" << id;
    }
    forAll(typeIds_, id)
    {
        os2 << " speciesRhoN" << id << nl;
    }
    os2 << nl;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        forAll(mesh_.boundaryMesh()[patchI], faceI)
        {
            os2 << rhoNBF_[patchI][faceI] << " ";
            os2 << rhoMBF_[patchI][faceI] << " ";
            os2 << linearKEBF_[patchI][faceI] << " ";
            os2 << rotationalEBF_[patchI][faceI] << " ";
            os2 << rotationalDofBF_[patchI][faceI] << " ";
            os2 << qBF_[patchI][faceI] << " ";
            os2 << speciesRhoNIntBF_[patchI][faceI] << " ";
            os2 << momentumBF_[patchI][faceI][0] << " ";
            os2 << momentumBF_[patchI][faceI][1] << " ";
            os2 << momentumBF_[patchI][faceI][2] << " ";
            os2 << fDBF_[patchI][faceI][0] << " ";
            os2 << fDBF_[patchI][faceI][1] << " ";
            os2 << fDBF_[patchI][faceI][2] << " ";
            forAll(typeIds_, id)
            {
                os2 << vibrationalEBF_[id][patchI][faceI] << " ";
            }
            forAll(typeIds_, id)
            {
                os2 << speciesRhoNBF_[id][patchI][faceI] << " ";
            }
            os2 << nl;
        }
    }
}

//- initial conditions
void hybridSubrelax::createField()
{
    Info << "Initialising hybridSubrelax field" << endl;
}


void hybridSubrelax::calculateField()
{  
    nTimeSteps_++;

    scalar relaxFactor = 1.0 / scalar(itBeforeRelaxing_);
    scalar factor1 = relaxFactor / (1.0 - relaxFactor);
    scalar factor2 = 1.0 - relaxFactor;

    List<label> cellsWithParticles;
    List<label> indexesUnique;

    // Add while subrelaxing
    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& p = iter();
        label iD = findIndex(typeIds_, p.typeId());

        if(iD != -1)
        {
            const label& cell = p.cell();
            const scalar& mass = cloud_.constProps(p.typeId()).mass();
            const scalar& rotationalDof = cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();

            rhoNMean_[cell] += 1.0 * factor1;
            rhoMMean_[cell] += mass * factor1;
            linearKEMean_[cell] += mass * (p.U() & p.U()) * factor1;
            momentumMean_[cell] += mass * p.U() * factor1;
            rotationalEMean_[cell] += p.ERot() * factor1;
            rotationalDofMean_[cell] += rotationalDof * factor1; 
            vibrationalETotal_[iD][cell] += p.EVib() * factor1;
            nParcels_[iD][cell] += 1.0 * factor1;
            
            muu_[cell] += mass * sqr(p.U().x()) * factor1;
            muv_[cell] += mass * ( (p.U().x()) * (p.U().y()) ) * factor1;
            muw_[cell] += mass * ( (p.U().x()) * (p.U().z()) ) * factor1;
            mvv_[cell] += mass * sqr(p.U().y()) * factor1;
            mvw_[cell] += mass * ( (p.U().y()) * (p.U().z()) ) * factor1;
            mww_[cell] += mass * sqr(p.U().z()) * factor1;
            
            mcc_[cell] += mass * mag(p.U()) * mag(p.U()) * factor1;
            mccu_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().x())
                * factor1;
            mccv_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().y())
                 * factor1;
            mccw_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().z())
                 * factor1;
            
            eu_[cell] += ( p.ERot() + p.EVib() ) * (p.U().x()) * factor1;
            ev_[cell] += ( p.ERot() + p.EVib() ) * (p.U().y()) * factor1;
            ew_[cell] += ( p.ERot() + p.EVib() ) * (p.U().z()) * factor1;
            e_[cell] += ( p.ERot() + p.EVib() ) * factor1;

            if(rotationalDof > VSMALL)
            {
                rhoNMeanInt_[cell] += 1.0 * factor1;
            }

            if(measureClassifications_)
            {
                label classification = p.classification();
                
                if(classification == 0)
                {
                    nClassI_[cell] += 1.0 * factor1;
                }
                
                if(classification == 1)
                {
                    nClassII_[cell] += 1.0 * factor1;
                }
                
                if(classification == 2)
                {
                    nClassIII_[cell] += 1.0 * factor1;
                }
            }
        }
    }

    // Correct subrelaxation
    uniqueOrder(cellsWithParticles, indexesUnique);
    forAll(indexesUnique, indexI)
    {
        label& cell = cellsWithParticles[indexesUnique[indexI]];

        rhoNMean_[cell] *= factor2;
        rhoMMean_[cell] *= factor2;
        linearKEMean_[cell] *= factor2;
        momentumMean_[cell] *= factor2;
        rotationalEMean_[cell] *= factor2;
        rotationalDofMean_[cell] *= factor2;
        muu_[cell] *= factor2;
        muv_[cell] *= factor2;
        muw_[cell] *= factor2;
        mvv_[cell] *= factor2;
        mvw_[cell] *= factor2;
        mww_[cell] *= factor2;
        mcc_[cell] *= factor2;
        mccu_[cell] *= factor2;
        mccv_[cell] *= factor2;
        mccw_[cell] *= factor2;
        eu_[cell] *= factor2;
        ev_[cell] *= factor2;
        ew_[cell] *= factor2;
        e_[cell] *= factor2;
        forAll(typeIds_, iD)
        {
            vibrationalETotal_[iD][cell] *= factor2;
            nParcels_[iD][cell] *= factor2;
        }
        rhoNMeanInt_[cell] = rhoNMeanInt_[cell] *  factor2;
        if(measureClassifications_)
        {
            nClassI_[cell] *= factor2;
            nClassII_[cell] *= factor2;
            nClassIII_[cell] *= factor2;
        }
    }
    volScalarField nD
    (
        IOobject
        (
            "nD",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField D
    (
        IOobject
        (
            "D",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField nDInt
    (
        IOobject
        (
            "nDInt",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField lKE
    (
        IOobject
        (
            "lKE",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField muu
    (
        IOobject
        (
            "muu",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField muv
    (
        IOobject
        (
            "muv",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField muw
    (
        IOobject
        (
            "muw",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField mvv
    (
        IOobject
        (
            "mvv",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField mvw
    (
        IOobject
        (
            "mvw",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField mww
    (
        IOobject
        (
            "mww",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volScalarField mcc
    (
        IOobject
        (
            "mcc",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volVectorField mccU
    (
        IOobject
        (
            "mccU",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    );
    volVectorField eU
    (
        IOobject
        (
            "eU",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    );
    volScalarField e
    (
        IOobject
        (
            "e",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    volVectorField m
    (
        IOobject
        (
            "m",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    );
    forAll(mesh_.cells(), cell)
    {
        nD[cell] = rhoNMean_[cell];
        nDInt[cell] = rhoNMeanInt_[cell];
        D[cell] = rhoMMean_[cell];
        lKE[cell] = linearKEMean_[cell];
        muu[cell] = muu_[cell];
        muv[cell] = muv_[cell];
        muw[cell] = muw_[cell];
        mvv[cell] = mvv_[cell];
        mvw[cell] = mvw_[cell];
        mww[cell] = mww_[cell];
        mcc[cell] = mcc_[cell];
        mccU[cell].x() = mccu_[cell];
        mccU[cell].y() = mccv_[cell];
        mccU[cell].z() = mccw_[cell];
        eU[cell].x() = eu_[cell];
        eU[cell].y() = ev_[cell];
        eU[cell].z() = ew_[cell];
        e[cell] = e_[cell];
        m[cell].x() = momentumMean_[cell].x();
        m[cell].y() = momentumMean_[cell].y();
        m[cell].z() = momentumMean_[cell].z();
    }

    nD.write();
    nDInt.write();
    D.write();
    lKE.write();
    muu.write();
    muv.write();
    muw.write();
    mvv.write();
    mvw.write();
    mww.write();
    mcc.write();
    mccU.write();
    eU.write();
    e.write();
    m.write();

    // obtain boundary measurements, while subrelaxing
    forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
    {
        label iD = findIndex(typeIds_, i);

        if(iD != -1)
        {            
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
            {                
                forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
                {
                    rhoNBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k]
                        * factor1;
                    rhoMBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k]
                        * factor1;
                    linearKEBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k]
                        * factor1;
                    momentumBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k]
                        * factor1;
                    rotationalEBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k]
                        * factor1;
                    rotationalDofBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k]
                        * factor1;
                    qBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().qBF()[i][j][k]
                        * factor1;
                    fDBF_[j][k]
                        += cloud_.boundaryFluxMeasurements().fDBF()[i][j][k]
                        * factor1;
                }
            }
        }
    }
    
    forAll(speciesRhoNBF_, i)
    {
        label iD = findIndex(typeIds_, i);

        if(iD != -1)
        { 
            forAll(speciesRhoNBF_[i], j)
            {                
                forAll(speciesRhoNBF_[i][j], k)
                {
                    speciesRhoNBF_[i][j][k]
                        += cloud_.boundaryFluxMeasurements().rhoNBF()[typeIds_[i]][j][k]
                        * factor1;
                    vibrationalEBF_[i][j][k]
                        += cloud_.boundaryFluxMeasurements().vibrationalEBF()[typeIds_[i]][j][k]
                        * factor1;
                }
            }
        }
    }
    
    forAll(speciesRhoNIntBF_, j)
    {
        forAll(speciesRhoNIntBF_[j], k)
        {
            speciesRhoNIntBF_[j][k]
                += cloud_.boundaryFluxMeasurements().rhoNIntBF()[j][k]
                * factor1;
        }
    }

    // Correcting subrelaxation
    forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
    {
        label iD = findIndex(typeIds_, i);
        if(iD != -1)
        {
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
            {                
                forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
                {
                    if(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k] > VSMALL)
                    {
                        rhoNBF_[j][k] *= factor2;
                        rhoMBF_[j][k] *= factor2;
                        linearKEBF_[j][k] *= factor2;
                        momentumBF_[j][k] *= factor2;
                        rotationalEBF_[j][k] *= factor2;
                        rotationalDofBF_[j][k] *= factor2;
                        qBF_[j][k] *= factor2;
                        fDBF_[j][k] *= factor2;
                    }
                }
            }
        }
    }
    forAll(speciesRhoNBF_, i)
    {
        label iD = findIndex(typeIds_, i);
        if(iD != -1)
        { 
            forAll(speciesRhoNBF_[i], j)
            {                
                forAll(speciesRhoNBF_[i][j], k)
                {
                    if(cloud_.boundaryFluxMeasurements().rhoNBF()[typeIds_[i]][j][k] > VSMALL)
                    {
                        speciesRhoNBF_[i][j][k] *= factor2;
                        vibrationalEBF_[i][j][k] *= factor2;
                    }
                }
            }
        }
    }
    forAll(speciesRhoNIntBF_, j)
    {
        forAll(speciesRhoNIntBF_[j], k)
        {
            if(cloud_.boundaryFluxMeasurements().rhoNIntBF()[j][k] > VSMALL)
            {
                speciesRhoNIntBF_[j][k] = speciesRhoNIntBF_[j][k] * factor2;
            }
        }
    }

    if(nTimeSteps_ % itBeforeRelaxing_ == 0)
    {
        scalar newFactor1 = 1.0 / (1.0 - pow(1.0 - relaxFactor, itBeforeRelaxing_));
        scalar newFactor2 = -newFactor1 * pow(1.0 - relaxFactor, itBeforeRelaxing_);

        // fields
        rhoNMean_ = newFactor1 * rhoNMean_ + newFactor2 * rhoNMean0_;
        rhoMMean_ = newFactor1 * rhoMMean_ + newFactor2 * rhoMMean0_;
        linearKEMean_ = newFactor1 * linearKEMean_ + newFactor2 * linearKEMean0_;
        rotationalEMean_ = newFactor1 * rotationalEMean_ + newFactor2
            * rotationalEMean0_;
        rotationalDofMean_ = newFactor1 * rotationalDofMean_ + newFactor2
            * rotationalDofMean0_;
        muu_ = newFactor1 * muu_ + newFactor2 * muu0_;
        muv_ = newFactor1 * muv_ + newFactor2 * muv0_;
        muw_ = newFactor1 * muw_ + newFactor2 * muw0_;
        mvv_ = newFactor1 * mvv_ + newFactor2 * mvv0_;
        mvw_ = newFactor1 * mvw_ + newFactor2 * mvw0_;
        mww_ = newFactor1 * mww_ + newFactor2 * mww0_;
        mcc_ = newFactor1 * mcc_ + newFactor2 * mcc0_;
        mccu_ = newFactor1 * mccu_ + newFactor2 * mccu0_;
        mccv_ = newFactor1 * mccv_ + newFactor2 * mccv0_;
        mccw_ = newFactor1 * mccw_ + newFactor2 * mccw0_;
        eu_ = newFactor1 * eu_ + newFactor2 * eu0_;
        ev_ = newFactor1 * ev_ + newFactor2 * ev0_;
        ew_ = newFactor1 * ew_ + newFactor2 * ew0_;
        e_ = newFactor1 * e_ + newFactor2 * e0_;
        nClassI_ = newFactor1 * nClassI_ + newFactor2 * nClassI0_;
        nClassII_ = newFactor1 * nClassII_ + newFactor2 * nClassII0_;
        nClassIII_ = newFactor1 * nClassIII_ + newFactor2 * nClassIII0_;
        momentumMean_ = newFactor1 * momentumMean_ + newFactor2 * momentumMean0_;
        forAll(typeIds_, iD)
        {
            vibrationalETotal_[iD] = newFactor1 * vibrationalETotal_[iD]
                + newFactor2 * vibrationalETotal0_[iD];
            nParcels_[iD] = newFactor1 * nParcels_[iD] + newFactor2
                * nParcels0_[iD];
        }
        forAll(rhoNMeanInt_, cell)
        {
            rhoNMeanInt_[cell] = newFactor1 * rhoNMeanInt_[cell]
                + newFactor2 * rhoNMeanInt0_[cell];
        }
        // boundaries
        forAll(rhoNBF_, j)
        {
            rhoNBF_[j] = newFactor1 * rhoNBF_[j] + newFactor2 * rhoNBF0_[j];
            rhoMBF_[j] = newFactor1 * rhoMBF_[j] + newFactor2 * rhoMBF0_[j];
            linearKEBF_[j] = newFactor1 * linearKEBF_[j] + newFactor2
                * linearKEBF0_[j];
            momentumBF_[j] = newFactor1 * momentumBF_[j] + newFactor2
                * momentumBF0_[j];
            rotationalEBF_[j] = newFactor1 * rotationalEBF_[j] + newFactor2
                * rotationalEBF0_[j];
            rotationalDofBF_[j] = newFactor1 * rotationalDofBF_[j] + newFactor2
                * rotationalDofBF0_[j];
            qBF_[j] = newFactor1 * qBF_[j] + newFactor2 * qBF0_[j];
            fDBF_[j] = newFactor1 * fDBF_[j] + newFactor2 * fDBF0_[j];
            speciesRhoNIntBF_[j] = newFactor1 * speciesRhoNIntBF_[j]
                + newFactor2 * speciesRhoNIntBF0_[j];
            forAll(typeIds_, i)
            {
                vibrationalEBF_[i][j] = newFactor1 * vibrationalEBF_[i][j]
                    + newFactor2 * vibrationalEBF0_[i][j];
                speciesRhoNBF_[i][j] = newFactor1 * speciesRhoNBF_[i][j]
                    + newFactor2 * speciesRhoNBF0_[i][j];
            }
        }

        // fields0
        rhoNMean0_ = rhoNMean_;
        rhoMMean0_ = rhoMMean_;
        rhoNMeanInt0_ = rhoNMeanInt_;
        linearKEMean0_ = linearKEMean_;
        rotationalEMean0_ = rotationalEMean_;
        rotationalDofMean0_ = rotationalDofMean_;
        muu0_ = muu_;
        muv0_ = muv_;
        muw0_ = muw_;
        mvv0_ = mvv_;
        mvw0_ = mvw_;
        mww0_ = mww_;
        mcc0_ = mcc_;
        mccu0_ = mccu_;
        mccv0_ = mccv_;
        mccw0_ = mccw_;
        eu0_ = eu_;
        ev0_ = ev_;
        ew0_ = ew_;
        e0_ = e_;
        nClassI0_ = nClassI_;
        nClassII0_ = nClassII_;
        nClassIII0_ = nClassIII_;
        momentumMean0_ = momentumMean_;
        forAll(typeIds_, iD)
        {
            vibrationalETotal0_[iD] = vibrationalETotal_[iD];
            nParcels0_[iD] = nParcels_[iD];
        }
        // boundaries0
        forAll(rhoNBF_, j)
        {
            rhoNBF0_[j] = rhoNBF_[j];
            rhoMBF0_[j] = rhoMBF_[j];
            linearKEBF0_[j] = linearKEBF_[j];
            momentumBF0_[j] = momentumBF_[j];
            rotationalEBF0_[j] = rotationalEBF_[j];
            rotationalDofBF0_[j] = rotationalDofBF_[j];
            qBF0_[j] = qBF_[j];
            fDBF0_[j] = fDBF_[j];
            speciesRhoNIntBF0_[j] = speciesRhoNIntBF_[j];
            forAll(typeIds_, i)
            {
                vibrationalEBF0_[i][j] = vibrationalEBF_[i][j];
                speciesRhoNBF0_[i][j] = speciesRhoNBF_[i][j];
            }
        }
    }

    if(time_.time().outputTime())
    {
        scalar nAvTimeSteps = 1.0;
        
        forAll(rhoNMean_, cell)
        {
            if(rhoNMean_[cell] > VSMALL)
            {
                scalar V = mesh_.cellVolumes()[cell];
                
                dsmcRhoN_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
        
                rhoN_[cell] = (rhoNMean_[cell]*cloud_.nParticle())/(nAvTimeSteps*V);
        
                rhoM_[cell] = (rhoMMean_[cell]*cloud_.nParticle())/(nAvTimeSteps*V);
                
                scalar rhoMMean = rhoMMean_[cell]*cloud_.nParticle()/(V*nAvTimeSteps);
                UMean_[cell] = momentumMean_[cell]*cloud_.nParticle() / (rhoMMean*V*nAvTimeSteps);
                scalar linearKEMean = 0.5*linearKEMean_[cell]*cloud_.nParticle() 
                                        / (V*nAvTimeSteps);
                scalar rhoNMean = rhoNMean_[cell]*cloud_.nParticle()/(V*nAvTimeSteps);

                translationalT_[cell] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                *(linearKEMean - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell]));
                                
                p_[cell] = rhoN_[cell]*physicoChemical::k.value()*translationalT_[cell];
            }
            else
            {
                dsmcRhoN_[cell] = 0.0;
                rhoN_[cell] = 0.0;
                rhoM_[cell] = 0.0;
                UMean_[cell] = vector::zero;
                translationalT_[cell] = 0.0;           
                p_[cell] = 0.0;
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
                pressureTensor_[cell].xx() = rhoN_[cell]*( muu_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()*UMean_[cell].x()) );
                pressureTensor_[cell].xy() = rhoN_[cell]*( muv_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell])))*UMean_[cell].x()*UMean_[cell].y() );
                pressureTensor_[cell].xz() = rhoN_[cell]*( muw_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()*UMean_[cell].z()) );
                pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                pressureTensor_[cell].yy() = rhoN_[cell]*( mvv_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell])))*UMean_[cell].y()*UMean_[cell].y() );
                pressureTensor_[cell].yz() = rhoN_[cell]*( mvw_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*UMean_[cell].y()*UMean_[cell].z()) );
                pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                pressureTensor_[cell].zz() = rhoN_[cell]*(mww_[cell]/(rhoNMean_[cell]) - ((rhoMMean_[cell]/(rhoNMean_[cell]))*UMean_[cell].z()*UMean_[cell].z()));
                
                scalar scalarPressure = (1.0/3.0)*(pressureTensor_[cell].xx() + pressureTensor_[cell].yy() + pressureTensor_[cell].zz());
                
                shearStressTensor_[cell] = -pressureTensor_[cell];
                shearStressTensor_[cell].xx() += scalarPressure;
                shearStressTensor_[cell].yy() += scalarPressure;
                shearStressTensor_[cell].zz() += scalarPressure;
                
                heatFluxVector_[cell].x() = rhoN_[cell]*(
                                        0.5*(mccu_[cell]/(rhoNMean_[cell]))
                                        - 0.5*(mcc_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()
                                        + eu_[cell]/(rhoNMean_[cell])
                                        - (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()
                                )
                                        - pressureTensor_[cell].xx()*UMean_[cell].x()
                                        - pressureTensor_[cell].xy()*UMean_[cell].y()
                                        - pressureTensor_[cell].xz()*UMean_[cell].z();
                                        
                //terms involving pressure tensor should not be multiplied by the number density (see Bird corrigendum)
                                        
                heatFluxVector_[cell].y() = rhoN_[cell]*(
                                        0.5*(mccv_[cell]/(rhoNMean_[cell]))
                                        - 0.5*(mcc_[cell]/(rhoNMean_[cell]))*UMean_[cell].y()
                                        + ev_[cell]/(rhoNMean_[cell])
                                        - (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].y()
                                )
                                        - pressureTensor_[cell].yx()*UMean_[cell].x()
                                        - pressureTensor_[cell].yy()*UMean_[cell].y()
                                        - pressureTensor_[cell].yz()*UMean_[cell].z();
                                        
                heatFluxVector_[cell].z() = rhoN_[cell]*(
                                        0.5*(mccw_[cell]/(rhoNMean_[cell]))
                                        - 0.5*(mcc_[cell]/(rhoNMean_[cell]))*UMean_[cell].z()
                                        + ew_[cell]/(rhoNMean_[cell])
                                        - (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].z()
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
                Ma_[cell] = mag(UMean_[cell])/speedOfSound;
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
                            scalar nDensQ = (cloud_.nParticle()*nParcels_[qspec][cell])/(mesh_.cellVolumes()[cell]*nAvTimeSteps);
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
                        scalar nDensP = (cloud_.nParticle()*nParcels_[iD][cell])/(mesh_.cellVolumes()[cell]*nAvTimeSteps);
                        
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
                if(dsmcRhoN_[cell] > VSMALL && gamma > VSMALL && Ma_[cell] > VSMALL && particleConstantVolumeSpecificHeat[cell] > VSMALL)
                {
                    densityError_[cell] = 1.0/sqrt(dsmcRhoN_[cell]*nTimeSteps_);
                    velocityError_[cell] = (1.0/sqrt(dsmcRhoN_[cell]*nTimeSteps_))*(1.0/(Ma_[cell]*sqrt(gamma)));
                    temperatureError_[cell] = (1.0/sqrt(dsmcRhoN_[cell]*nTimeSteps_))
                        *sqrt(physicoChemical::k.value()/particleConstantVolumeSpecificHeat[cell]);
                    pressureError_[cell] = sqrt(gamma)/sqrt(dsmcRhoN_[cell]*nTimeSteps_);
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
                            translationalT_.boundaryFieldRef()[j][k] = translationalT_[boundaryCells_[j][k]];
                            rotationalT_.boundaryFieldRef()[j][k] = rotationalT_[boundaryCells_[j][k]];
                            vibrationalT_.boundaryFieldRef()[j][k] = vibrationalT_[boundaryCells_[j][k]];
                            overallT_.boundaryFieldRef()[j][k] = overallT_[boundaryCells_[j][k]];
                            dsmcRhoN_.boundaryFieldRef()[j][k] = dsmcRhoN_[boundaryCells_[j][k]];
                            rhoN_.boundaryFieldRef()[j][k] = rhoN_[boundaryCells_[j][k]];
                            rhoM_.boundaryFieldRef()[j][k] = rhoM_[boundaryCells_[j][k]];
                            p_.boundaryFieldRef()[j][k] = p_[boundaryCells_[j][k]];
                            Ma_.boundaryFieldRef()[j][k] = Ma_[boundaryCells_[j][k]];
                            UMean_.boundaryFieldRef()[j][k] = UMean_[boundaryCells_[j][k]];
                        }
                    }
                }
            }
//             if(measureMeanFreePath_)
//             {
                if(!isA<emptyPolyPatch>(patch))
                {
                    if(!isA<cyclicPolyPatch>(patch))
                    {
                        forAll(boundaryCells_[j], k)
                        {
                            if(measureMeanFreePath_)
                            {
                                meanFreePath_.boundaryFieldRef()[j][k] = meanFreePath_[boundaryCells_[j][k]];
                                mfpCellRatio_.boundaryFieldRef()[j][k] = mfpCellRatio_[boundaryCells_[j][k]];
                                meanCollisionRate_.boundaryFieldRef()[j][k] = meanCollisionRate_[boundaryCells_[j][k]];
                                meanCollisionTime_.boundaryFieldRef()[j][k] = meanCollisionTime_[boundaryCells_[j][k]];
                                meanCollisionTimeTimeStepRatio_.boundaryFieldRef()[j][k] = meanCollisionTimeTimeStepRatio_[boundaryCells_[j][k]];
                            }
                            if(measureErrors_)
                            {
                                densityError_.boundaryFieldRef()[j][k] = densityError_[boundaryCells_[j][k]];
                                pressureError_.boundaryFieldRef()[j][k] = pressureError_[boundaryCells_[j][k]];
                                velocityError_.boundaryFieldRef()[j][k] = velocityError_[boundaryCells_[j][k]];
                                temperatureError_.boundaryFieldRef()[j][k] = temperatureError_[boundaryCells_[j][k]];
                            }
                            if(measureClassifications_)
                            {
                                classIDistribution_.boundaryFieldRef()[j][k] = classIDistribution_[boundaryCells_[j][k]];
                                classIIDistribution_.boundaryFieldRef()[j][k] = classIIDistribution_[boundaryCells_[j][k]];
                                classIIIDistribution_.boundaryFieldRef()[j][k] = classIIIDistribution_[boundaryCells_[j][k]];
                            }
                        }
                    }
                }
//             }
//             if(measureClassifications_)
//             {
//                 if(isA<polyPatch>(patch))
//                 {
//                     if(!isA<emptyPolyPatch>(patch))
//                     {
//                         if(!isA<cyclicPolyPatch>(patch))
//                         {
//                             forAll(boundaryCells_[j], k)
//                             {
//                                 classIDistribution_.boundaryFieldRef()[j][k] = classIDistribution_[boundaryCells_[j][k]];
//                                 classIIDistribution_.boundaryFieldRef()[j][k] = classIIDistribution_[boundaryCells_[j][k]];
//                                 classIIIDistribution_.boundaryFieldRef()[j][k] = classIIIDistribution_[boundaryCells_[j][k]];
//                             }
//                         }
//                     }
//                 }
//             }
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
        
        dsmcRhoN_.boundaryFieldRef() = dsmcRhoN_.boundaryField().boundaryInternalField();
        
        // computing boundary measurements
        forAll(rhoNBF_, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            if(isA<wallPolyPatch>(patch))
            {
                forAll(rhoN_.boundaryField()[j], k)
                {
                    rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                    rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                }
            }
        }
        
        rhoN_.correctBoundaryConditions();
        rhoM_.correctBoundaryConditions();
        
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
        forAll(rhoNBF_, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            if(isA<wallPolyPatch>(patch))
            {
                forAll(rhoN_.boundaryField()[j], k)
                {
                    if(rhoM_.boundaryField()[j][k] > VSMALL)
                    {
                        UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]*cloud_.nParticle()/(rhoM_.boundaryField()[j][k]*nAvTimeSteps);
                    }
                    else
                    {
                        UMean_.boundaryFieldRef()[j][k] = vector::zero;
                    }
                        
                    scalar rhoMMean = rhoMBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                    scalar linearKEMean = linearKEBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                    scalar rhoNMean = rhoNBF_[j][k]*cloud_.nParticle()/nAvTimeSteps;
                    
                    if(rhoNMean > VSMALL)
                    {
                        translationalT_.boundaryFieldRef()[j][k] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                        *(linearKEMean - 0.5*rhoMMean*(UMean_.boundaryField()[j][k] & UMean_.boundaryField()[j][k]));
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
                    
                    /**************************************************************************************************************/
                    
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
                    
                    /**************************************************************************************************************/
                    
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
                        Ma_.boundaryFieldRef()[j][k] = mag(UMean_.boundaryField()[j][k])/speedOfSound;
                    }
                    else
                    {
                        Ma_.boundaryFieldRef()[j][k] = 0.0;
                    }
                    
                    q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;
                    
                    fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;
                    
                    const fvMesh& mesh = fD_.mesh();

                    forAll(mesh.boundaryMesh(), i)
                    {
                        const polyPatch& patch = mesh_.boundaryMesh()[i];
                        const vectorField& fC = patch.faceCentres();
                        
                        if (isA<wallPolyPatch>(patch))
                        {
                            p_.boundaryFieldRef()[i] =
                                fD_.boundaryField()[i]
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
                                t1[f] = fC[f] - mesh_.points()[mesh_.faces()[faces[f]][0]];
                            }
                            
                            t1 /= mag(t1);

                            // Other tangential unit vector.  Rescaling in case face is not
                            // flat and n and t1 aren't perfectly orthogonal
                            vectorField t2 = (patch.faceAreas()/mag(patch.faceAreas()))^t1;
                            t2 /= mag(t2);
                            
                            tau_.boundaryFieldRef()[i] = sqrt(
                                sqr(fD_.boundaryField()[i] & t1)
                                + sqr(fD_.boundaryField()[i] & t2));
                        }
                    }
                }
            }
        }
        
//         UMean_.correctBoundaryConditions();
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
        if(time_.resetFieldsAtOutput() and not stopSamplingReset_)
        {
            nTimeSteps_ = 0;
            
            forAll(rhoNMean_, c)
            {
                rhoNMean_[c] = scalar(0.0);
                rhoMMean_[c] = scalar(0.0);
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
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                speciesRhoNIntBF_[j] = 0.0;
                rotationalEBF_[j] = 0.0;
                rotationalDofBF_[j] = 0.0;
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
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

void hybridSubrelax::setReset
(
    bool& stopReset
)
{
    stopSamplingReset_ = stopReset;
}

void hybridSubrelax::setInitialFields
(
    const List<word>& molecules
)
{
    dimensionedScalar massDim("massDim", dimMass, 1.0);
    dimensionedScalar tempDim("tempDim", dimTemperature, 1.0);

    volScalarField p0_
    (
        IOobject
        (
            pName_,
            time_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField Tt0_
    (
        IOobject
        (
            TtName_,
            time_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField Tr0_
    (
        IOobject
        (
            TrName_,
            time_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField mutr0_
    (
        IOobject
        (
            mutrName_,
            time_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField kappatr0_
    (
        IOobject
        (
            kappatrName_,
            time_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volVectorField U0_
    (
        IOobject
        (
            UName_,
            time_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    PtrList<volScalarField> rhoNs0_(typeIds_.size());
    PtrList<volScalarField> Tvs0_(typeIds_.size());
    PtrList<volScalarField> kappavs0_(typeIds_.size());

    forAll(typeIds_, i)
    {
        rhoNs0_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    rhoNNames_[i],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        Tvs0_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    TvNames_[i],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        kappavs0_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    kappavNames_[i],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    volScalarField rhoNTot
    (
        IOobject
        (
            "rhoNTot",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.0 * rhoNs0_[0]
    );
    forAll(typeIds_, i)
    {
        rhoNTot += rhoNs0_[i];
    }

    volScalarField massEq
    (
        IOobject
        (
            "massEq",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    );
    forAll(typeIds_, id)
    {
        massEq += rhoNs0_[id] * mass_[id] * massDim / rhoNTot;
    }


    volTensorField stressTensor
    (
        IOobject
        (
            "stressTensor",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mutr0_ * (fvc::grad(U0_) + dev2(Foam::T(fvc::grad(U0_))))
    );
    stressTensor.write();
    volScalarField rotDOF
    (
        IOobject
        (
            "rotDOF",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(typeIds_, id)
    {
        rotDOF += rotDOF_[id] * rhoNs0_[id] / rhoNTot;
    }

    volVectorField qt
    (
        IOobject
        (
            "heatFluxTra",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -(15.0 / (15.0 + 2.0 * rotDOF)) * (kappatr0_ / p0_)
            * fvc::grad(Tt0_) * sqrt(2.0 * massEq
            / (physicoChemical::k * Tt0_))
    );
    qt.write();
    volVectorField qr
    (
        IOobject
        (
            "heatFluxRot",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -(2.0 * rotDOF / (15.0 + 2.0 * rotDOF)) * (kappatr0_ / p0_)
            * fvc::grad(Tt0_) * sqrt(2.0 * massEq
            / (physicoChemical::k * Tt0_))
    );
    qr.write();
    PtrList<volVectorField> qv(typeIds_.size());
    forAll(typeIds_, id)
    {
        dimensionedScalar thetaV = max(thetaV_[id], 1.0) * tempDim;
        qv.set
        (
            id,
            new volVectorField
            (
                IOobject
                (
                    "qv0_" + molecules[id],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -(kappavs0_[id] / p0_) * fvc::grad(Tvs0_[id])
                    * sqrt(2.0 * mass_[id] * massDim / (physicoChemical::k
                    * Tt0_)) * (Tt0_ * Tvs0_[id] / pow(thetaV, 2))
            )
        );

        forAll(mesh_.cells(), cell)
        {
            if(Tvs0_[id][cell] > VSMALL)
            {
                qv[id][cell] *= ((exp(thetaV.value() / Tvs0_[id][cell]) - 1.0)
                    * (exp(thetaV.value() / Tvs0_[id][cell]) - 1.0)
                    / exp(thetaV.value() / Tvs0_[id][cell]));
            }
            else
            {
                qv[id][cell] *= 0.0;
            }
        }
        forAll(Tt0_.boundaryField(), boundaryI)
        {
            forAll(Tt0_.boundaryField()[boundaryI], faceI)
            {
                if(Tvs0_[id].boundaryField()[boundaryI][faceI] > VSMALL)
                {
                    qv[id][boundaryI][faceI] *= ((exp(thetaV.value()
                        / Tvs0_[id].boundaryField()[boundaryI][faceI]) - 1.0)
                        * (exp(thetaV.value()
                        / Tvs0_[id].boundaryField()[boundaryI][faceI]) - 1.0)
                        / exp(thetaV.value()
                        / Tvs0_[id].boundaryField()[boundaryI][faceI]));
                }
                else
                {
                    qv[id][boundaryI][faceI] *= 0.0;
                }
            }
        }
        qv[id].write();
    }

    // set scalar fields
    rhoNMean0_ = rhoNTot.internalField() * mesh_.V() / nEq_;
    rhoMMean0_ = rhoNMean0_ * massEq.internalField();
    momentumMean0_ = rhoMMean0_ * U0_.internalField();
    linearKEMean0_ = 3.0 * physicoChemical::k.value() * (Tt0_.internalField()
        * rhoNMean0_) + rhoMMean0_ * (U0_.internalField()
        & U0_.internalField());
    rotationalDofMean0_ = rotDOF.internalField();
    rotationalEMean0_ = 0.5 * physicoChemical::k.value() * rotationalDofMean0_
        * Tr0_.internalField();
    e0_ = rotationalEMean0_;
    forAll(typeIds_, id)
    {
        if(rotDOF_[id] > VSMALL)
        {
            rhoNMeanInt0_ += rhoNs0_[id].internalField();
        }
        nParcels0_[id] = rhoNs0_[id].internalField() * mesh_.V()
            / nEq_;
        vibrationalETotal0_[id] = nParcels0_[id] * physicoChemical::k.value()
            * thetaV_[id];
        forAll(vibrationalETotal0_[id], cell)
        {
            if(Tvs0_[id][cell] > VSMALL)
            {
                vibrationalETotal0_[id][cell] /= (exp(thetaV_[id]
                    / Tvs0_[id][cell]) - 1.0);
            }
            else
            {
                vibrationalETotal0_[id][cell] *= 0.0;
            }
        }

        e0_ += vibrationalETotal0_[id];
    }
    rhoNMeanInt0_ *= mesh_.V() / nEq_;
    mcc0_ = 3.0 * p0_.internalField() * mesh_.V() / rhoNMean0_ / nEq_;
    forAll(muu_, cell)
    {
        muu0_[cell] = mesh_.V()[cell] * (p0_[cell] - stressTensor[cell].xx())
            / nEq_ + rhoMMean0_[cell] * U0_[cell].x()
            * U0_[cell].x();
        muv0_[cell] = -mesh_.V()[cell] * stressTensor[cell].xy()
            / nEq_ + rhoMMean0_[cell] * U0_[cell].x()
            * U0_[cell].y();
        muw0_[cell] = -mesh_.V()[cell] * stressTensor[cell].xz()
            / nEq_ + rhoMMean0_[cell] * U0_[cell].x()
            * U0_[cell].z();
        mvv0_[cell] = mesh_.V()[cell] * (p0_[cell] - stressTensor[cell].yy())
            / nEq_ + rhoMMean0_[cell] * U0_[cell].y()
            * U0_[cell].y();
        mvw0_[cell] = -mesh_.V()[cell] * stressTensor[cell].yz()
            / nEq_ + rhoMMean0_[cell] * U0_[cell].y()
            * U0_[cell].z();
        mww0_[cell] = mesh_.V()[cell] * (p0_[cell] - stressTensor[cell].zz())
            / nEq_ + rhoMMean0_[cell] * U0_[cell].z()
            * U0_[cell].z();
        mccu0_[cell] = mcc0_[cell] * U0_[cell].x() + 2.0 * mesh_.V()[cell]
            * (qt[cell].x() + (p0_[cell] - stressTensor[cell].xx())
            * U0_[cell].x() - stressTensor[cell].xy() * U0_[cell].y()
            - stressTensor[cell].xz() * U0_[cell].z()) / nEq_;
        mccv0_[cell] = mcc0_[cell] * U0_[cell].y() + 2.0 * mesh_.V()[cell]
            * (qt[cell].y() + (p0_[cell] - stressTensor[cell].yy())
            * U0_[cell].y() - stressTensor[cell].xy() * U0_[cell].x()
            - stressTensor[cell].yz() * U0_[cell].z()) / nEq_;
        mccw0_[cell] = mcc0_[cell] * U0_[cell].z() + 2.0 * mesh_.V()[cell]
            * (qt[cell].z() + (p0_[cell] - stressTensor[cell].zz())
            * U0_[cell].z() - stressTensor[cell].xz() * U0_[cell].x()
            - stressTensor[cell].yz() * U0_[cell].y()) / nEq_;
        forAll(typeIds_, id)
        {
            eu0_[cell] += qv[id][cell].x();
            ev0_[cell] += qv[id][cell].y();
            ew0_[cell] += qv[id][cell].z();
        }
        eu0_[cell] += qr[cell].x();
        ev0_[cell] += qr[cell].y();
        ew0_[cell] += qr[cell].z();
        eu0_[cell] *= (mesh_.V()[cell] / nEq_);
        ev0_[cell] *= (mesh_.V()[cell] / nEq_);
        ew0_[cell] *= (mesh_.V()[cell] / nEq_);
        eu0_[cell] += e0_[cell] * U0_[cell].x();
        ev0_[cell] += e0_[cell] * U0_[cell].y();
        ew0_[cell] += e0_[cell] * U0_[cell].z();
    }
    rhoNMean_ = rhoNMean0_;
    rhoNMeanInt_ = rhoNMeanInt0_;
    rhoMMean_ = rhoMMean0_;
    linearKEMean_ = linearKEMean0_;
    rotationalDofMean_ = rotationalDofMean0_;
    rotationalEMean_ = rotationalEMean0_;
    muu_ = muu0_;
    muv_ = muv0_;
    muw_ = muw0_;
    mvv_ = mvv0_;
    mvw_ = mvw0_;
    mww_ = mww0_;
    mcc_ = mcc0_;
    mccu_ = mccu0_;
    mccv_ = mccv0_;
    mccw_ = mccw0_;
    eu_ = eu0_;
    ev_ = ev0_;
    ew_ = ew0_;
    e_ = e0_;
    momentumMean_ = momentumMean0_;
    forAll(typeIds_, id)
    {
        nParcels_[id] = nParcels0_[id];
        vibrationalETotal_[id] = vibrationalETotal0_[id];
    }

    //set boundary fields
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        if(isA<polyPatch>(patch))
        {
            if(!isA<emptyPolyPatch>(patch))
            {
                if(!isA<cyclicPolyPatch>(patch))
                {
                    rhoNBF0_[patchI] = rhoNTot.boundaryField()[patchI]
                        / nEq_;
                    rhoMBF0_[patchI] = rhoNBF0_[patchI]
                        * massEq.boundaryField()[patchI];
                    linearKEBF0_[patchI] = 3.0
                        * physicoChemical::k.value()
                        * (Tt0_.boundaryField()[patchI]
                        * rhoNBF0_[patchI]) + rhoMBF0_[patchI]
                        * (U0_.boundaryField()[patchI]
                        & U0_.boundaryField()[patchI]);
                    momentumBF0_[patchI] = massEq.boundaryField()[patchI]
                        * U0_.boundaryField()[patchI];
                    rotationalDofBF0_[patchI] = rotDOF.boundaryField()[patchI];
                    rotationalEBF0_[patchI] = 0.5 * physicoChemical::k.value()
                        * rotationalDofBF0_[patchI]
                        * Tr0_.boundaryField()[patchI];
                    vectorField hf = qt.boundaryField()[patchI]
                        + qr.boundaryField()[patchI];
                    forAll(typeIds_, id)
                    {
                        dimensionedScalar thetaV("thetaV", dimTemperature,
                            max(1.0, thetaV_[id]));
                        hf += qv[id].boundaryField()[patchI];
                        if(rotDOF_[id] > VSMALL)
                        {
                            speciesRhoNIntBF0_[patchI]
                                += rhoNs0_[id].boundaryField()[patchI];
                        }
                        speciesRhoNBF0_[id][patchI]
                            = rhoNs0_[id].boundaryField()[patchI] / nEq_;
                        vibrationalEBF0_[id][patchI]
                            = speciesRhoNBF0_[id][patchI]
                            * physicoChemical::k.value() * thetaV.value();
                        forAll(mesh_.boundaryMesh()[patchI], faceI)
                        {
                            if(Tvs0_[id].boundaryField()[patchI][faceI]
                                > VSMALL)
                            {
                                vibrationalEBF0_[id][patchI][faceI]
                                    /= (exp(thetaV.value()  
                                    / Tvs0_[id].boundaryField()[patchI][faceI])
                                    - 1.0);
                            }
                            else
                            {
                                vibrationalEBF0_[id][patchI][faceI] *= 0.0;
                            }
                        }
                    }
                    qBF0_[patchI] = (hf
                            & mesh_.Sf().boundaryField()[patchI])
                            / mesh_.magSf().boundaryField()[patchI];
                    speciesRhoNIntBF0_[patchI] /= cloud_.nParticle();
                    forAll(mesh_.boundaryMesh()[patchI], faceI)
                    {
                        fDBF0_[patchI][faceI]
                            = (stressTensor.boundaryField()[patchI][faceI]
                            & mesh_.Sf().boundaryField()[patchI][faceI])
                            / mesh_.magSf().boundaryField()[patchI][faceI];
                    }
                }
            }
        }
    }
    rhoNBF_ = rhoNBF0_;
    rhoMBF_ = rhoMBF0_;
    linearKEBF_ = linearKEBF0_;
    rotationalEBF_ = rotationalEBF0_;
    rotationalDofBF_ = rotationalDofBF0_;
    qBF_ = qBF0_;
    speciesRhoNIntBF_ = speciesRhoNIntBF0_;
    momentumBF_ = momentumBF0_;
    fDBF_ = fDBF0_;
    vibrationalEBF_ = vibrationalEBF0_;
    speciesRhoNBF_ = speciesRhoNBF0_;
}


//- write field
void hybridSubrelax::writeField()
{}

void hybridSubrelax::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //

