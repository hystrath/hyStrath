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
    dsmcRhoNMean_
    (
        IOobject
        (
            "dsmcRhoNMean_"+ fieldName_,
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
        dimensionedScalar("zero",  dimensionSet(0, -3, -1, 0, 0), 0.0)
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
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNInstantaneous_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    molsElec_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
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
    collisionSeparation_(mesh_.nCells(), 0.0),
    nColls_(mesh_.nCells(), 0.0),
    momentumMean_(mesh.nCells(), vector::zero),
    momentumMeanXnParticle_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    vibT_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
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
    nGroundElectronicLevel_.setSize(typeIds_.size());
    nFirstElectronicLevel_.setSize(typeIds_.size());
    vDof_.setSize(typeIds_.size());
    vibrationalETotal_.setSize(typeIds_.size());
    electronicETotal_.setSize(typeIds_.size());
    nParcels_.setSize(typeIds_.size());
    nParcelsXnParticle_.setSize(typeIds_.size());
    mccSpecies_.setSize(typeIds_.size());
    mfp_.setSize(typeIds_.size());
    mcr_.setSize(typeIds_.size());
    
    boundaryCells_.setSize(mesh_.boundaryMesh().size());
        
    forAll(typeIds_, i)
    {
        vibT_[i].setSize(mesh_.nCells());
        nGroundElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
        nFirstElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
        vDof_[i].setSize(mesh_.nCells(), 0.0);
        electronicETotal_[i].setSize(mesh_.nCells(), 0.0);
        nParcels_[i].setSize(mesh_.nCells(), 0.0);
        nParcelsXnParticle_[i].setSize(mesh_.nCells(), 0.0);
        mccSpecies_[i].setSize(mesh_.nCells(), 0.0);
        mfp_[i].setSize(mesh_.nCells(), 0.0);
        mcr_[i].setSize(mesh_.nCells(), 0.0);

        vibrationalETotal_[i].setSize
        (
            cloud_.constProps(typeIds_[i]).vibrationalDegreesOfFreedom()
        );
        
        forAll(vibrationalETotal_[i], j)
        {
            vibrationalETotal_[i][j].setSize(mesh_.nCells(), 0.0);
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
    
    rhoNInstantaneous_ = 0.0;
    
    const scalar kB = physicoChemical::k.value();
    const scalar NAvo = physicoChemical::NA.value();
    
    if (sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;
        
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

                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    
                    rhoNMeanXnParticle_[cell] += nParticles;
                    rhoMMeanXnParticle_[cell] += mass*nParticles;
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
                    const label cell = p.cell();
                    const scalar nParticles = cloud_.nParticles(cell);
                    const scalar mass = cloud_.constProps(p.typeId()).mass();
                    const scalar massBySqMagU = mass*(p.U() & p.U());
                    const scalar rotationalDof = cloud_.
                        constProps(p.typeId()).rotationalDegreesOfFreedom();
                    const scalarList& electronicEnergies = cloud_.
                        constProps(p.typeId()).electronicEnergyList();

                    scalarList EVib
                    (
                        cloud_.constProps(p.typeId()).thetaV().size(), 0.0
                    );
                    
                    scalar vibEn = 0.0;
                    
                    forAll(EVib, mode)
                    {
                        EVib[mode] = p.vibLevel()[mode]*kB
                            *cloud_.constProps(p.typeId()).thetaV()[mode];
                        
                        vibrationalETotal_[iD][mode][cell] += EVib[mode];
                        vibEn += EVib[mode];
                    }
                                    
                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    rhoMMean_[cell] += mass;
                    linearKEMean_[cell] += massBySqMagU;
                    momentumMean_[cell] += mass*p.U();
                    rotationalEMean_[cell] += p.ERot();
                    rotationalDofMean_[cell] += rotationalDof;
                    electronicETotal_[iD][cell] += 
                        electronicEnergies[p.ELevel()];
                    nParcels_[iD][cell] += 1.0;
                    mccSpecies_[iD][cell] += massBySqMagU;
                    
                    nParcelsXnParticle_[iD][cell] += nParticles;
                    rhoNMeanXnParticle_[cell] += nParticles;
                    rhoMMeanXnParticle_[cell] += mass*nParticles;
                    momentumMeanXnParticle_[cell] += mass*(p.U())*nParticles;
                    linearKEMeanXnParticle_[cell] += massBySqMagU*nParticles;
                                                        
                    
                    muu_[cell] += mass*sqr(p.U().x());
                    muv_[cell] += mass*p.U().x()*p.U().y();
                    muw_[cell] += mass*p.U().x()*p.U().z();
                    mvv_[cell] += mass*sqr(p.U().y());
                    mvw_[cell] += mass*p.U().y()*p.U().z();
                    mww_[cell] += mass*sqr(p.U().z());
                    
                    mcc_[cell] += massBySqMagU;
                    mccu_[cell] += massBySqMagU*(p.U().x());
                    mccv_[cell] += massBySqMagU*(p.U().y());
                    mccw_[cell] += massBySqMagU*(p.U().z());
                    
                    eu_[cell] += (p.ERot() + vibEn)*p.U().x();
                    ev_[cell] += (p.ERot() + vibEn)*p.U().y();
                    ew_[cell] += (p.ERot() + vibEn)*p.U().z();
                    e_[cell] += (p.ERot() + vibEn);
                    
                    if (rotationalDof > VSMALL)
                    {
                        rhoNMeanInt_[cell] += 1.0;
                    }
                    
                    const label& nElecLevels =                 
                       cloud_.constProps(p.typeId()).numberOfElectronicLevels();
                    
                    if (nElecLevels > 1)
                    {
                        molsElec_[cell] += 1.0;
                        
                        if (p.ELevel() == 0)
                        {
                            nGroundElectronicLevel_[iD][cell]++;
                        }
                        if (p.ELevel() == 1)
                        {
                            nFirstElectronicLevel_[iD][cell]++;
                        }
                    }
                    
                    if (measureClassifications_)
                    {
                        const label& classification = p.classification();
                        
                        if (classification == 0)
                        {
                            nClassI_[cell] += 1.0;
                        }
                        
                        if (classification == 1)
                        {
                            nClassII_[cell] += 1.0;
                        }
                        
                        if (classification == 2)
                        {
                            nClassIII_[cell] += 1.0;
                        }
                    }
                }
            }
            
            //Info<< "fields myCal" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl; // TODO VINCENT
            
            // obtain collision quality measurements
            forAll(cloud_.cellPropMeasurements().collisionSeparation(), cell)
            {
                collisionSeparation_[cell] += 
                    cloud_.cellPropMeasurements().collisionSeparation()[cell];
                nColls_[cell] += cloud_.cellPropMeasurements().nColls()[cell];
            }
            
            //- obtain boundary measurements
            forAll(typeIds_, i)
            {
                const label iD = typeIds_[i];
                
                forAll(mesh_.boundaryMesh(), j)
                {                
                    forAll(mesh_.boundaryMesh()[j], k)
                    {
                        rhoNBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[iD][j][k];
                        rhoMBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoMBF()[iD][j][k];
                        linearKEBF_[j][k] += cloud_.boundaryFluxMeasurements().linearKEBF()[iD][j][k];
                        momentumBF_[j][k] += cloud_.boundaryFluxMeasurements().momentumBF()[iD][j][k];
                        rotationalEBF_[j][k] += cloud_.boundaryFluxMeasurements().rotationalEBF()[iD][j][k];
                        rotationalDofBF_[j][k] += cloud_.boundaryFluxMeasurements().rotationalDofBF()[iD][j][k];
                        qBF_[j][k] += cloud_.boundaryFluxMeasurements().qBF()[iD][j][k];
                        fDBF_[j][k] += cloud_.boundaryFluxMeasurements().fDBF()[iD][j][k];
                        speciesRhoNBF_[i][j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[iD][j][k];
                        vibrationalEBF_[i][j][k] += cloud_.boundaryFluxMeasurements().vibrationalEBF()[iD][j][k];
                        electronicEBF_[i][j][k] += cloud_.boundaryFluxMeasurements().electronicEBF()[iD][j][k];
                        mccSpeciesBF_[i][j][k] += cloud_.boundaryFluxMeasurements().mccSpeciesBF()[iD][j][k];
                        speciesRhoNIntBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNIntBF()[iD][j][k];
                        speciesRhoNElecBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNElecBF()[iD][j][k];
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
                                    .evmsBF()[iD][mode][j][k];
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
            forAll(rhoNMean_, cell)
            {
                if (rhoNMean_[cell] > SMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/nAvTimeSteps;
                    
                    rhoN_[cell] = rhoNMeanXnParticle_[cell]
                        /(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = rhoMMeanXnParticle_[cell]
                        /(nAvTimeSteps*cellVolume);
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcRhoNMean_[cell] = 0.001;
                    //dsmcRhoN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                }
                
                if (rhoNInstantaneous_[cell] > SMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    dsmcRhoN_[cell] = 0.001;
                }
            }
        }
        else
        {                  
            const label nCells = mesh_.nCells();
            
            scalarField vibT(nCells, 0.0);
            scalarField vibTForOverallT(nCells, 0.0);
            
            //- Heat capacity at constant volume/(0.5*kB), trans-rotational
            scalarField molarCv_transrot(nCells, 0.0);
            //- Heat capacity at constant pressure/(0.5*kB), trans-rotational
            scalarField molarCp_transrot(nCells, 0.0);
            scalarField molecularMass(nCells, 0.0);
            scalarField particleCv(nCells, 0.0);
            scalarField totalvDof(nCells, 0.0);
            scalarField totalvDofOverall(nCells, 0.0);
            
            forAll(rhoNMean_, cell)
            {                
                if (rhoNMean_[cell] > 1e-3) // NEW VINCENT
                {                  
                    const scalar cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/nAvTimeSteps;
                    
                    //dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                    
                    const scalar rhoNMean = rhoNMeanXnParticle_[cell]
                        /(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean = rhoMMeanXnParticle_[cell]
                        /(nAvTimeSteps*cellVolume);
                    
                    rhoN_[cell] = rhoNMean;
                    rhoM_[cell] = rhoMMean;
                    
                    UMean_[cell] = momentumMeanXnParticle_[cell]
                        /rhoMMeanXnParticle_[cell];

                    const scalar linearKEMean = 0.5
                        *linearKEMeanXnParticle_[cell] 
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
                    dsmcRhoNMean_[cell] = 0.001; 
                    //dsmcRhoN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;           
                    p_[cell] = 0.0;
                }
                
                if (rhoNInstantaneous_[cell] > SMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    dsmcRhoN_[cell] = 0.001;
                }

                //- Rotational energy mode
                const scalar totalrDof 
                (
                    rhoNMean_[cell] > SMALL
                  ? rotationalDofMean_[cell]/rhoNMean_[cell]
                  : 0.0
                );
                
                rotationalT_[cell] =
                (
                    rotationalDofMean_[cell] > SMALL
                  ? 2.0*rotationalEMean_[cell]/(kB*rotationalDofMean_[cell])
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
                            .vibrationalDegreesOfFreedom(), 0.0
                    );
                    vibTMode[i].setSize
                    (
                        cloud_.constProps(typeIds_[i])
                            .vibrationalDegreesOfFreedom(), 0.0
                    );
                }
               
                forAll(typeIds_, i)
                {
                    forAll(vibrationalETotal_[i], mode)
                    {
                        if (vibrationalETotal_[i][mode][cell] > VSMALL
                            && nParcels_[i][cell] > VSMALL
                            && degreesOfFreedomMode.size() > VSMALL)
                        {        
                            const scalar thetaV = 
                                cloud_.constProps(typeIds_[i]).thetaV()[mode];
                            
                            const scalar vibrationalEMean = 
                                vibrationalETotal_[i][mode][cell]
                               /nParcels_[i][cell];
                            
                            const scalar iMean = vibrationalEMean/(kB*thetaV);
                            
                            vibTMode[i][mode] = thetaV/log(1.0 + 1.0/iMean);

                            degreesOfFreedomMode[i][mode] = 
                                (2.0*thetaV/vibTMode[i][mode]) 
                              / (exp(thetaV/vibTMode[i][mode]) - 1.0);
                        }
                    }
                    
                    forAll(degreesOfFreedomMode[i], mode)
                    {
                        degreesOfFreedomSpecies[i] += 
                            degreesOfFreedomMode[i][mode];
                    }
                    
                    forAll(degreesOfFreedomMode[i], mode)
                    {
                        if (degreesOfFreedomSpecies[i] > SMALL)
                        {
                            vibTID[i] += vibTMode[i][mode]
                                *(degreesOfFreedomMode[i][mode]
                                /degreesOfFreedomSpecies[i]);
                        }
                    }
                    
                    totalvDof[cell] += degreesOfFreedomSpecies[i];
                    
                    if (rhoNMeanInt_[cell] > VSMALL 
                        && rhoNMean_[cell] > VSMALL 
                        && nParcels_[i][cell] > VSMALL)
                    {
                        const scalar fraction = nParcels_[i][cell]
                            /rhoNMeanInt_[cell];
                        
                        const scalar fractionOverall = nParcels_[i][cell]
                            /rhoNMean_[cell];
                        
                        totalvDofOverall[cell] += 
                            totalvDof[cell]
                            *(fractionOverall/fraction);
                        
                        vibT[cell] += vibTID[i]*fraction;
                    }
                }

                vibrationalT_[cell] = vibT[cell];
                
                //- Electronic energy mode
                scalar totalelDof = 0.0;
                scalar elecT = 0.0;
                    
                forAll(nParcels_, i) // TODO
                {
                    const scalarList& electronicEnergies = 
                        cloud_.constProps(typeIds_[i]).electronicEnergyList();
                    const labelList& degeneracies = 
                        cloud_.constProps(typeIds_[i]).degeneracyList();

                    if (   nGroundElectronicLevel_[i][cell] > SMALL 
                        && nFirstElectronicLevel_[i][cell] > SMALL 
                        && nFirstElectronicLevel_[i][cell]*degeneracies[0] 
                            != nGroundElectronicLevel_[i][cell]*degeneracies[1]
                       )
                    {
                        const scalar elecTID = 
                            (electronicEnergies[1] - electronicEnergies[0])
                           /(
                                kB*log
                                (
                                    nGroundElectronicLevel_[i][cell]*degeneracies[1]
                                   /(nFirstElectronicLevel_[i][cell]*degeneracies[0])
                                )
                            );

                        const scalar fraction = nParcels_[i][cell]
                            /molsElec_[cell];
                            
                        if (elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }
                        
                        const scalar eDof = 2.0*electronicETotal_[i][cell]
                            /nParcels_[i][cell]/(kB*elecTID);

                        totalelDof += fraction*eDof;
                    }
                    
    //                 label nElectronicLevels = cloud_.constProps(typeIds_[i]).numberOfElectronicLevels();
    //                 
    //                 if (nElectronicLevels > 1 && nParcels_[i][cell] > VSMALL && molsElec_[cell] > VSMALL)
    //                 {
    //                     const scalarList& electronicEnergies = cloud_.constProps(typeIds_[i]).electronicEnergyList();
    //                     const labelList& degeneracies = cloud_.constProps(typeIds_[i]).degeneracyList();
    //                     
    //                     scalar speciesTransT = (1.0/(3.0*kB))
    //                                             *(
    //                                                 (mccSpecies_[i][cell]/(nParcels_[i][cell]))
    //                                                 - (
    //                                                     cloud_.constProps(typeIds_[i]).mass()*mag(UMean_[cell])*mag(UMean_[cell])
    //                                                 )
    //                                             );
    //                     
    //                     scalar fraction = nParcels_[i][cell]/molsElec_[cell];
    //                     
    //                     if (speciesTransT > SMALL && electronicETotal_[i][cell] > VSMALL)
    //                     {
    //                         scalar sum1 = 0.0;
    //                         scalar sum2 = 0.0;
    //                         
    //                         forAll(electronicEnergies, ii)
    //                         {
    //                             sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(kB*speciesTransT));
    //                             sum2 += degeneracies[ii]*(electronicEnergies[ii]/(kB*speciesTransT))
    //                                         *exp(-electronicEnergies[ii]/(kB*speciesTransT));
    //                         }
    //                         
    //                         if (sum2 > VSMALL && sum1 > VSMALL)
    //                         {
    //                             scalar elecTID = (electronicETotal_[i][cell]/(kB*nParcels_[i][cell]))*(sum1/sum2);
    //                             
    //                             if (elecTID > SMALL && elecTID < GREAT)
    //                             {
    //                                 elecT += fraction*elecTID;
    //                                 
    //                                 scalar eDof = (2.0*(electronicETotal_[i][cell]/nParcels_[i][cell]))/(kB*speciesTransT);
    //                                 
    //                                 totalelDof += fraction*eDof;
    //                             }
    //                         }
    //                     }
    //                 }
                }

                electronicT_[cell] = elecT; // TODO

                overallT_[cell] = 
                    ( 
                        3.0*translationalT_[cell]
                        + totalrDof*rotationalT_[cell]
                        + totalvDof_[cell]*vibrationalT_[cell]
                        + totalelDof*electronicT_[cell]
                    ) /
                    (3.0 + totalrDof + totalvDof_[cell] + totalelDof); // TODO
                                    
                if (measureHeatFluxShearStress_)
                {                    
                    if (rhoNMean_[cell] > SMALL)
                    {
                        pressureTensor_[cell].xx() = 
                            rhoN_[cell]*
                            (
                                muu_[cell]/rhoNMean_[cell] 
                              - rhoMMean_[cell]/rhoNMean_[cell]
                                  *sqr(UMean_[cell].x())
                            );
                        pressureTensor_[cell].xy() = 
                            rhoN_[cell]*
                            (
                                muv_[cell]/rhoNMean_[cell]
                              - rhoMMean_[cell]/rhoNMean_[cell]
                                    *UMean_[cell].x()*UMean_[cell].y() 
                            );
                        pressureTensor_[cell].xz() = 
                            rhoN_[cell]*
                            (
                                muw_[cell]/rhoNMean_[cell] 
                              - rhoMMean_[cell]/rhoNMean_[cell]
                                  *UMean_[cell].x()*UMean_[cell].z()
                            );
                        
                        pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                        pressureTensor_[cell].yy() = 
                            rhoN_[cell]*
                            (
                                mvv_[cell]/(rhoNMean_[cell]) 
                              - rhoMMean_[cell]/rhoNMean_[cell]
                                    *sqr(UMean_[cell].y())
                            );
                        pressureTensor_[cell].yz() = 
                            rhoN_[cell]*
                            (
                                mvw_[cell]/rhoNMean_[cell] 
                              - rhoMMean_[cell]/rhoNMean_[cell]
                                    *UMean_[cell].y()*UMean_[cell].z()
                            );
                        
                        pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                        pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                        pressureTensor_[cell].zz() = 
                            rhoN_[cell]*
                            (
                                mww_[cell]/rhoNMean_[cell] 
                              - rhoMMean_[cell]/rhoNMean_[cell]
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
                                0.5*mccu_[cell]/rhoNMean_[cell]
                              - 0.5*mcc_[cell]/rhoNMean_[cell]*UMean_[cell].x()
                              + eu_[cell]/rhoNMean_[cell]
                              - e_[cell]/rhoNMean_[cell]*UMean_[cell].x()
                            )
                          - pressureTensor_[cell].xx()*UMean_[cell].x()
                          - pressureTensor_[cell].xy()*UMean_[cell].y()
                          - pressureTensor_[cell].xz()*UMean_[cell].z();
                                                
                        heatFluxVector_[cell].y() = 
                            rhoN_[cell]*
                            (
                                0.5*mccv_[cell]/rhoNMean_[cell]
                              - 0.5*mcc_[cell]/rhoNMean_[cell]*UMean_[cell].y()
                              + ev_[cell]/rhoNMean_[cell]
                              - e_[cell]/rhoNMean_[cell]*UMean_[cell].y()
                            )
                          - pressureTensor_[cell].yx()*UMean_[cell].x()
                          - pressureTensor_[cell].yy()*UMean_[cell].y()
                          - pressureTensor_[cell].yz()*UMean_[cell].z();
                                                
                        heatFluxVector_[cell].z() = 
                            rhoN_[cell]*
                            (
                                0.5*mccw_[cell]/rhoNMean_[cell]
                              - 0.5*mcc_[cell]/rhoNMean_[cell]*UMean_[cell].z()
                              + ew_[cell]/rhoNMean_[cell]
                              - e_[cell]/rhoNMean_[cell]*UMean_[cell].z()
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

                forAll(nParcels_, i)  
                {
                    const label iD = typeIds_[i];

                    if (rhoNMean_[cell] > SMALL)
                    {
                        const scalar molarFrac = nParcels_[i][cell]
                            /rhoNMean_[cell];
                        
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
                    molecularMass[cell] > 0.0
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
                            const scalar symmFactor = (i == qspec ? 1.0 : 2.0);
                            
                            if (
                                  nParcels_[qspec][cell] > SMALL 
                                  && translationalT_[cell] > SMALL
                               )
                            {
                                const scalar nDensQ = 
                                    nParcelsXnParticle_[qspec][cell]
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
                                
                                mcr_[i][cell] += symmFactor*sqrt(pi)*sqr(dPQ)
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
                        nColls_[cell] > SMALL
                      ? collisionSeparation_[cell]/nColls_[cell]
                      : GREAT
                    );

                    const scalar deltaT = cloud_.deltaTValue(cell);
                    
                    if (rhoN_[cell] > SMALL)
                    {
                        cr_[cell] = nColls_[cell]*cloud_.nParticles(cell)
                            /(rhoN_[cell]*mesh_.cellVolumes()[cell]
                                *nAvTimeSteps*deltaT);
                    }
                    
                    meanFreePath_[cell] = 0.0;
                    meanCollisionRate_[cell] = 0.0;

                    forAll(typeIds_, i)
                    {
                        if (rhoN_[cell] > SMALL)
                        {                    
                            const scalar nDensP = nParcelsXnParticle_[i][cell]
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
                    if (dsmcRhoN_[cell] >= 4.0) 
                    {
                        cellMfpRatio_[cell] = 1.0/mfpCellRatio_[cell];
                    }
                }
                
                if (measureClassifications_)
                {            
                    if (rhoNMean_[cell] > SMALL)
                    {
                        classIDistribution_[cell] = nClassI_[cell]
                            /rhoNMean_[cell];
                        classIIDistribution_[cell] = nClassII_[cell]
                            /rhoNMean_[cell];
                        classIIIDistribution_[cell] = nClassIII_[cell]
                            /rhoNMean_[cell];
                    }
                }
                
                if (measureErrors_)
                {
                    // TODO
                    if (
                           dsmcRhoNMean_[cell] > SMALL && Ma_[cell] > SMALL 
                        && gamma > SMALL && particleCv[cell] > SMALL
                       )
                    {
                        densityError_[cell] = 1.0
                            /sqrt(dsmcRhoNMean_[cell]*nAvTimeSteps);
                        velocityError_[cell] = (1.0/sqrt(dsmcRhoNMean_[cell]*nAvTimeSteps))*(1.0/(Ma_[cell]*sqrt(gamma)));
                        temperatureError_[cell] = 
                            (1.0/sqrt(dsmcRhoNMean_[cell]*nAvTimeSteps))
                           *sqrt(kB/particleCv[cell]);
                        pressureError_[cell] = sqrt(gamma)
                            /sqrt(dsmcRhoNMean_[cell]*nAvTimeSteps);
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
                        const scalar nParticles = cloud_.nParticles(j, k);
                        
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
                            
    //                     forAll(electronicEBF_, i)
    //                     {
    //                         label nElectronicLevels = cloud_.constProps(typeIds_[i]).numberOfElectronicLevels();
    //                         
    //                         if (nElectronicLevels > 1 && speciesRhoNBF_[i][j][k] > VSMALL && speciesRhoNElecBF_[j][k] > VSMALL)
    //                         {
    //                             const scalarList& electronicEnergies = cloud_.constProps(typeIds_[i]).electronicEnergyList();
    //                             const labelList& degeneracies = cloud_.constProps(typeIds_[i]).degeneracyList();
    //                             
    //                             scalar speciesTransT = (1.0/(3.0*kB))
    //                                                     *(
    //                                                         ((mccSpeciesBF_[i][j][k]/(speciesRhoNBF_[i][j][k])))
    //                                                         - (
    //                                                             (cloud_.constProps(typeIds_[i]).mass()
    //                                                             )*mag(UMean_.boundaryField()[j][k])*mag(UMean_.boundaryField()[j][k]))
    //                                                     );
    //                             
    //                             scalar fraction = speciesRhoNBF_[i][j][k]/speciesRhoNElecBF_[j][k];
    //                             
    //                             if (speciesTransT > VSMALL)
    //                             {
    //                                 scalar sum1 = 0.0;
    //                                 scalar sum2 = 0.0;
    //                                 
    //                                 forAll(electronicEnergies, ii)
    //                                 {
    //                                     sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(kB*speciesTransT));
    //                                     sum2 += degeneracies[ii]*(electronicEnergies[ii]/(kB*speciesTransT))
    //                                                 *exp(-electronicEnergies[ii]/(kB*speciesTransT));
    //                                 }
    //                                 
    //                                 if (sum2 > VSMALL && sum1> VSMALL)
    //                                 {
    //                                     scalar elecTID = (electronicEBF_[i][j][k]/(kB*speciesRhoNBF_[i][j][k]))*(sum1/sum2);
    //                                     
    //                                     if (elecTID > SMALL && elecTID < GREAT)
    //                                     {
    //                                         elecT += fraction*elecTID;
    //                                         
    //                                         scalar eDof = (2.0*(electronicEBF_[i][j][k]/speciesRhoNBF_[i][j][k]))/(kB*speciesTransT);
    //                                         
    //                                         totalelDof += fraction*eDof;
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
                        
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
                            /
                            (
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
                        dsmcRhoNMean_.boundaryFieldRef()[j][k] = 
                            dsmcRhoNMean_[celli];
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
                meanCollisionRate_.write();
                meanCollisionTime_.write();
                meanCollisionTimeTimeStepRatio_.write();
                meanCollisionSeparation_.write();
                cr_.write();
                SOF_.write();
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
            
            forAll(rhoNMean_, celli)
            {
                rhoNMean_[celli] = 0.0;
                rhoNInstantaneous_[celli] = 0.0;
                rhoMMean_[celli] = 0.0;
                linearKEMean_[celli] = 0.0;
                momentumMean_[celli] = vector::zero;
                rotationalEMean_[celli] = 0.0;
                rotationalDofMean_[celli] = 0.0;
                rhoNMeanInt_[celli] = 0.0;
                molsElec_[celli] = scalar(0.0),
                nClassI_[celli] = 0.0;
                nClassII_[celli] = 0.0;
                nClassIII_[celli] = 0.0;
                collisionSeparation_[celli] = 0.0;
                nColls_[celli] = 0.0;
                cr_[celli] = 0.0;
                muu_[celli] = 0.0;
                muv_[celli] = 0.0;
                muw_[celli] = 0.0;
                mvv_[celli] = 0.0;
                mvw_[celli] = 0.0;
                mww_[celli] = 0.0;
                mcc_[celli] = 0.0;
                mccu_[celli] = 0.0;
                mccv_[celli] = 0.0;
                mccw_[celli] = 0.0;
                eu_[celli] = 0.0;
                ev_[celli] = 0.0;
                ew_[celli] = 0.0;
                e_[celli] = 0.0;
                totalvDof_[celli] = 0.0;
                rhoNMeanXnParticle_[celli] = 0.0;
                rhoMMeanXnParticle_[celli] = 0.0;
                momentumMeanXnParticle_[celli] = vector::zero;
                linearKEMeanXnParticle_[celli] = 0.0;
            }

            forAll(vibrationalETotal_, i)
            {
                forAll(vibrationalETotal_[i], cell)
                {
                    vibT_[i][cell] = 0.0;
                    vDof_[i][cell] = 0.0;
                }
                
                forAll(vibrationalETotal_[i], mode)
                {
                    forAll(vibrationalETotal_[i][mode], cell)
                    {
                        vibrationalETotal_[i][mode][cell] = 0.0; 
                    }
                }
            }
            
            forAll(nParcels_, i)
            {
                forAll(nParcels_[i], cell)
                {
                    electronicETotal_[i][cell] = 0.0;
                    mccSpecies_[i][cell] = 0.0;
                    nParcels_[i][cell] = 0.0;
                    nGroundElectronicLevel_[i][cell] = 0.0;
                    nFirstElectronicLevel_[i][cell] = 0.0;
                    nParcelsXnParticle_[i][cell] = 0.0;
                    
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
    
    rhoNMean_.clear();
    rhoNInstantaneous_.clear();
    rhoMMean_.clear();
    linearKEMean_.clear();
    momentumMean_.clear();
    rotationalEMean_.clear();
    rotationalDofMean_.clear();
    rhoNMeanInt_.clear();
    molsElec_.clear();
    nClassI_.clear();
    nClassII_.clear();
    nClassIII_.clear();
    collisionSeparation_.clear();
    nColls_.clear();
    muu_.clear();
    muv_.clear();
    muw_.clear();
    mvv_.clear();
    mvw_.clear();
    mww_.clear();
    mcc_.clear();
    mccu_.clear();
    mccv_.clear();
    mccw_.clear();
    eu_.clear();
    ev_.clear();
    ew_.clear();
    e_.clear();
    totalvDof_.clear();
    rhoNMeanXnParticle_.clear();
    rhoMMeanXnParticle_.clear();
    momentumMeanXnParticle_.clear();
    linearKEMeanXnParticle_.clear();
    

    rhoNMean_.setSize(mesh_.nCells(), 0.0);
    rhoNInstantaneous_.setSize(mesh_.nCells(), 0.0);
    rhoMMean_.setSize(mesh_.nCells(), 0.0);
    linearKEMean_.setSize(mesh_.nCells(), 0.0);
    momentumMean_.setSize(mesh_.nCells(), vector::zero);
    rotationalEMean_.setSize(mesh_.nCells(), 0.0);
    rotationalDofMean_.setSize(mesh_.nCells(), 0.0);
    rhoNMeanInt_.setSize(mesh_.nCells(), 0.0);
    molsElec_.setSize(mesh_.nCells(), 0.0);
    nClassI_.setSize(mesh_.nCells(), 0.0);
    nClassII_.setSize(mesh_.nCells(), 0.0);
    nClassIII_.setSize(mesh_.nCells(), 0.0);
    collisionSeparation_.setSize(mesh_.nCells(), 0.0);
    nColls_.setSize(mesh_.nCells(), 0.0);
    cr_.setSize(mesh_.nCells(), 0.0);
    muu_.setSize(mesh_.nCells(), 0.0);
    muv_.setSize(mesh_.nCells(), 0.0);
    muw_.setSize(mesh_.nCells(), 0.0);
    mvv_.setSize(mesh_.nCells(), 0.0);
    mvw_.setSize(mesh_.nCells(), 0.0);
    mww_.setSize(mesh_.nCells(), 0.0);
    mcc_.setSize(mesh_.nCells(), 0.0);
    mccu_.setSize(mesh_.nCells(), 0.0);
    mccv_.setSize(mesh_.nCells(), 0.0);
    mccw_.setSize(mesh_.nCells(), 0.0);
    eu_.setSize(mesh_.nCells(), 0.0);
    ev_.setSize(mesh_.nCells(), 0.0);
    ew_.setSize(mesh_.nCells(), 0.0);
    e_.setSize(mesh_.nCells(), 0.0);
    totalvDof_.setSize(mesh_.nCells(), 0.0);
    rhoNMeanXnParticle_.setSize(mesh_.nCells(), 0.0);
    rhoMMeanXnParticle_.setSize(mesh_.nCells(), 0.0);
    momentumMeanXnParticle_.setSize(mesh_.nCells(), vector::zero);
    linearKEMeanXnParticle_.setSize(mesh_.nCells(), 0.0);
    
    forAll(nParcels_, i)
    {
        electronicETotal_[i].clear();
        mccSpecies_[i].clear();
        nParcels_[i].clear();
        nGroundElectronicLevel_[i].clear();
        nFirstElectronicLevel_[i].clear();
        nParcelsXnParticle_[i].clear();
        
        mfp_[i].clear();
        mcr_[i].clear();
        
        electronicETotal_[i].setSize(mesh_.nCells(), 0.0);
        mccSpecies_[i].setSize(mesh_.nCells(), 0.0);
        nParcels_[i].setSize(mesh_.nCells(), 0.0);
        nGroundElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
        nFirstElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
        nParcelsXnParticle_[i].setSize(mesh_.nCells(), 0.0);
        
        mfp_[i].setSize(mesh_.nCells(), 0.0);
        mcr_[i].setSize(mesh_.nCells(), 0.0);
    }    
        
    
    forAll(vibrationalETotal_, i)
    {
        forAll(vibrationalETotal_[i], mode)
        {
           vibrationalETotal_[i][mode].clear();
           
           vibrationalETotal_[i][mode].setSize(mesh_.nCells(), 0.0);
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

