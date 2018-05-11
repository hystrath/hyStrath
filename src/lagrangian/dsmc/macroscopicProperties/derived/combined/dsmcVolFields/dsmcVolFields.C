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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dsmcVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];
        
        if(isA<wallPolyPatch>(pPatch))
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
    n_(),
    t1_(),
    t2_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false)
{

    // standard to reading typeIds ------------ 
    const List<word>& molecules (propsDict_.lookup("typeIds"));

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
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
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
    
    nGroundElectronicLevel_.setSize(typeIds_.size());
        
    forAll(nGroundElectronicLevel_, i)
    {
        nGroundElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
    }
    
    nFirstElectronicLevel_.setSize(typeIds_.size());
        
    forAll(nFirstElectronicLevel_, i)
    {
        nFirstElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
    }

    vDof_.setSize(typeIds_.size());
    
    forAll(vDof_, i)
    {
        //vDof_[i].setSize(mesh_.nCells());
        vDof_[i].setSize(mesh_.nCells(), 0.0); // NEW VINCENT
    }
    
    vibrationalETotal_.setSize(typeIds_.size());
    
//     forAll(vibrationalETotal_, i)
//     {
//         vibrationalETotal_[i].setSize(mesh_.nCells());
//     }
//     
    electronicETotal_.setSize(typeIds_.size());
    
    forAll(electronicETotal_, i)
    {
        electronicETotal_[i].setSize(mesh_.nCells(), 0.0);
    }
    
    nParcels_.setSize(typeIds_.size());
    
    forAll(nParcels_, i)
    {
        //nParcels_[i].setSize(mesh_.nCells());
        nParcels_[i].setSize(mesh_.nCells(), 0.0); // NEW VINCENT
    }
    
    nParcelsXnParticle_.setSize(typeIds_.size());
    
    forAll(nParcelsXnParticle_, i)
    {
        //nParcelsXnParticle_[i].setSize(mesh_.nCells());
        nParcelsXnParticle_[i].setSize(mesh_.nCells(), 0.0); // NEW VINCENT
    }
    
    mccSpecies_.setSize(typeIds_.size());
    
    forAll(nParcels_, i)
    {
        //mccSpecies_[i].setSize(mesh_.nCells());
        mccSpecies_[i].setSize(mesh_.nCells(), 0.0); // NEW VINCENT
    }
    
    mfp_.setSize(typeIds_.size());
    
    forAll(mfp_, i)
    {
        //mfp_[i].setSize(mesh_.nCells());
        mfp_[i].setSize(mesh_.nCells(), 0.0); // NEW VINCENT
    }
    
    mcr_.setSize(typeIds_.size());
    
    forAll(mcr_, i)
    {
        //mcr_[i].setSize(mesh_.nCells());
        mcr_[i].setSize(mesh_.nCells(), 0.0); // NEW VINCENT
    }
    
    boundaryCells_.setSize(mesh.boundaryMesh().size());
    
//     const polyPatch& patch = mesh.boundaryMesh()[p];
    
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
    
    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());
        
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
    }
    
    if (propsDict_.found("sampleInterval"))
    {
        sampleInterval_ = 
                        readLabel(propsDict_.lookup("sampleInterval"));
    }    
    if (propsDict_.found("measureClassifications"))
    {
        measureClassifications_ = Switch(propsDict_.lookup("measureClassifications"));
    }
    
    if (propsDict_.found("measureErrors"))
    {
        measureErrors_ = Switch(propsDict_.lookup("measureErrors"));
    }
    
    if (propsDict_.found("densityOnly"))
    {
        densityOnly_ = Switch(propsDict_.lookup("densityOnly"));
    }
    
    if (propsDict_.found("measureHeatFluxShearStress"))
    {
        measureHeatFluxShearStress_ = Switch(propsDict_.lookup("measureHeatFluxShearStress"));
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
}


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

//- initial conditions
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;
    
    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize(cloud_.constProps(typeIds_[i]).vibrationalDegreesOfFreedom());
        
        forAll(vibrationalETotal_[i], j)
        {
            vibrationalETotal_[i][j].setSize(mesh_.nCells(), 0.0);
        }
    }
}


void dsmcVolFields::calculateField()
{ 
    sampleCounter_++;
    
    rhoNInstantaneous_ = 0.0;
    
    if(sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;
        
        if(densityOnly_)
        {
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label iD = findIndex(typeIds_, p.typeId());

                if(iD != -1 && p.isFree())
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
            scalar timer = mesh_.time().elapsedCpuTime();
            
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label& iD = findIndex(typeIds_, p.typeId());

                if(iD != -1 && p.isFree())
                {
                    const label cell = p.cell();
                    const scalar nParticles = cloud_.nParticles(cell);
                    const scalar mass = cloud_.constProps(p.typeId()).mass();
                    const scalar massBySqMagU = mass*(p.U() & p.U());
                    const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const scalar rotationalDof = cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();

                    scalarList EVib
                    (
                        cloud_.constProps(typeIds_[iD])
                            .vibrationalDegreesOfFreedom()
                    );
                    
                    if(EVib.size() > 0)
                    {
                        forAll(EVib, i)
                        {
                            EVib[i] = p.vibLevel()[i]
                                 *physicoChemical::k.value()                
                                 *cloud_.constProps(p.typeId()).thetaV()[i];
                            
                            vibrationalETotal_[iD][i][cell] += EVib[i];
                        }
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
                    muv_[cell] += mass*( (p.U().x()) * (p.U().y()) );
                    muw_[cell] += mass*( (p.U().x()) * (p.U().z()) );
                    mvv_[cell] += mass*sqr(p.U().y());
                    mvw_[cell] += mass*( (p.U().y()) * (p.U().z()) );
                    mww_[cell] += mass*sqr(p.U().z());
                    
                    mcc_[cell] += massBySqMagU;
                    mccu_[cell] += massBySqMagU*(p.U().x());
                    mccv_[cell] += massBySqMagU*(p.U().y());
                    mccw_[cell] += massBySqMagU*(p.U().z());
                    
                    scalar vibEn = 0.0;

                    if(EVib.size() > 0)
                    {
                        forAll(EVib, i)
                        {
                            vibEn += EVib[i];
                        }
                    }
                    
                    eu_[cell] += ( p.ERot() + vibEn )*(p.U().x());
                    ev_[cell] += ( p.ERot() + vibEn )*(p.U().y());
                    ew_[cell] += ( p.ERot() + vibEn )*(p.U().z());
                    e_[cell] += ( p.ERot() + vibEn );
                    
                    if(rotationalDof > VSMALL)
                    {
                        rhoNMeanInt_[cell] += 1.0;
                    }
                    
                    const label& nElecLevels =                 
                       cloud_.constProps(p.typeId()).numberOfElectronicLevels();
                    
                    if(nElecLevels > 1)
                    {
                        molsElec_[cell] += 1.0;
                        
                        if(p.ELevel() == 0)
                        {
                            nGroundElectronicLevel_[iD][cell]++;
                        }
                        if(p.ELevel() == 1)
                        {
                            nFirstElectronicLevel_[iD][cell]++;
                        }
                    }
                    
                    if(measureClassifications_)
                    {
                        const label& classification = p.classification();
                        
                        if(classification == 0)
                        {
                            nClassI_[cell] += 1.0;
                        }
                        
                        if(classification == 1)
                        {
                            nClassII_[cell] += 1.0;
                        }
                        
                        if(classification == 2)
                        {
                            nClassIII_[cell] += 1.0;
                        }
                    }
                }
            }
            
            Info<< "fields myCal" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl;
            
            // obtain collision quality measurements
            forAll(cloud_.cellPropMeasurements().collisionSeparation(), cell)
            {
                collisionSeparation_[cell] += 
                    cloud_.cellPropMeasurements().collisionSeparation()[cell];
                nColls_[cell] += cloud_.cellPropMeasurements().nColls()[cell];
            }
            
            // obtain boundary measurements
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
            {
                const label iD = findIndex(typeIds_, i);
            
                forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
                {                
                    forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
                    {
                        if(iD != -1)
                        { 
                            rhoNBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                            rhoMBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                            linearKEBF_[j][k] += cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                            momentumBF_[j][k] += cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                            rotationalEBF_[j][k] += cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                            rotationalDofBF_[j][k] += cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                            qBF_[j][k] += cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                            fDBF_[j][k] += cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
                            speciesRhoNBF_[iD][j][k] += cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                            vibrationalEBF_[iD][j][k] += cloud_.boundaryFluxMeasurements().vibrationalEBF()[i][j][k];
                            electronicEBF_[iD][j][k] += cloud_.boundaryFluxMeasurements().electronicEBF()[i][j][k];
                            mccSpeciesBF_[iD][j][k] += cloud_.boundaryFluxMeasurements().mccSpeciesBF()[i][j][k];
                            speciesRhoNIntBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNIntBF()[i][j][k];
                            speciesRhoNElecBF_[j][k] += cloud_.boundaryFluxMeasurements().rhoNElecBF()[i][j][k];
                        }
                    }
                }
            }
        }
        
        sampleCounter_ = 0;
    }
    
    
    if(time_.time().outputTime())
    {
        const scalar nAvTimeSteps = nTimeSteps_;
        
        if(densityOnly_)
        {
            forAll(rhoNMean_, cell)
            {
                if(rhoNMean_[cell] > VSMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
                    
                    rhoN_[cell] = (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);

                }
                else
                {
                    dsmcRhoNMean_[cell] = 0.001; // not zero so that weighted decomposition still works
                    //dsmcRhoN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                }
                
                if(rhoNInstantaneous_[cell] > VSMALL)
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
            scalarField molarconstantPressureSpecificHeat(nCells, 0.0);
            scalarField molarconstantVolumeSpecificHeat(nCells, 0.0);
            scalarField molecularMass(nCells, 0.0);
            scalarField particleConstantVolumeSpecificHeat(nCells, 0.0);
            scalarField totalvDof(nCells, 0.0);
            scalarField totalvDofOverall(nCells, 0.0);
            
            forAll(rhoNMean_, cell)
            {                
                if(rhoNMean_[cell] > VSMALL)
                {                  
                    const scalar& cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
                    
                    //dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                    
//                     Info << "dsmcRhoN_[cell] = " << dsmcRhoN_[cell] << endl;
                    
                    rhoN_[cell] = (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
                    
                    scalar rhoMMean = rhoMMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);
                    //Info << "rhoMMean is zero = crash " << tab << rhoMMean << endl;
                    UMean_[cell] = momentumMeanXnParticle_[cell] / (rhoMMean*cellVolume*nAvTimeSteps);

                    scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell] 
                                            / (cellVolume*nAvTimeSteps);
                    scalar rhoNMean = rhoNMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);

                    translationalT_[cell] = 2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                    *(linearKEMean - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell]));
                                    
                    p_[cell] = rhoN_[cell]*physicoChemical::k.value()*translationalT_[cell];
                }
                else
                {
                    dsmcRhoNMean_[cell] = 0.001; // not zero so that weighted decomposition still works
                    //dsmcRhoN_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;           
                    p_[cell] = 0.0;
                }
                
                if(rhoNInstantaneous_[cell] > VSMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    dsmcRhoN_[cell] = 0.001;
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
               
                scalarList degreesOfFreedomSpecies(typeIds_.size(),0.0);
                scalarList vibTID(vibrationalETotal_.size(), 0.0);
                
                List<scalarList> degreesOfFreedomMode;
                List<scalarList> vibTMode;
                
                degreesOfFreedomMode.setSize(typeIds_.size());
                vibTMode.setSize(typeIds_.size());
                
                forAll(degreesOfFreedomMode, iD)
                {
                    degreesOfFreedomMode[iD].setSize(cloud_.constProps(typeIds_[iD]).vibrationalDegreesOfFreedom(), 0.0);
                    vibTMode[iD].setSize(cloud_.constProps(typeIds_[iD]).vibrationalDegreesOfFreedom(), 0.0);
                }
               
                forAll(vibrationalETotal_, iD)
                {
                    forAll(vibrationalETotal_[iD], v)
                    {
                        if(vibrationalETotal_[iD][v][cell] > VSMALL
                                    && nParcels_[iD][cell] > VSMALL
                                    && degreesOfFreedomMode.size() > VSMALL)
                        {        
                            scalar thetaV = 
                                cloud_.constProps(typeIds_[iD]).thetaV()[v];
                            
                            scalar vibrationalEMean = 
                                    vibrationalETotal_[iD][v][cell]
                                    /nParcels_[iD][cell];
                            
                            scalar iMean = 
                                    vibrationalEMean
                                    /(physicoChemical::k.value()*thetaV);
                            
                            vibTMode[iD][v] = thetaV / log(1.0 + (1.0/iMean));

                            degreesOfFreedomMode[iD][v] = 
                                (2.0*thetaV/vibTMode[iD][v]) 
                                / (exp(thetaV/vibTMode[iD][v]) - 1.0);
                        }
                    }
                    
                    forAll(degreesOfFreedomMode[iD], v)
                    {
                        degreesOfFreedomSpecies[iD] += 
                                    degreesOfFreedomMode[iD][v];
                    }
                    
                    forAll(degreesOfFreedomMode[iD], v)
                    {
                        if(degreesOfFreedomSpecies[iD] > VSMALL)
                        {
                            vibTID[iD] += 
                                vibTMode[iD][v]
                                *(degreesOfFreedomMode[iD][v]
                                /degreesOfFreedomSpecies[iD]);
                        }
                    }

                    
                    totalvDof[cell] += degreesOfFreedomSpecies[iD];
                    
                    if(rhoNMeanInt_[cell] > VSMALL 
                        && rhoNMean_[cell] > VSMALL 
                        && nParcels_[iD][cell] > VSMALL)
                    {
                        scalar fraction = 
                                nParcels_[iD][cell]
                                /rhoNMeanInt_[cell];
                        
                        scalar fractionOverall = 
                                nParcels_[iD][cell]
                                /rhoNMean_[cell];
                        
                        totalvDofOverall[cell] += 
                                totalvDof[cell]
                                *(fractionOverall/fraction);
                        
                        vibT[cell] += vibTID[iD]*fraction;
                    }
                }

                vibrationalT_[cell] = vibT[cell];
                
//                 forAll(vibrationalETotal_, iD)
//                 {
//                     if(vibrationalETotal_[iD][cell] > VSMALL && nParcels_[iD][cell] > VSMALL)
//                     {        
//                         const scalar& thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
//                         
//                         scalar vibrationalEMean = (vibrationalETotal_[iD][cell]/nParcels_[iD][cell]);
//                         
//                         scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
//                         
//                         vibT_[iD][cell] = thetaV / log(1.0 + (1.0/iMean));
//                         
//                         scalar fraction = nParcels_[iD][cell]/rhoNMeanInt_[cell];
//                         
//                         vibT[cell] += vibT_[iD][cell]*fraction;
//                         
//                         vDof_[iD][cell] = fraction*(2.0*thetaV/vibT_[iD][cell]) / (exp(thetaV/vibT_[iD][cell]) - 1.0);
//                         
//                         totalvDof_[cell] += vDof_[iD][cell];
//                     }
//                 }
// 
//                 vibrationalT_[cell] = vibT[cell];
                
                // electronic temperature
                scalar totalEDof = 0.0;
                scalar elecT = 0.0;
                    
                forAll(nParcels_, iD)
                {
                    const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const labelList& degeneracies = cloud_.constProps(typeIds_[iD]).degeneracyList();

                    if(nGroundElectronicLevel_[iD][cell] > VSMALL && nFirstElectronicLevel_[iD][cell] > VSMALL && nFirstElectronicLevel_[iD][cell]*degeneracies[0] != nGroundElectronicLevel_[iD][cell]*degeneracies[1])
                    {
                        
                        scalar elecTID = (electronicEnergies[1]-electronicEnergies[0])/
                            (physicoChemical::k.value()*log((nGroundElectronicLevel_[iD][cell]*degeneracies[1])/(nFirstElectronicLevel_[iD][cell]*degeneracies[0])));

                    
                        scalar fraction = nParcels_[iD][cell]/molsElec_[cell];
                            
                        if(elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }
                        
                        
                        scalar eDof = (2.0*(electronicETotal_[iD][cell]/nParcels_[iD][cell]))/(physicoChemical::k.value()*elecTID);

                        
                        totalEDof += fraction*eDof;
                    }
                    
    //                 label nElectronicLevels = cloud_.constProps(typeIds_[iD]).numberOfElectronicLevels();
    //                 
    //                 if(nElectronicLevels > 1 && nParcels_[iD][cell] > VSMALL && molsElec_[cell] > VSMALL)
    //                 {
    //                     const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();
    //                     const labelList& degeneracies = cloud_.constProps(typeIds_[iD]).degeneracyList();
    //                     
    //                     scalar speciesTransT = (1.0/(3.0*physicoChemical::k.value()))
    //                                             *(
    //                                                 (mccSpecies_[iD][cell]/(nParcels_[iD][cell]))
    //                                                 - (
    //                                                     cloud_.constProps(typeIds_[iD]).mass()*mag(UMean_[cell])*mag(UMean_[cell])
    //                                                 )
    //                                             );
    //                     
    //                     scalar fraction = nParcels_[iD][cell]/molsElec_[cell];
    //                     
    //                     if(speciesTransT > SMALL && electronicETotal_[iD][cell] > VSMALL)
    //                     {
    //                         scalar sum1 = 0.0;
    //                         scalar sum2 = 0.0;
    //                         
    //                         forAll(electronicEnergies, ii)
    //                         {
    //                             sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
    //                             sum2 += degeneracies[ii]*(electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT))
    //                                         *exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
    //                         }
    //                         
    //                         if(sum2 > VSMALL && sum1 > VSMALL)
    //                         {
    //                             scalar elecTID = (electronicETotal_[iD][cell]/(physicoChemical::k.value()*nParcels_[iD][cell]))*(sum1/sum2);
    //                             
    //                             if(elecTID > SMALL && elecTID < GREAT)
    //                             {
    //                                 elecT += fraction*elecTID;
    //                                 
    //                                 scalar eDof = (2.0*(electronicETotal_[iD][cell]/nParcels_[iD][cell]))/(physicoChemical::k.value()*speciesTransT);
    //                                 
    //                                 totalEDof += fraction*eDof;
    //                             }
    //                         }
    //                     }
    //                 }
                }

                electronicT_[cell] = elecT;

                scalar nRotDof = 0.0;
                
                if(rhoNMean_[cell] > VSMALL)
                {
                    nRotDof = rotationalDofMean_[cell] / rhoNMean_[cell];
                }
                
                overallT_[cell] = ( 
                                        (3.0*translationalT_[cell]) 
                                        + (nRotDof*rotationalT_[cell]) 
                                        + (totalvDof_[cell]*vibrationalT_[cell])
                                        + (totalEDof*electronicT_[cell])
                                    ) /
                                    (3.0 + nRotDof + totalvDof_[cell] + totalEDof);
                                    
                if(measureHeatFluxShearStress_)
                {                    
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
                }

                totalvDof_[cell] = 0.0;
                 
                forAll(nParcels_, iD)  
                {
                    label typeId = typeIds_[iD];

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

                if(gamma > 0.1 && gasConstant > 1.0 && translationalT_[cell] > 1.0)
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
                        forAll(typeIds_, qspec)
                        {
                            const scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());
                            const scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());
                            const scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();
                            const scalar symmFactor = (iD == qspec ? 1.0 : 2.0);
                            
                            if(nParcels_[qspec][cell] > SMALL && translationalT_[cell] > SMALL)
                            {
                                const scalar nDensQ = nParcelsXnParticle_[qspec][cell]/(mesh_.cellVolumes()[cell]*nTimeSteps_);
                                const scalar reducedMass = (cloud_.constProps(typeIds_[iD]).mass()*cloud_.constProps(typeIds_[qspec]).mass())
                                    / (cloud_.constProps(typeIds_[iD]).mass()+cloud_.constProps(typeIds_[qspec]).mass());
                                
                                mfp_[iD][cell] += pi*dPQ*dPQ*nDensQ*pow(mfpReferenceTemperature_/translationalT_[cell],omegaPQ-0.5)*sqrt(1.0+massRatio); //Bird, eq (4.76)
                                
                                mcr_[iD][cell] += symmFactor*sqrt(pi)*dPQ*dPQ*nDensQ*pow(translationalT_[cell]/mfpReferenceTemperature_,1.0-omegaPQ)
                                                    *sqrt(2.0*physicoChemical::k.value()*mfpReferenceTemperature_/reducedMass); // Bird, eq (4.74)
                            }
                        }
                        
                        if(mfp_[iD][cell] > VSMALL)
                        {
                            mfp_[iD][cell] = 1.0/mfp_[iD][cell];
                        }
                    }
                    
                    meanFreePath_[cell] = 0.0;
                    mfpCellRatio_[cell] = 0.0;
                    cellMfpRatio_[cell] = 0.0;
                    meanCollisionRate_[cell] = 0.0;
                    meanCollisionTime_[cell] = 0.0;
                    meanCollisionTimeTimeStepRatio_[cell] = 0.0;
                    meanCollisionSeparation_[cell] = 0.0;
                    
                    if(nColls_[cell] > VSMALL)
                    {
                        meanCollisionSeparation_[cell] = collisionSeparation_[cell]/nColls_[cell];
                    }
                    else
                    {
                       meanCollisionSeparation_[cell] = GREAT; 
                    }
                    
                    const scalar deltaT = cloud_.deltaTValue(cell);
                    
                    cr_[cell] = nColls_[cell]*cloud_.nParticles(cell)
                        /(rhoN_[cell]*mesh_.cellVolumes()[cell]*nTimeSteps_*deltaT); // NEW VINCENT 16/04/2018
                    
                    forAll(mfp_, iD)
                    {
                        if(rhoN_[cell] > VSMALL)
                        {                    
                            scalar nDensP = (nParcelsXnParticle_[iD][cell])/(mesh_.cellVolumes()[cell]*nTimeSteps_);
                            
                            meanFreePath_[cell] += mfp_[iD][cell]*nDensP/rhoN_[cell]; //Bird, eq (4.77)
                            
                            meanCollisionRate_[cell] += mcr_[iD][cell]*nDensP/rhoN_[cell]; //Bird, eq (1.38)
                        }
                    }

                    if(meanFreePath_[cell] < VSMALL)
                    {
                        meanFreePath_[cell] = GREAT;
                    }
    
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
                        mfp_[iD][cell] = 0.0;
                        mcr_[iD][cell] = 0.0;
                        
                    }
                    
                    if(meanFreePath_[cell] != GREAT)
                    {
                        scalar largestCellDimension = 0.0;

                        const labelList& pLabels(mesh_.cells()[cell].labels(mesh_.faces()));
                        pointField pLocal(pLabels.size(), vector::zero);

                        forAll (pLabels, pointi)
                        {
                            pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                        }
                        
                        scalarField dimension;
                        
                        dimension.setSize(2, 0.0);

                        dimension[0] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
                        dimension[1] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
//                         dimension[2] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
                        
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
                        
                        if(meanFreePath_[cell] > VSMALL && meanCollisionSeparation_[cell] > VSMALL)
                        {
                            SOF_[cell] = meanCollisionSeparation_[cell]/meanFreePath_[cell];
                        }
                    }
                    else
                    {
                        mfpCellRatio_[cell] = GREAT;
                        SOF_[cell] = GREAT;
                    }
                    
                    // when no p in cell, refines when it should not
                    // the condition should eliminates this
                    if (dsmcRhoN_[cell] >= 4.0) 
                    {
                        cellMfpRatio_[cell] = 1.0/mfpCellRatio_[cell];
                    }
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
                    if(dsmcRhoNMean_[cell] > VSMALL && Ma_[cell] > VSMALL && gamma > VSMALL && particleConstantVolumeSpecificHeat[cell] > VSMALL)
                    {
                        densityError_[cell] = 1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_);
                        velocityError_[cell] = (1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_))*(1.0/(Ma_[cell]*sqrt(gamma)));
                        temperatureError_[cell] = (1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_))
                            *sqrt(physicoChemical::k.value()/particleConstantVolumeSpecificHeat[cell]);
                        pressureError_[cell] = sqrt(gamma)/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_);
                    }
                    
                }
            }
            
//             dsmcRhoN_.boundaryField() = dsmcRhoN_.boundaryField().boundaryInternalField();
            
            
//             rhoN_.correctBoundaryConditions();
//             rhoM_.correctBoundaryConditions();
            
            List<scalarField> vibTBF(mesh_.boundaryMesh().size());
            List<scalarField> molecularMassBoundary(mesh_.boundaryMesh().size());
            List<scalarField> molarconstantPressureSpecificHeatBoundary(mesh_.boundaryMesh().size());
            List<scalarField> molarconstantVolumeSpecificHeatBoundary(mesh_.boundaryMesh().size());
            List<scalarField> particleConstantVolumeSpecificHeatBoundary(mesh_.boundaryMesh().size());

            // computing boundary measurements
            forAll(rhoNBF_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                //const vectorField& fC = patch.faceCentres();

                vibTBF[j].setSize(patch.size(), 0.0);
                molecularMassBoundary[j].setSize(patch.size(), 0.0);
                molarconstantPressureSpecificHeatBoundary[j].setSize(patch.size(), 0.0);
                molarconstantVolumeSpecificHeatBoundary[j].setSize(patch.size(), 0.0);
                particleConstantVolumeSpecificHeatBoundary[j].setSize(patch.size(), 0.0);
                
                if(isA<wallPolyPatch>(patch))
                {                               
                    forAll(rhoN_.boundaryField()[j], k)
                    {                        
                        const scalar nParticles = cloud_.nParticles(j, k);
                        
                        rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*nParticles/nAvTimeSteps;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*nParticles/nAvTimeSteps;
                        
                        if(rhoM_.boundaryFieldRef()[j][k] > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]*nParticles/(rhoM_.boundaryField()[j][k]*nAvTimeSteps);
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                        }
                            
                        scalar rhoMMean = rhoMBF_[j][k]*nParticles/nAvTimeSteps;
                        scalar linearKEMean = linearKEBF_[j][k]*nParticles/nAvTimeSteps;
                        scalar rhoNMean = rhoNBF_[j][k]*nParticles/nAvTimeSteps;
                        
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
                        
//                         forAll(vibrationalEBF_, i)
//                         {
//                             if(rhoNBF_[j][k] > VSMALL)
//                             {                       
//                                 molecularMassBoundary[j][k] +=  cloud_.constProps(typeIds_[i]).mass()
//                                             *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
//                                             
//                                 molarconstantPressureSpecificHeatBoundary[j][k] += (5.0 + cloud_.constProps(typeIds_[i]).rotationalDegreesOfFreedom())
//                                             *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
//                                             
//                                 molarconstantVolumeSpecificHeatBoundary[j][k] += (3.0 + cloud_.constProps(typeIds_[i]).rotationalDegreesOfFreedom())
//                                             *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
//                             }
//                             
//                             if(vibrationalEBF_[i][j][k] > VSMALL && speciesRhoNBF_[i][j][k] > VSMALL && speciesRhoNIntBF_[j][k] > VSMALL)
//                             {        
//                                 const scalar& thetaV = cloud_.constProps(typeIds_[i]).thetaV();
//                                 
//                                 scalar vibrationalEMean = (vibrationalEBF_[i][j][k]/speciesRhoNBF_[i][j][k]);
//                                 
//                                 scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
//                                 
//                                 vibTBF_[i][j][k] = thetaV / log(1.0 + (1.0/iMean));
//                                 
//                                 scalar fraction = speciesRhoNBF_[i][j][k]/speciesRhoNIntBF_[j][k];
//                                 
//                                 vDofBF_[i][j][k] = fraction*(2.0*thetaV/vibTBF_[i][j][k]) / (exp(thetaV/vibTBF_[i][j][k]) - 1.0);
//                                 
//                                 vibTBF[j][k] += fraction*vibTBF_[i][j][k];
//                                 
//                                 totalvDofBF_[j][k] += vDofBF_[i][j][k];
//                             }
// 
//                         }
//                         
//                         if(totalvDofBF_[j][k] > VSMALL)
//                         {
//                             vibrationalT_.boundaryFieldRef()[j][k] = vibTBF[j][k];
//                         }
//                         else
//                         {
//                             vibrationalT_.boundaryFieldRef()[j][k] = 0.0;
//                         }
                        
                        // electronic temperature
                        scalar totalEDof = 0.0;
                        scalar elecT = 0.0;
                            
    //                     forAll(electronicEBF_, i)
    //                     {
    //                         label nElectronicLevels = cloud_.constProps(typeIds_[i]).numberOfElectronicLevels();
    //                         
    //                         if(nElectronicLevels > 1 && speciesRhoNBF_[i][j][k] > VSMALL && speciesRhoNElecBF_[j][k] > VSMALL)
    //                         {
    //                             const scalarList& electronicEnergies = cloud_.constProps(typeIds_[i]).electronicEnergyList();
    //                             const labelList& degeneracies = cloud_.constProps(typeIds_[i]).degeneracyList();
    //                             
    //                             scalar speciesTransT = (1.0/(3.0*physicoChemical::k.value()))
    //                                                     *(
    //                                                         ((mccSpeciesBF_[i][j][k]/(speciesRhoNBF_[i][j][k])))
    //                                                         - (
    //                                                             (cloud_.constProps(typeIds_[i]).mass()
    //                                                             )*mag(UMean_.boundaryField()[j][k])*mag(UMean_.boundaryField()[j][k]))
    //                                                     );
    //                             
    //                             scalar fraction = speciesRhoNBF_[i][j][k]/speciesRhoNElecBF_[j][k];
    //                             
    //                             if(speciesTransT > VSMALL)
    //                             {
    //                                 scalar sum1 = 0.0;
    //                                 scalar sum2 = 0.0;
    //                                 
    //                                 forAll(electronicEnergies, ii)
    //                                 {
    //                                     sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
    //                                     sum2 += degeneracies[ii]*(electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT))
    //                                                 *exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
    //                                 }
    //                                 
    //                                 if(sum2 > VSMALL && sum1> VSMALL)
    //                                 {
    //                                     scalar elecTID = (electronicEBF_[i][j][k]/(physicoChemical::k.value()*speciesRhoNBF_[i][j][k]))*(sum1/sum2);
    //                                     
    //                                     if(elecTID > SMALL && elecTID < GREAT)
    //                                     {
    //                                         elecT += fraction*elecTID;
    //                                         
    //                                         scalar eDof = (2.0*(electronicEBF_[i][j][k]/speciesRhoNBF_[i][j][k]))/(physicoChemical::k.value()*speciesTransT);
    //                                         
    //                                         totalEDof += fraction*eDof;
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
                        
                        electronicT_.boundaryFieldRef()[j][k] = elecT;
                        
                        scalar nRotDof = 0.0;
                        
                        if(rhoNBF_[j][k] > VSMALL)
                        {
                            nRotDof = rotationalDofBF_[j][k] / rhoNBF_[j][k];
                        }
                        
                        overallT_.boundaryFieldRef()[j][k] = ( 
                                                (3.0*translationalT_.boundaryField()[j][k]) 
                                                + (nRotDof*rotationalT_.boundaryField()[j][k]) 
                                                + (totalvDofBF_[j][k]*vibrationalT_.boundaryField()[j][k])
                                                + (totalEDof*elecT)
                                            ) /
                                            (3.0 + nRotDof + totalvDofBF_[j][k] + totalEDof);
                                            
                        totalvDofBF_[j][k] = 0.0;
                        
                        /**************************************************************************************************************/
                        
                        particleConstantVolumeSpecificHeatBoundary[j][k] = molarconstantVolumeSpecificHeatBoundary[j][k]/6.02214129e23;
                        
                        scalar gasConstant = 0.0;
                        scalar gamma = 0.0;
                        scalar speedOfSound = 0.0;
                        
                        if(molecularMassBoundary[j][k] > VSMALL)
                        {
                            gasConstant = physicoChemical::k.value()/molecularMassBoundary[j][k]; // R = k/m
                        }
                        
                        if(molarconstantVolumeSpecificHeatBoundary[j][k] > VSMALL)
                        {
                            gamma = molarconstantPressureSpecificHeatBoundary[j][k]/molarconstantVolumeSpecificHeatBoundary[j][k]; // gamma = cP/cV
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
                        
//                         rhoN_.boundaryFieldRef()[j][k] = rhoN_[boundaryCells_[j][k]];
//                         rhoM_.boundaryFieldRef()[j][k] = rhoM_[boundaryCells_[j][k]];
                    }
                    
                    /*p_.boundaryFieldRef()[j] = fD_.boundaryField()[j]
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
                    const scalarField magt2 = mag(t2);
                    forAll(t2, f)
                    {
                        if(magt2[f] > 0.0)
                        {
                            t2[f] /= magt2[f];
                        }
                    }
                    
                    tau_.boundaryFieldRef()[j] = sqrt(
                        sqr(fD_.boundaryField()[j] & t1)
                        + sqr(fD_.boundaryField()[j] & t2));
                    */ // DELETED VINCENT 10/04/2018
                    
                    forAll(p_.boundaryField()[j], facei)
                    {
                        p_.boundaryFieldRef()[j][facei] = 
                            fD_.boundaryField()[j][facei] & n_[j][facei];
                        
                        tau_.boundaryFieldRef()[j][facei] = sqrt
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
                
                if(isA<polyPatch>(patch))
                {
                    if(!isA<emptyPolyPatch>(patch))
                    {
                        if(!isA<cyclicPolyPatch>(patch)) 
                        {
                            if(!isA<wallPolyPatch>(patch))
                            {
                                forAll(boundaryCells_[j], k)
                                {       
                                    translationalT_.boundaryFieldRef()[j][k] = translationalT_[boundaryCells_[j][k]];
                                    rotationalT_.boundaryFieldRef()[j][k] = rotationalT_[boundaryCells_[j][k]];
                                    vibrationalT_.boundaryFieldRef()[j][k] = vibrationalT_[boundaryCells_[j][k]];
                                    overallT_.boundaryFieldRef()[j][k] = overallT_[boundaryCells_[j][k]];
                                    dsmcRhoNMean_.boundaryFieldRef()[j][k] = dsmcRhoNMean_[boundaryCells_[j][k]];
                                    rhoN_.boundaryFieldRef()[j][k] = rhoN_[boundaryCells_[j][k]];
                                    rhoM_.boundaryFieldRef()[j][k] = rhoM_[boundaryCells_[j][k]];
                                    p_.boundaryFieldRef()[j][k] = p_[boundaryCells_[j][k]];
                                    Ma_.boundaryFieldRef()[j][k] = Ma_[boundaryCells_[j][k]];
                                    UMean_.boundaryFieldRef()[j][k] = UMean_[boundaryCells_[j][k]];
                                }
                            }
                            
                            if(measureMeanFreePath_)
                            {
                                forAll(boundaryCells_[j], k)
                                {
                                    meanFreePath_.boundaryFieldRef()[j][k] = meanFreePath_[boundaryCells_[j][k]];
                                    SOF_.boundaryFieldRef()[j][k] = SOF_[boundaryCells_[j][k]];
                                    mfpCellRatio_.boundaryFieldRef()[j][k] = mfpCellRatio_[boundaryCells_[j][k]];
                                    meanCollisionRate_.boundaryFieldRef()[j][k] = meanCollisionRate_[boundaryCells_[j][k]];
                                    meanCollisionTime_.boundaryFieldRef()[j][k] = meanCollisionTime_[boundaryCells_[j][k]];
                                    meanCollisionTimeTimeStepRatio_.boundaryFieldRef()[j][k] = meanCollisionTimeTimeStepRatio_[boundaryCells_[j][k]];
                                }
                            }
                            
                            if(measureHeatFluxShearStress_)
                            {
                                forAll(boundaryCells_[j], k)
                                {
                                    shearStressTensor_.boundaryFieldRef()[j][k] = shearStressTensor_[boundaryCells_[j][k]];
                                    heatFluxVector_.boundaryFieldRef()[j][k] = heatFluxVector_[boundaryCells_[j][k]];
                                    pressureTensor_.boundaryFieldRef()[j][k] = pressureTensor_[boundaryCells_[j][k]];
                                }
                            }
                            
                            if(measureClassifications_)
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
                
                /*if(measureMeanFreePath_)
                {
                    if(!isA<emptyPolyPatch>(patch))
                    {
                        if(!(isA<cyclicPolyPatch>(patch) or isA<wallPolyPatch>(patch)))
                        {
                            forAll(boundaryCells_[j], k)
                            {
                                meanFreePath_.boundaryFieldRef()[j][k] = meanFreePath_[boundaryCells_[j][k]];
                                SOF_.boundaryFieldRef()[j][k] = SOF_[boundaryCells_[j][k]];
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
                            if(!(isA<cyclicPolyPatch>(patch) or isA<wallPolyPatch>(patch)))
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
                
                if(measureHeatFluxShearStress_)
                { 
                    if(isA<polyPatch>(patch))
                    {
                        if(!isA<emptyPolyPatch>(patch))
                        {
                            if(!(isA<cyclicPolyPatch>(patch) or isA<wallPolyPatch>(patch)))
                            {
                                forAll(boundaryCells_[j], k)
                                {
                                    shearStressTensor_.boundaryFieldRef()[j][k] = shearStressTensor_[boundaryCells_[j][k]];
                                    heatFluxVector_.boundaryFieldRef()[j][k] = heatFluxVector_[boundaryCells_[j][k]];
                                    pressureTensor_.boundaryFieldRef()[j][k] = pressureTensor_[boundaryCells_[j][k]];
                                }
                            }
                        }
                    }
                }*/ // DELETED VINCENT 10/04/2018
            } // end loop boundary cells for generalPatches
            
            if(measureMeanFreePath_)
            {
                meanFreePath_.write();
                mfpCellRatio_.write();
                meanCollisionRate_.write();
                meanCollisionTime_.write();
                meanCollisionTimeTimeStepRatio_.write();
                meanCollisionSeparation_.write();
                SOF_.write();
                cr_.write();
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
            
            if(measureHeatFluxShearStress_)
            {
                heatFluxVector_.write();
                pressureTensor_.write();
                shearStressTensor_.write();
            }
            
            p_.write();
            translationalT_.write();
            rotationalT_.write();
            vibrationalT_.write();
            electronicT_.write();
            overallT_.write();
            q_.write();
            tau_.write();
            Ma_.write();
            UMean_.write();
            fD_.write();
        }
        
        //- reset
        if(time_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0.0;
            
            forAll(rhoNMean_, c)
            {
                rhoNMean_[c] = 0.0;
                rhoNInstantaneous_[c] = 0.0; // NEW VINCENT 03/06/2017
                rhoMMean_[c] = 0.0;
                linearKEMean_[c] = 0.0;
                momentumMean_[c] = vector::zero;
                rotationalEMean_[c] = 0.0;
                rotationalDofMean_[c] = 0.0;
                rhoNMeanInt_[c] = 0.0;
                molsElec_[c] = scalar(0.0),
                nClassI_[c] = 0.0;
                nClassII_[c] = 0.0;
                nClassIII_[c] = 0.0;
                collisionSeparation_[c] = 0.0;
                nColls_[c] = 0.0;
                muu_[c] = 0.0;
                muv_[c] = 0.0;
                muw_[c] = 0.0;
                mvv_[c] = 0.0;
                mvw_[c] = 0.0;
                mww_[c] = 0.0;
                mcc_[c] = 0.0;
                mccu_[c] = 0.0;
                mccv_[c] = 0.0;
                mccw_[c] = 0.0;
                eu_[c] = 0.0;
                ev_[c] = 0.0;
                ew_[c] = 0.0;
                e_[c] = 0.0;
                totalvDof_[c] = 0.0; // NEW VINCENT 03/06/2017
                rhoNMeanXnParticle_[c] = 0.0;
                rhoMMeanXnParticle_[c] = 0.0;
                momentumMeanXnParticle_[c] = vector::zero;
                linearKEMeanXnParticle_[c] = 0.0;
                
                cr_[c] = 0.0; // NEW VINCENT 16/04/2018
            }

            // NEW VINCENT
            forAll(vibrationalETotal_, iDv)
            {
                forAll(vibrationalETotal_[iDv], cell)
                {
                    vibT_[iDv][cell] = 0.0; // NEW VINCENT 03/06/2017
                    vDof_[iDv][cell] = 0.0; // NEW VINCENT 03/06/2017
                }
                
                forAll(vibrationalETotal_[iDv], v)
                {
                    forAll(vibrationalETotal_[iDv][v], cell)
                    {
                        vibrationalETotal_[iDv][v][cell] = 0.0; 
                    }
                }
            }
            
            forAll(nParcels_, iD)
            {
                forAll(nParcels_[iD], cell)
                {
                    electronicETotal_[iD][cell] = 0.0;
                    mccSpecies_[iD][cell] = 0.0;
                    nParcels_[iD][cell] = 0.0;
                    nGroundElectronicLevel_[iD][cell] = 0.0;
                    nFirstElectronicLevel_[iD][cell] = 0.0;
                    nParcelsXnParticle_[iD][cell] = 0.0;
                    
                    mfp_[iD][cell] = 0.0; // NEW VINCENT 03/06/2017
                    mcr_[iD][cell] = 0.0; // NEW VINCENT 03/06/2017
                }
            }
            // END NEW VINCENT
            
            
            // reset boundary information
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
                vibTxvDofBF_[j] = 0.0; // NEW VINCENT 03/06/2017
                totalvDofBF_[j] = 0.0; // NEW VINCENT 03/06/2017
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
            
            forAll(vibTBF_, i) // NEW VINCENT 03/06/2017
            {       
                forAll(vibTBF_[i], j)
                {
                    vibTBF_[i][j] = 0.0;
                    vDofBF_[i][j] = 0.0;
                }    
            }  
        }
        
        if(averagingAcrossManyRuns_)
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
    
    cr_.setSize(mesh_.nCells(), 0.0); // NEW VINCENT 16/04/2018
    
    forAll(nParcels_, iD)
    {
        electronicETotal_[iD].clear();
        mccSpecies_[iD].clear();
        nParcels_[iD].clear();
        nGroundElectronicLevel_[iD].clear();
        nFirstElectronicLevel_[iD].clear();
        nParcelsXnParticle_[iD].clear();
        
        mfp_[iD].clear();
        mcr_[iD].clear();
        
        electronicETotal_[iD].setSize(mesh_.nCells(), 0.0);
        mccSpecies_[iD].setSize(mesh_.nCells(), 0.0);
        nParcels_[iD].setSize(mesh_.nCells(), 0.0);
        nGroundElectronicLevel_[iD].setSize(mesh_.nCells(), 0.0);
        nFirstElectronicLevel_[iD].setSize(mesh_.nCells(), 0.0);
        nParcelsXnParticle_[iD].setSize(mesh_.nCells(), 0.0);
        
        mfp_[iD].setSize(mesh_.nCells(), 0.0);
        mcr_[iD].setSize(mesh_.nCells(), 0.0);
    }    
        
    
    forAll(vibrationalETotal_, iD)
    {
        forAll(vibrationalETotal_[iD], v)
        {
           vibrationalETotal_[iD][v].clear();
           
           vibrationalETotal_[iD][v].setSize(mesh_.nCells(), 0.0);
        }
    }
    
    forAll(vDof_, iDv)
    {
        vibT_[iDv].clear();
        vDof_[iDv].clear();
        
        vibT_[iDv].setSize(mesh_.nCells(), 0.0);
        vDof_[iDv].setSize(mesh_.nCells(), 0.0);
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
    }
    
    
    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[p];
        
        boundaryCells_[p].clear();
        boundaryCells_[p].setSize(patch.size());
        
        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }
    
    //Info << "Reset volFields!" << endl;
}


//- write field
void dsmcVolFields::writeField()
{}

void dsmcVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //

