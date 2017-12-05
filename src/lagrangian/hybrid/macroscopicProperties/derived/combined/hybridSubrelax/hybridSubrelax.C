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
    zoneToModify_(propsDict_.lookup("dsmcZoneName")),
    itBeforeRelaxing_(propsDict_.lookupOrDefault<label>(
        "itBeforeRelaxing", 1000)),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
//    translationalTs_(),
//    rotationalTs_(),
    vibrationalTs_(),
    rhoNType_(),
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
    translationalT0_
    (
        IOobject
        (
            "Tt0_",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    rotationalT0_
    (
        IOobject
        (
            "Tr0_",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    vibrationalTs0_(),
    rhoNType0_(),
    UMean0_
    (
        IOobject
        (
            "U0_",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    momentumMean_(mesh_.nCells(), vector::zero),
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
    vibrationalETotal_(),
    nParcels_(),

    nTimeSteps_(0),
    typeIds_(),

    startingTime_(propsDict_.lookupOrDefault<fileName>("fileName", "0")),

/*    pName_(propsDict_.lookup("pName")),
    UName_(propsDict_.lookup("UName")),
    TtName_(propsDict_.lookup("TtName")),
    TrName_(propsDict_.lookup("TrName")),
    rhoNNames_(propsDict_.lookup("rhoNames")),
    TvNames_(propsDict_.lookup("TvNames")),
    mass_(readList<scalar>(propsDict_.lookup("masses"))),*/

//    meshChanging_(true),
    subResetting_(true),
    subrelaxation_(true)
{

    label zoneID = mesh_.cellZones().findZoneID(zoneToModify_);

    if(zoneID == -1)
    {
        FatalErrorIn("zoneToModify")
            << "Cannot find region: " << zoneToModify_ << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

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
    rhoNType_.setSize(typeIds_.size());
    rhoNType0_.setSize(typeIds_.size());
    vibrationalTs_.setSize(typeIds_.size());
    vibrationalTs0_.setSize(typeIds_.size());

    vibrationalETotal_.setSize(typeIds_.size());
    nParcels_.setSize(typeIds_.size());

    forAll(molecules, moleculeI)
    {
        rhoNType_.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "rhoN_" + molecules[moleculeI],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0.0", dimless/dimVolume, 0.0)
            )
        );
        rhoNType0_.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "rhoN0_" + molecules[moleculeI],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0.0", dimless/dimVolume, 0.0)
            )
        );
        vibrationalTs_.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "vibrationalT_" + molecules[moleculeI],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0.0", dimTemperature, 0.0)
            )
        );
        vibrationalTs0_.set
        (
            moleculeI,
            new volScalarField
            (
                IOobject
                (
                    "vibrationalT0_" + molecules[moleculeI],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0.0", dimTemperature, 0.0)
            )
        );

        vibrationalETotal_[moleculeI].setSize(mesh_.nCells(), 0.0);
        nParcels_[moleculeI].setSize(mesh_.nCells(), 0.0);
    }

    setInitialFields(molecules, startingTime_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hybridSubrelax::~hybridSubrelax()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//- initial conditions
void hybridSubrelax::createField()
{
    Info << "Initialising hybridSubrelax field" << endl;
}


void hybridSubrelax::calculateField()
{
    if(subrelaxation_)
    {
        if(subResetting_)
        {
            label localSteps = nTimeSteps_ % itBeforeRelaxing_;
            nTimeSteps_ -= localSteps;
        }

        nTimeSteps_++;

        scalar relaxFactor = 1.0 / scalar(itBeforeRelaxing_);
        scalar factor1 = (1.0 - relaxFactor);

        rhoNMean_ *= 0.0;
        rhoMMean_ *= 0.0;
        rhoNMeanInt_ *= 0.0;
        linearKEMean_ *= 0.0;
        momentumMean_ *= 0.0;
        rotationalEMean_ *= 0.0;
        rotationalDofMean_ *= 0.0;
        muu_ *= 0.0;
        muv_ *= 0.0;
        muw_ *= 0.0;
        mvv_ *= 0.0;
        mvw_ *= 0.0;
        mww_ *= 0.0;
        mcc_ *= 0.0;
        mccu_ *= 0.0;
        mccv_ *= 0.0;
        mccw_ *= 0.0;
        eu_ *= 0.0;
        ev_ *= 0.0;
        ew_ *= 0.0;
        e_ *= 0.0;
        forAll(typeIds_, iD)
        {
            vibrationalETotal_[iD] *= 0.0;
            nParcels_[iD] *= 0.0;
        }

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

                rhoNMean_[cell] += 1.0;
                rhoMMean_[cell] += mass;
                linearKEMean_[cell] += mass * (p.U() & p.U());
                momentumMean_[cell] += mass * p.U();
                rotationalEMean_[cell] += p.ERot();
                rotationalDofMean_[cell] += rotationalDof; 
                vibrationalETotal_[iD][cell] += p.EVib();
                nParcels_[iD][cell] += 1.0;
            
                muu_[cell] += mass * sqr(p.U().x());
                muv_[cell] += mass * ( (p.U().x()) * (p.U().y()) );
                muw_[cell] += mass * ( (p.U().x()) * (p.U().z()) );
                mvv_[cell] += mass * sqr(p.U().y());
                mvw_[cell] += mass * ( (p.U().y()) * (p.U().z()) );
                mww_[cell] += mass * sqr(p.U().z());
            
                mcc_[cell] += mass * mag(p.U()) * mag(p.U());
                mccu_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().x());
                mccv_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().y());
                mccw_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().z());
            
                eu_[cell] += ( p.ERot() + p.EVib() ) * (p.U().x());
                ev_[cell] += ( p.ERot() + p.EVib() ) * (p.U().y());
                ew_[cell] += ( p.ERot() + p.EVib() ) * (p.U().z());
                e_[cell] += ( p.ERot() + p.EVib() );

                if(rotationalDof > VSMALL)
                {
                    rhoNMeanInt_[cell] += 1.0;
                }
            }
        }

        List<scalarField> rhoNBF_(mesh_.boundaryMesh().size());
        List<scalarField> rhoMBF_(mesh_.boundaryMesh().size());
        List<scalarField> linearKEBF_(mesh_.boundaryMesh().size());
        List<vectorField> momentumBF_(mesh_.boundaryMesh().size());
        List<scalarField> rotationalEBF_(mesh_.boundaryMesh().size());
        List<scalarField> rotationalDofBF_(mesh_.boundaryMesh().size());
        List<scalarField> qBF_(mesh_.boundaryMesh().size());
        List<vectorField> fDBF_(mesh_.boundaryMesh().size());
        List<scalarField> speciesRhoNIntBF_(mesh_.boundaryMesh().size());
        List<List<scalarField> > speciesRhoNBF_(typeIds_.size());
        List<List<scalarField> > vibrationalEBF_(typeIds_.size());
        forAll(typeIds_, iD)
        {
            speciesRhoNBF_[iD].setSize(mesh_.boundaryMesh().size());
            vibrationalEBF_[iD].setSize(mesh_.boundaryMesh().size());
        }

        forAll(rhoNBF_, i)
        {
            label sizePatch = mesh_.boundaryMesh()[i].size();
            rhoNBF_[i].setSize(sizePatch, 0.0);
            rhoMBF_[i].setSize(sizePatch, 0.0);
            linearKEBF_[i].setSize(sizePatch, 0.0);
            momentumBF_[i].setSize(sizePatch, vector::zero);
            rotationalEBF_[i].setSize(sizePatch, 0.0);
            rotationalDofBF_[i].setSize(sizePatch, 0.0);
            qBF_[i].setSize(sizePatch, 0.0);
            fDBF_[i].setSize(sizePatch, vector::zero);
            speciesRhoNIntBF_[i].setSize(sizePatch, 0.0);
            forAll(typeIds_, iD)
            {
                speciesRhoNBF_[iD][i].setSize(sizePatch, 0.0);
                vibrationalEBF_[iD][i].setSize(sizePatch, 0.0);
            }
        }
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
                            += cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                        rhoMBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                        linearKEBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                        momentumBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                        rotationalEBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                        rotationalDofBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                        qBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                        fDBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
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
                            += cloud_.boundaryFluxMeasurements().rhoNBF()[typeIds_[i]][j][k];
                        vibrationalEBF_[i][j][k]
                            += cloud_.boundaryFluxMeasurements().vibrationalEBF()[typeIds_[i]][j][k];
                    }
                }
            }
        }

        forAll(speciesRhoNIntBF_, j)
        {
            forAll(speciesRhoNIntBF_[j], k)
            {
                speciesRhoNIntBF_[j][k]
                    += cloud_.boundaryFluxMeasurements().rhoNIntBF()[j][k];
            }
        }

        label zoneID = mesh_.cellZones().findZoneID(zoneToModify_);
        const cellZone& zone = mesh_.cellZones()[zoneID];

        scalar k_ = physicoChemical::k.value();

        forAll(zone, idx)
        {
            const label cell = zone[idx];

            translationalT_[cell] *= factor1;
            rotationalT_[cell] *= factor1;
            UMean_[cell] *= factor1;
            forAll(typeIds_, iD)
            {
                rhoNType_[iD][cell] *= factor1;
                vibrationalTs_[iD][cell] *= factor1;

                scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
                rhoNType_[iD][cell] += relaxFactor * (nParcels_[iD][cell]
                    * cloud_.nParticle() / mesh_.V()[cell]);
                if(vibrationalETotal_[iD][cell] > VSMALL && nParcels_[iD][cell] > VSMALL)
                {
                    scalar iMean = vibrationalETotal_[iD][cell] / k_ / thetaV
                        / nParcels_[iD][cell];
                    vibrationalTs_[iD][cell] += relaxFactor * thetaV / log(1.0 + (1.0
                        / iMean));
                }
            }

            if(rhoNMean_[cell] > VSMALL)
            {
                vector UMeanLocal = momentumMean_[cell] / rhoMMean_[cell];
                UMean_[cell] += relaxFactor * UMeanLocal;

                scalar translationalTLocal = 1.0 / (3.0 * k_
                    * rhoNMean_[cell]) * (linearKEMean_[cell] - rhoMMean_[cell]
                    * (UMeanLocal & UMeanLocal));
                translationalT_[cell] += relaxFactor * translationalTLocal;
            }
            if(rotationalDofMean_[cell] > VSMALL)
            {
                rotationalT_[cell] += relaxFactor * (2.0 / k_)
                    * (rotationalEMean_[cell] / rotationalDofMean_[cell]);
            }
        }

        if(nTimeSteps_ % itBeforeRelaxing_ == 0)
        {
            Info << "Reaching Subrelaxing Step" << endl;
            scalar newFactor1 = 1.0 / (1.0 - pow(1.0 - relaxFactor, itBeforeRelaxing_));
            scalar newFactor2 = -newFactor1 * pow(1.0 - relaxFactor, itBeforeRelaxing_);

            forAll(zone, idx)
            {
                const label cell = zone[idx];

                translationalT_[cell] = newFactor1 * translationalT_[cell]
                    + newFactor2 * translationalT0_[cell];
                rotationalT_[cell] = newFactor1 * rotationalT_[cell] + newFactor2
                    * rotationalT0_[cell];
                UMean_[cell] = newFactor1 * UMean_[cell] + newFactor2 * UMean_[cell];

                forAll(typeIds_, iD)
                {
                    rhoNType_[iD][cell] = newFactor1 * rhoNType_[iD][cell]
                        + newFactor2 * rhoNType0_[iD][cell];
                    vibrationalTs_[iD][cell] = newFactor1
                        * vibrationalTs_[iD][cell] + newFactor2
                        * vibrationalTs0_[iD][cell];
                }
            }

            translationalT0_ = translationalT_;
            rotationalT0_ = rotationalT_;
            UMean0_ = UMean_;
            forAll(typeIds_, iD)
            {
                rhoNType0_[iD] = rhoNType_[iD];
                vibrationalTs0_[iD] = vibrationalTs_[iD];
            }
        }

        dimensionedScalar massDim("massDim", dimMass, 1.0);

        rhoN_ *= 0.0;
        rhoM_ *= 0.0;
        forAll(typeIds_, iD)
        {
            rhoN_ += rhoNType_[iD];
            rhoM_ += rhoNType_[iD] * cloud_.constProps(typeIds_[iD]).mass()
                * massDim;
        }
        p_ = rhoN_ * physicoChemical::k * translationalT_;

/*        if(time_.time().outputTime())
        {
            nTimeSteps_ -= itBeforeRelaxing_;
        }*/
    }
    else
    {
        nTimeSteps_++;

        forAllConstIter(dsmcCloud, cloud_, iter)
        {
            const dsmcParcel& p = iter();
            label iD = findIndex(typeIds_, p.typeId());

            if(iD != -1)
            {
                const label& cell = p.cell();
                const scalar& mass = cloud_.constProps(p.typeId()).mass();
                const scalar& rotationalDof = cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();

                rhoNMean_[cell] += 1.0;
                rhoMMean_[cell] += mass;
                linearKEMean_[cell] += mass * (p.U() & p.U());
                momentumMean_[cell] += mass * p.U();
                rotationalEMean_[cell] += p.ERot();
                rotationalDofMean_[cell] += rotationalDof; 
                vibrationalETotal_[iD][cell] += p.EVib();
                nParcels_[iD][cell] += 1.0;
            
                muu_[cell] += mass * sqr(p.U().x());
                muv_[cell] += mass * ( (p.U().x()) * (p.U().y()) );
                muw_[cell] += mass * ( (p.U().x()) * (p.U().z()) );
                mvv_[cell] += mass * sqr(p.U().y());
                mvw_[cell] += mass * ( (p.U().y()) * (p.U().z()) );
                mww_[cell] += mass * sqr(p.U().z());
            
                mcc_[cell] += mass * mag(p.U()) * mag(p.U());
                mccu_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().x());
                mccv_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().y());
                mccw_[cell] += mass * mag(p.U()) * mag(p.U()) * (p.U().z());
            
                eu_[cell] += ( p.ERot() + p.EVib() ) * (p.U().x());
                ev_[cell] += ( p.ERot() + p.EVib() ) * (p.U().y());
                ew_[cell] += ( p.ERot() + p.EVib() ) * (p.U().z());
                e_[cell] += ( p.ERot() + p.EVib() );

                if(rotationalDof > VSMALL)
                {
                    rhoNMeanInt_[cell] += 1.0;
                }
            }
        }

        List<scalarField> rhoNBF_(mesh_.boundaryMesh().size());
        List<scalarField> rhoMBF_(mesh_.boundaryMesh().size());
        List<scalarField> linearKEBF_(mesh_.boundaryMesh().size());
        List<vectorField> momentumBF_(mesh_.boundaryMesh().size());
        List<scalarField> rotationalEBF_(mesh_.boundaryMesh().size());
        List<scalarField> rotationalDofBF_(mesh_.boundaryMesh().size());
        List<scalarField> qBF_(mesh_.boundaryMesh().size());
        List<vectorField> fDBF_(mesh_.boundaryMesh().size());
        List<scalarField> speciesRhoNIntBF_(mesh_.boundaryMesh().size());
        List<List<scalarField> > speciesRhoNBF_(typeIds_.size());
        List<List<scalarField> > vibrationalEBF_(typeIds_.size());
        forAll(typeIds_, iD)
        {
            speciesRhoNBF_[iD].setSize(mesh_.boundaryMesh().size());
            vibrationalEBF_[iD].setSize(mesh_.boundaryMesh().size());
        }

        forAll(rhoNBF_, i)
        {
            label sizePatch = mesh_.boundaryMesh()[i].size();
            rhoNBF_[i].setSize(sizePatch, 0.0);
            rhoMBF_[i].setSize(sizePatch, 0.0);
            linearKEBF_[i].setSize(sizePatch, 0.0);
            momentumBF_[i].setSize(sizePatch, vector::zero);
            rotationalEBF_[i].setSize(sizePatch, 0.0);
            rotationalDofBF_[i].setSize(sizePatch, 0.0);
            qBF_[i].setSize(sizePatch, 0.0);
            fDBF_[i].setSize(sizePatch, vector::zero);
            speciesRhoNIntBF_[i].setSize(sizePatch, 0.0);
            forAll(typeIds_, iD)
            {
                speciesRhoNBF_[iD][i].setSize(sizePatch, 0.0);
                vibrationalEBF_[iD][i].setSize(sizePatch, 0.0);
            }
        }
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
                            += cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                        rhoMBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                        linearKEBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                        momentumBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                        rotationalEBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                        rotationalDofBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                        qBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                        fDBF_[j][k]
                            += cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
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
                            += cloud_.boundaryFluxMeasurements().rhoNBF()[typeIds_[i]][j][k];
                        vibrationalEBF_[i][j][k]
                            += cloud_.boundaryFluxMeasurements().vibrationalEBF()[typeIds_[i]][j][k];
                    }
                }
            }
        }

        forAll(speciesRhoNIntBF_, j)
        {
            forAll(speciesRhoNIntBF_[j], k)
            {
                speciesRhoNIntBF_[j][k]
                    += cloud_.boundaryFluxMeasurements().rhoNIntBF()[j][k];
            }
        }


        if(time_.time().outputTime())
        {
            label zoneID = mesh_.cellZones().findZoneID(zoneToModify_);
            const cellZone& zone = mesh_.cellZones()[zoneID];

            scalar k_ = physicoChemical::k.value();

            forAll(zone, idx)
            {
                const label cell = zone[idx];

                forAll(typeIds_, iD)
                {
                    scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
                    rhoNType_[iD][cell] = nParcels_[iD][cell]
                        * cloud_.nParticle() / mesh_.V()[cell] / nTimeSteps_;
                    if(vibrationalETotal_[iD][cell] > VSMALL
                        && nParcels_[iD][cell] > VSMALL)
                    {
                        scalar iMean = vibrationalETotal_[iD][cell] / k_ /
                            thetaV / nParcels_[iD][cell];
                        vibrationalTs_[iD][cell] = thetaV / log(1.0 + (1.0
                            / iMean));
                    }
                }

                if(rhoNMean_[cell] > VSMALL)
                {
                    UMean_[cell] = momentumMean_[cell] / rhoMMean_[cell];
                    translationalT_[cell] = 1.0 / (3.0 * k_ * rhoNMean_[cell])
                        * (linearKEMean_[cell] - rhoMMean_[cell] * (UMean_[cell]
                        & UMean_[cell]));
                }
                if(rotationalDofMean_[cell] > VSMALL)
                {
                    rotationalT_[cell] =(2.0 / k_) * (rotationalEMean_[cell]
                        / rotationalDofMean_[cell]);
                }
            }

            dimensionedScalar massDim("massDim", dimMass, 1.0);

            rhoN_ *= 0.0;
            rhoM_ *= 0.0;
            forAll(typeIds_, iD)
            {
                rhoN_ += rhoNType_[iD];
                rhoM_ += rhoNType_[iD] * cloud_.constProps(typeIds_[iD]).mass()
                    * massDim;
            }

            p_ = rhoN_ * physicoChemical::k * translationalT_;

        }
    }
}

void hybridSubrelax::setStartTime
(
    fileName& startTimeName
)
{
    startingTime_ = startTimeName;
}

void hybridSubrelax::setSubReset
(
    bool& subReset
)
{
    subResetting_ = subReset;
}

void hybridSubrelax::setSubrelax
(
    bool& subrelax
)
{
    subrelaxation_ = subrelax;
}

void hybridSubrelax::setCellZone
(
    word& zoneToModify
)
{
    zoneToModify_ = zoneToModify;
}

void hybridSubrelax::getFields
(
    label& nIter,
    volScalarField& rho,
    volScalarField& Ttr,
    volVectorField& U,
    volScalarField& Tt,
    volScalarField& Tr,
    PtrList<volScalarField>& Tv,
    PtrList<volScalarField>& Y
)
{
    dimensionedScalar massDim("massDim", dimMass, 1.0);
    dimensionedScalar volDim("volDim", dimVolume, 1.0);

    assignFields(rho.mesh(), rho, rhoM_);
    assignFields(U.mesh(), U, UMean_);
    scalarField rotDOF(mesh_.nCells(), 0.0);
    forAll(typeIds_, iD)
    {
        rotDOF += cloud_.constProps(typeIds_[iD]).rotationalDegreesOfFreedom()
            * rhoNType_[iD].internalField();
        assignFields(Tv[iD].mesh(), Tv[iD], vibrationalTs_[iD]);
        volScalarField dummyField = cloud_.constProps(typeIds_[iD]).mass() * rhoNType_[iD]
            * volDim; // dimless * (1 / dimVolume) * dimVolume
        assignFields(Y[iD].mesh(), Y[iD], dummyField);
    }
    forAll(mesh_.cells(), cell)
    {
        if(rhoN_[cell] > VSMALL)
        {
            rotDOF[cell] /= rhoN_[cell];
            forAll(typeIds_, iD)
            {
                Y[iD][cell] /= rhoM_[cell]; // cell values are scalar
            }
        }
    }
    scalarField dummyTField = (3.0 * translationalT_.internalField() + rotDOF
        * rotationalT_.internalField()) / (3.0 + rotDOF);
    assignFields(Ttr.mesh(), Ttr.internalField(), dummyTField);
    // One Specie Only
    assignFields(Tt.mesh(), Tt, translationalT_);
    assignFields(Tr.mesh(), Tr, rotationalT_);

    nIter = nTimeSteps_;

}

void hybridSubrelax::set0Fields
(
    volScalarField& rho,
    volScalarField& Ttr,
    volVectorField& U,
    PtrList<volScalarField>& Tv,
    PtrList<volScalarField>& Y
)
{
    dimensionedScalar massDim("massDim", dimMass, 1.0);
    dimensionedScalar volDim("volDim", dimVolume, 1.0);

    rhoN_ *= 0.0;
    assignFields(rho.mesh(), rhoM_, rho);
    assignFields(U.mesh(), UMean_, U);
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
        dimensionedScalar("0.0", dimless, 0.0)
    );

    forAll(typeIds_, iD)
    {
        assignFields(rho.mesh(), vibrationalTs_[iD], Tv[iD]);
        volScalarField rhoN = Y[iD] * rho / cloud_.constProps(typeIds_[iD]).mass()
            / massDim; // dimless * (dimMass / dimVolume) / dimMass
        assignFields(rho.mesh(), rhoNType_[iD], rhoN);
        rhoN_ += rhoNType_[iD];
        rotDOF += cloud_.constProps(typeIds_[iD]).rotationalDegreesOfFreedom()
            * rhoNType_[iD] * volDim;
    }

    assignFields(Ttr.mesh(), translationalT_, Ttr);
    assignFields(Ttr.mesh(), rotationalT_, Ttr);
    forAll(mesh_.cells(), cell)
    {
        if(rotDOF[cell] < VSMALL) rotationalT_[cell] *= 0.0;
    }
}

void hybridSubrelax::setInitialFields
(
    const List<word>& molecules,
    fileName& timeToReadFrom
)
{
    volScalarField Tt0
    (
        IOobject
        (
            "translationalT_" + fieldName_,
            timeToReadFrom,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField Tr0
    (
        IOobject
        (
            "rotationalT_" + fieldName_,
            timeToReadFrom,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volVectorField U0
    (
        IOobject
        (
            "UMean_" + fieldName_,
            timeToReadFrom,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    PtrList<volScalarField> rhoNs0(typeIds_.size());
    PtrList<volScalarField> Tvs0(typeIds_.size());

    forAll(typeIds_, i)
    {
        rhoNs0.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoN_" + molecules[i],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        Tvs0.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "vibrationalT_" + molecules[i],
                    time_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    translationalT0_ = Tt0;
    rotationalT0_ = Tr0;
    UMean0_ = U0;
    forAll(typeIds_, iD)
    {
        vibrationalTs0_[iD] = Tvs0[iD];
        rhoNType0_[iD] = rhoNs0[iD];
    }

    translationalT_ = translationalT0_;
    rotationalT_ = rotationalT0_;
    UMean_ = UMean0_;
    forAll(typeIds_, iD)
    {
        vibrationalTs_[iD] = vibrationalTs0_[iD];
        rhoNType_[iD] = rhoNType0_[iD];
    }
}

void hybridSubrelax::setItNumber(label& itValue)
{
    nTimeSteps_ = itValue;
}

//- write field
void hybridSubrelax::writeField()
{
    rhoN_.write();
    rhoM_.write();
    p_.write();
    translationalT_.write();
    rotationalT_.write();
    UMean_.write();
    forAll(typeIds_, iD)
    {
        vibrationalTs_[iD].write();
        rhoNType_[iD].write();
    }
}

void hybridSubrelax::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}

// volScalarField, volVectorField, scalarField, vectorField
template<typename T>
void hybridSubrelax::assignFields
(
    const fvMesh& mesh1,
    T& field1,
    T& field2
)
{
    if(mesh1 == mesh_)
    {
//        Info << "Identical mesh" << endl;
        field1 = field2;
    }
    else
    {
//        Info << "Non-identical mesh" << endl;
        // assuming still same mesh
        forAll(mesh1.cells(), cell)
        {
            field1[cell] = field2[cell];
        }
    }
}

} // End namespace Foam

// ************************************************************************* //

