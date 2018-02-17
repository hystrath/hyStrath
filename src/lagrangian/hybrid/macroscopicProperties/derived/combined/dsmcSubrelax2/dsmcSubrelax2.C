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
#include "dsmcSubrelax2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcSubrelax2, 0);

addToRunTimeSelectionTable(dsmcField, dsmcSubrelax2, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcSubrelax2::dsmcSubrelax2
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    rhoN0_
    (
        IOobject
        (
            "rhoN0_",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM0_
    (
        IOobject
        (
            "rhoM0_",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    ),
    p0_
    (
        IOobject
        (
            "p0_",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimPressure, 0.0)
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

    nTimeSteps_(0),
    typeIds_(),

    startingTime_(propsDict_.lookupOrDefault<fileName>("fileName", "0")),

    pName_(propsDict_.lookup("pName")),
    UName_(propsDict_.lookup("UName")),
    TtName_(propsDict_.lookup("TtName")),
    TrName_(propsDict_.lookup("TrName")),
    rhoNNames_(propsDict_.lookup("rhoNames")),
    TvNames_(propsDict_.lookup("TvNames")),
    mass_(readList<scalar>(propsDict_.lookup("masses"))),

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
            FatalErrorIn("dsmcSubrelax2::dsmcSubrelax2()")
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
                    IOobject::AUTO_WRITE
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
                    IOobject::AUTO_WRITE
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
    }

    setInitialFields(molecules, startingTime_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSubrelax2::~dsmcSubrelax2()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//- initial conditions
void dsmcSubrelax2::createField()
{
    Info << "Initialising dsmcSubrelax2 field" << endl;
}


void dsmcSubrelax2::calculateField()
{  
    nTimeSteps_++;

/*    if(nTimeSteps_ < 1000)
    {
        itBeforeRelaxing_ = 100;
    }
    else if(nTimeSteps_ > 5000)
    {
        itBeforeRelaxing_ = 1000;
    }
    else
    {
        itBeforeRelaxing_ = 500;
    }*/
    scalar relaxFactor = 1.0 / scalar(itBeforeRelaxing_);
    scalar factor1 = (1.0 - relaxFactor);

    // Start subrelaxation
    scalarField rhoNMean_(mesh_.nCells(), 0.0);
    scalarField rhoMMean_(mesh_.nCells(), 0.0);
    scalarField rhoNMeanInt_(mesh_.nCells(), 0.0);
    scalarField linearKEMean_(mesh_.nCells(), 0.0);
    vectorField momentumMean_(mesh_.nCells(), vector::zero);
    scalarField rotationalEMean_(mesh_.nCells(), 0.0);
    scalarField rotationalDofMean_(mesh_.nCells(), 0.0);
    scalarField muu_(mesh_.nCells(), 0.0);
    scalarField muv_(mesh_.nCells(), 0.0);
    scalarField muw_(mesh_.nCells(), 0.0);
    scalarField mvv_(mesh_.nCells(), 0.0);
    scalarField mvw_(mesh_.nCells(), 0.0);
    scalarField mww_(mesh_.nCells(), 0.0);
    scalarField mcc_(mesh_.nCells(), 0.0);
    scalarField mccu_(mesh_.nCells(), 0.0);
    scalarField mccv_(mesh_.nCells(), 0.0);
    scalarField mccw_(mesh_.nCells(), 0.0);
    scalarField eu_(mesh_.nCells(), 0.0);
    scalarField ev_(mesh_.nCells(), 0.0);
    scalarField ew_(mesh_.nCells(), 0.0);
    scalarField e_(mesh_.nCells(), 0.0);
    List<scalarField> vibrationalETotal_(typeIds_.size());
    List<scalarField> nParcels_(typeIds_.size());

    forAll(typeIds_, iD)
    {
        vibrationalETotal_[iD].setSize(mesh_.nCells(), 0.0);
        nParcels_[iD].setSize(mesh_.nCells(), 0.0);
    }

    volScalarField dummyVal
    (
        IOobject
        (
            "cellValue",
            "0",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

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

    scalar k_ = physicoChemical::k.value();
    rhoN_ *= factor1;
    rhoM_ *= factor1;
    translationalT_ *= factor1;
    rotationalT_ *= factor1;
    p_ *= factor1;
    UMean_ *= factor1;
    forAll(typeIds_, iD)
    {
        rhoNType_[iD] *= factor1;
        vibrationalTs_[iD] *= factor1;
    }
    scalarField rhoNLocal = rhoNMean_ * cloud_.nParticle() / mesh_.V();
    scalarField translationalTLocal(mesh_.nCells(),0.0);
    rhoN_.internalField() += relaxFactor * rhoNLocal;
    rhoM_.internalField() += relaxFactor * rhoMMean_ * cloud_.nParticle() / mesh_.V();

    forAll(mesh_.cells(), cell)
    {
        if(rhoNMean_[cell] > VSMALL)
        {
            vector UMeanLocal = momentumMean_[cell] / rhoMMean_[cell];
            UMean_[cell] += relaxFactor * UMeanLocal;
            translationalTLocal[cell] = 1.0 / (3.0 * k_
                * rhoNMean_[cell]) * (linearKEMean_[cell] - rhoMMean_[cell]
                * (UMeanLocal & UMeanLocal));
        }
        if(rotationalDofMean_[cell] > VSMALL)
        {
            rotationalT_[cell] += relaxFactor * (2.0 / k_)
                * (rotationalEMean_[cell] / rotationalDofMean_[cell]);
        }
    }
    translationalT_.internalField() += relaxFactor * translationalTLocal;
    p_.internalField() += relaxFactor * rhoNLocal * k_
        * translationalTLocal;
    forAll(typeIds_, iD)
    {
        scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
        rhoNType_[iD].internalField() += relaxFactor * (nParcels_[iD] * cloud_.nParticle()
            / mesh_.V());
        forAll(mesh_.cells(), cell)
        {
            if(vibrationalETotal_[iD][cell] > VSMALL && nParcels_[iD][cell] > VSMALL)
            {
                scalar iMean = vibrationalETotal_[iD][cell] / k_ / thetaV
                    / nParcels_[iD][cell];
                vibrationalTs_[iD][cell] += relaxFactor * thetaV / log(1.0 + (1.0
                    / iMean));
            }
        }
    }


/*    forAll(mesh_.cells(), cell)
    {
        if(rhoNMean_[cell] > VSMALL)
        {
            tau_[cellI].xx() += relaxFactor * (p_[cell] - (cloud_nParticle()
                * muu_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].x()
                * UMean_[cell].x()));
            tau_[cellI].xy() += relaxFactor * (- (cloud_nParticle()
                * muv_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].x()
                * UMean_[cell].y()));
            tau_[cellI].xz() += relaxFactor * (- (cloud_nParticle()
                * muw_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].x()
                * UMean_[cell].z()));
            tau_[cellI].yx() += relaxFactor * (- (cloud_nParticle()
                * muv_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].y()
                * UMean_[cell].x()));
            tau_[cellI].yy() += relaxFactor * (p_[cell] - (cloud_nParticle()
                * mvv_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].y()
                * UMean_[cell].y()));
            tau_[cellI].yz() += relaxFactor * (- (cloud_nParticle()
                * mvw_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].y()
                * UMean_[cell].z()));
            tau_[cellI].zx() += relaxFactor * (- (cloud_nParticle()
                * muw_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].z()
                * UMean_[cell].z()));
            tau_[cellI].zy() += relaxFactor * (- (cloud_nParticle()
                * mvw_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].z()
                * UMean_[cell].y()));
            tau_[cellI].zz() += relaxFactor * (p_[cell] - (cloud_nParticle()
                * mww_[cell] / mesh_.V() - rhoM_[cell] * UMean_[cell].z()
                * UMean_[cell].z()));

            heatFluxVector_[cell].x() += relaxFactor * (cloud_.nParticle()
                * (0.5 * mccu_[cell] - 0.5 * mcc_[cell] * UMean_[cell].x()
                + eu_[cell] - e_[cell] * UMean_[cell].x())/ mesh_.V()
                + (tau_[cell].xx() - p_[cell]) * UMean_[cell].x()
                + tau_[cell].xy() * UMean_[cell].y() + tau_[cell].xz()
                * UMean_[cell].z() );
            heatFluxVector_[cell].y() += relaxFactor * (cloud_.nParticle()
                * (0.5 * mccv_[cell] - 0.5 * mcc_[cell] * UMean_[cell].y()
                + ev_[cell] - e_[cell] * UMean_[cell].y())/ mesh_.V()
                + (tau_[cell].yx() - p_[cell]) * UMean_[cell].x()
                + tau_[cell].yy() * UMean_[cell].y() + tau_[cell].yz()
                * UMean_[cell].z() );
            heatFluxVector_[cell].z() += relaxFactor * (cloud_.nParticle()
                * (0.5 * mccw_[cell] - 0.5 * mcc_[cell] * UMean_[cell].z()
                + ew_[cell] - e_[cell] * UMean_[cell].z())/ mesh_.V()
                + (tau_[cell].zx() - p_[cell]) * UMean_[cell].x()
                + tau_[cell].zy() * UMean_[cell].y() + tau_[cell].zz()
                * UMean_[cell].z() );
        }
    }

    scalarField molarCP(mesh_.nCells(), 0.0);
    scalarField molarCV(mesh_.nCells(), 0.0);
    scalarField molecularMass(mesh_.nCells(), 0.0);
    forAll(nParcels_, iD)
    {
        forAll(mesh_.cells(), cell)
        {
            molecularMass[cell] += cloud_.constProps(iD).
        }
    }*/

    if(nTimeSteps_ % itBeforeRelaxing_ == 0)
    {
        scalar newFactor1 = 1.0 / (1.0 - pow(1.0 - relaxFactor, itBeforeRelaxing_));
        scalar newFactor2 = -newFactor1 * pow(1.0 - relaxFactor, itBeforeRelaxing_);

        rhoN_ = newFactor1 * rhoN_ + newFactor2 * rhoN0_;
        rhoM_ = newFactor1 * rhoM_ + newFactor2 * rhoM0_;
        translationalT_ = newFactor1 * translationalT_ + newFactor2
            * translationalT0_;
        rotationalT_ = newFactor1 * rotationalT_ + newFactor2 * rotationalT0_;
        p_ = newFactor1 * p_ + newFactor2 * p0_;
        UMean_ = newFactor1 * UMean_ + newFactor2 * UMean_;

        forAll(typeIds_, iD)
        {
            rhoNType_[iD] = newFactor1 * rhoNType_[iD] + newFactor2
                * rhoNType0_[iD];
            vibrationalTs_[iD] = newFactor1 * vibrationalTs_[iD] + newFactor2
                * vibrationalTs0_[iD];
        }


        rhoN0_ = rhoN_;
        rhoM0_ = rhoM_;
        translationalT0_ = translationalT_;
        rotationalT0_ = rotationalT_;
        p0_ = p_;
        UMean0_ = UMean_;
        forAll(typeIds_, iD)
        {
            rhoNType0_[iD] = rhoNType_[iD];
            vibrationalTs0_[iD] = vibrationalTs_[iD];
        }
    }
    scalar dummyScalar1 = 0.0;
    scalar dummyScalar2 = 0.0;
    dummyScalar1 = max(dummyVal.internalField() * rhoNMean_);
    dummyScalar2 = max(dummyVal.internalField() * rhoN_.internalField()
        * mesh_.V() / cloud_.nParticle());
    if(Pstream::parRun())
    {
        reduce(dummyScalar1, sumOp<scalar>());
        reduce(dummyScalar2, sumOp<scalar>());
    }
    Info << "valueLocal " << dummyScalar1 << endl;
    Info << "valueAv " << dummyScalar2 << endl;

    if(time_.time().outputTime())
    {        
        if(time_.resetFieldsAtOutput() and not stopSamplingReset_)
        {
            nTimeSteps_ = 0;

            const List<word> molecules (propsDict_.lookup("typeIds"));
            setInitialFields(molecules, startingTime_);
        }
        else
        {
            startingTime_ = time_.time().timeName();
        }
    }
}

void dsmcSubrelax2::setReset
(
    bool& stopReset
)
{
    stopSamplingReset_ = stopReset;
}

void dsmcSubrelax2::setInitialFields
(
    const List<word>& molecules,
    fileName& timeToReadFrom
)
{
    dimensionedScalar massDim("massDim", dimMass, 1.0);
    dimensionedScalar rhoDim("rhoDim", dimMass / dimVolume, 1.0);

    volScalarField p0
    (
        IOobject
        (
            "p_" + fieldName_,
            timeToReadFrom,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

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
    rhoN0_ = 0.0 * rhoDim / massDim;
    rhoM0_ = 0.0 * rhoDim;
    p0_ = p0;
    translationalT0_ = Tt0;
    rotationalT0_ = Tr0;
    UMean0_ = U0;
    forAll(typeIds_, iD)
    {
        vibrationalTs0_[iD] = Tvs0[iD];
        rhoNType0_[iD] = rhoNs0[iD];
        rhoN0_ += rhoNs0[iD];
        rhoM0_ += rhoNs0[iD] * mass_[iD] * massDim;
    }

    rhoN_ = rhoN0_;
    rhoM_ = rhoM0_;
    p_ = p0_;
    translationalT_ = translationalT0_;
    rotationalT_ = rotationalT0_;
    UMean_ = UMean0_;
    forAll(typeIds_, iD)
    {
        vibrationalTs_[iD] = vibrationalTs0_[iD];
        rhoNType_[iD] = rhoNType0_[iD];
    }
}


//- write field
void dsmcSubrelax2::writeField()
{}

void dsmcSubrelax2::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //

