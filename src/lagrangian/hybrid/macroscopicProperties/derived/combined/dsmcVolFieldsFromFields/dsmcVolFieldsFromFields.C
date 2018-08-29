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

#include "dsmcVolFieldsFromFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVolFieldsFromFields, 0);

addToRunTimeSelectionTable(dsmcVolFields, dsmcVolFieldsFromFields, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVolFieldsFromFields::dsmcVolFieldsFromFields
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcVolFields(t, mesh, cloud, dict),
    heatFluxTraVector_
    (
        IOobject
        (
            "heatFluxTraVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
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
    heatFluxRotVector_
    (
        IOobject
        (
            "heatFluxRotVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
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
    heatFluxVibVector_
    (
        IOobject
        (
            "heatFluxVibVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
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
    eRotU_(mesh_.nCells(), 0.0),
    eRotV_(mesh_.nCells(), 0.0),
    eRotW_(mesh_.nCells(), 0.0),
    eVibU_(mesh_.nCells(), 0.0),
    eVibV_(mesh_.nCells(), 0.0),
    eVibW_(mesh_.nCells(), 0.0),
    eVib_(mesh_.nCells(), 0.0)
{
    // initialisation
    nTimeSteps_ = propsDict_.lookupOrDefault<scalar>("nPrevTimeSteps", 0.0);

    if (nTimeSteps_ > SMALL)
    {
        scalar nTBynEq = nTimeSteps_/cloud_.nParticle();
        const scalarField& V_ = mesh_.cellVolumes();

        forAll(typeIds_, iD)
        {
            if (iD != -1)
            {
                scalarList EVib
                (
                    cloud_.constProps(typeIds_[iD])
                        .vibrationalDegreesOfFreedom()
                );
                    
                forAll(mesh_.cells(), cell)
                {
                    if (EVib.size() > 0)
                    {
                        forAll(EVib, i)
                        {
                            eVib_[cell] += vibrationalETotal_[iD][i][cell];
                        }
                    }
                }
            }
        }

        forAll(mesh_.cells(), cell)
        {
            const scalar nTVBynEq = nTBynEq*V_[cell];
            
            eRotU_[cell] = heatFluxRotVector_[cell].x() * nTVBynEq
                + rotationalEMean_[cell] * UMean_[cell].x();
            eRotV_[cell] = heatFluxRotVector_[cell].y() * nTVBynEq
                + rotationalEMean_[cell] * UMean_[cell].y();
            eRotW_[cell] = heatFluxRotVector_[cell].z() * nTVBynEq
                + rotationalEMean_[cell] * UMean_[cell].z();
                
            eVibU_[cell] = heatFluxVibVector_[cell].x() * nTVBynEq
                + eVib_[cell] * UMean_[cell].x();
            eVibV_[cell] = heatFluxVibVector_[cell].y() * nTVBynEq
                + eVib_[cell] * UMean_[cell].y();
            eVibW_[cell] = heatFluxVibVector_[cell].z() * nTVBynEq
                + eVib_[cell] * UMean_[cell].z();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVolFieldsFromFields::~dsmcVolFieldsFromFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVolFieldsFromFields::calculateField()
{  
    dsmcVolFields::calculateField();

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& p = iter();
        label iD = findIndex(typeIds_, p.typeId());

        if (iD != -1)
        {
            const label cell = p.cell();
                
            scalarList EVib
            (
                cloud_.constProps(typeIds_[iD])
                    .vibrationalDegreesOfFreedom()
            );
            
            scalar vibEn = 0.0;
            
            if (EVib.size() > 0)
            {
                forAll(EVib, i)
                {
                    EVib[i] = p.vibLevel()[i]
                         *physicoChemical::k.value()                
                         *cloud_.constProps(p.typeId()).thetaV()[i];
                    
                    vibEn += EVib[i];
                }
            }
            
            eRotU_[cell] += p.ERot()*p.U().x();
            eRotV_[cell] += p.ERot()*p.U().y();
            eRotW_[cell] += p.ERot()*p.U().z();
            eVibU_[cell] += vibEn*p.U().x();
            eVibV_[cell] += vibEn*p.U().y();
            eVibW_[cell] += vibEn*p.U().z();
            eVib_[cell] += vibEn;
        }
    }
    
    if (time_.time().outputTime())
    {
        forAll(rhoNMean_, cell)
        {
            if (rhoNMean_[cell] > SMALL)
            {
                heatFluxTraVector_[cell].x() = 0.5 * rhoN_[cell] * (mccu_[cell]
                    - mcc_[cell] * UMean_[cell].x()) / rhoNMean_[cell]
                    - pressureTensor_[cell].xx() * UMean_[cell].x()
                    - pressureTensor_[cell].xy() * UMean_[cell].y()
                    - pressureTensor_[cell].xz() * UMean_[cell].z();

                heatFluxRotVector_[cell].x() = rhoN_[cell] * (eRotU_[cell]
                    - rotationalEMean_[cell] * UMean_[cell].x()) / rhoNMean_[cell];

                heatFluxVibVector_[cell].x() = rhoN_[cell] * (eVibU_[cell]
                    - eVib_[cell] * UMean_[cell].x()) / rhoNMean_[cell];


                heatFluxVector_[cell].x() = heatFluxTraVector_[cell].x()
                    + heatFluxRotVector_[cell].x()
                    + heatFluxVibVector_[cell].x();


                heatFluxTraVector_[cell].y() = 0.5 * rhoN_[cell] * (mccu_[cell]
                    - mcc_[cell] * UMean_[cell].y()) / rhoNMean_[cell]
                    - pressureTensor_[cell].yx() * UMean_[cell].x()
                    - pressureTensor_[cell].yy() * UMean_[cell].y()
                    - pressureTensor_[cell].yz() * UMean_[cell].z();

                heatFluxRotVector_[cell].y() = rhoN_[cell] * (eRotV_[cell]
                    - rotationalEMean_[cell] * UMean_[cell].y()) / rhoNMean_[cell];

                heatFluxVibVector_[cell].y() = rhoN_[cell] * (eVibV_[cell]
                    - eVib_[cell] * UMean_[cell].y()) / rhoNMean_[cell];

                heatFluxVector_[cell].y() = heatFluxTraVector_[cell].y()
                    + heatFluxRotVector_[cell].y()
                    + heatFluxVibVector_[cell].y();


                heatFluxTraVector_[cell].z() = 0.5 * rhoN_[cell] * (mccw_[cell]
                    - mcc_[cell] * UMean_[cell].z()) / rhoNMean_[cell]
                    - pressureTensor_[cell].zx() * UMean_[cell].x()
                    - pressureTensor_[cell].zy() * UMean_[cell].y()
                    - pressureTensor_[cell].zz() * UMean_[cell].z();

                heatFluxRotVector_[cell].z() = rhoN_[cell] * (eRotW_[cell]
                    - rotationalEMean_[cell] * UMean_[cell].z()) / rhoNMean_[cell];

                heatFluxVibVector_[cell].z() = rhoN_[cell] * (eVibW_[cell]
                    - eVib_[cell] * UMean_[cell].z()) / rhoNMean_[cell];

                heatFluxVector_[cell].z() = heatFluxTraVector_[cell].z()
                    + heatFluxRotVector_[cell].z()
                    + heatFluxVibVector_[cell].z();
            }
            else
            {
                heatFluxTraVector_[cell] = vector::zero;
                heatFluxRotVector_[cell] = vector::zero;
                heatFluxVibVector_[cell] = vector::zero;
                heatFluxVector_[cell] = vector::zero;
            }
            
            heatFluxTraVector_.write();
            heatFluxRotVector_.write();
            heatFluxVibVector_.write();
        }

        //- reset
        if (time_.resetFieldsAtOutput())
        {
            forAll(rhoNMean_, celli)
            {
                eRotU_[celli] = 0.0;
                eRotV_[celli] = 0.0;
                eRotW_[celli] = 0.0;
                eVibU_[celli] = 0.0;
                eVibV_[celli] = 0.0;
                eVibW_[celli] = 0.0;
                eVib_[celli] = 0.0;
            }
        }

        if (averagingAcrossManyRuns_)
        {
            writeOut();
        }
    }
}


//- reset fields when mesh is edited
void dsmcVolFieldsFromFields::resetField()
{
    dsmcVolFields::resetField();
    
    eRotU_.clear();
    eRotV_.clear();
    eRotW_.clear();
    eVibU_.clear();
    eVibV_.clear();
    eVibW_.clear();
    eVib_.clear();
    
    eRotU_.setSize(mesh_.nCells(), 0.0);
    eRotV_.setSize(mesh_.nCells(), 0.0);
    eRotW_.setSize(mesh_.nCells(), 0.0);
    eVibU_.setSize(mesh_.nCells(), 0.0);
    eVibV_.setSize(mesh_.nCells(), 0.0);
    eVibW_.setSize(mesh_.nCells(), 0.0);
    eVib_.setSize(mesh_.nCells(), 0.0);
}


//- write field
void dsmcVolFieldsFromFields::writeField()
{
    dsmcVolFields::writeField();
}

void dsmcVolFieldsFromFields::updateProperties(const dictionary& newDict)
{
    dsmcVolFields::updateProperties(newDict);
}

} // End namespace Foam

// ************************************************************************* //

