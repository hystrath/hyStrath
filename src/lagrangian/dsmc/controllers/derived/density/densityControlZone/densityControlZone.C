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

Controls translational temperature in bins in a user-defined zone in a mesh.
Allows for a mixture of cell-by-cell and zone treatment, improved statistics
compared to cell-by-cell and increased resolution compared to zone treatment.

\*---------------------------------------------------------------------------*/

#include "densityControlZone.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(densityControlZone, 0);

addToRunTimeSelectionTable(dsmcStateController, densityControlZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
densityControlZone::densityControlZone
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    volume_(0.0),
    randomSquare_()

{
    density_ = readScalar(propsDict_.lookup("numberDensity"));
    
    velocity_ = propsDict_.lookup("velocity");
    
    temperature_ = readScalar(propsDict_.lookup("temperature"));


	
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
            FatalErrorIn("densityControlZone")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"controllersDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    // ---------------------------------------------------
    
    if(molecules.size() != 1)
    {
        FatalErrorIn("densityControlZone")
                << "This controller allows only one typeId: " << molecules 
                << nl << "in: "
                << t.system()/"controllersDict"
                << exit(FatalError);
    }
     
    Info << "densityControlZone" << nl << endl;
    
    forAll(controlZone(), c)
    {
        const label& cellI = controlZone()[c];
        volume_ += mesh_.cellVolumes()[cellI];
    }
    
    if(Pstream::parRun())
    {
        reduce(volume_, sumOp<scalar>());
    }

    randomSquare_.setInitialConfiguration
    (
        mesh_,
        regionId_
    );
    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

densityControlZone::~densityControlZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void densityControlZone::initialConfiguration()
{
}

void densityControlZone::calculateProperties()
{

}

void densityControlZone::controlParcelsBeforeMove()
{
	
}


void densityControlZone::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

    }
}


void densityControlZone::controlParcelsBeforeCollisions()
{
    
}

void densityControlZone::controlParcelsAfterCollisions()
{
    label nParcels = 0;
    
    //1. measure number of molecules
    {
        const labelList& cells = mesh_.cellZones()[regionId_];

        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cell];

            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];

                const label& typeId = p->typeId();
                
                if(findIndex(typeIds_, typeId) != -1) 
                {
                    const scalar& RWF = cloud_.getRWF_cell(cell);
                    
                    nParcels += RWF;
                }
            }
        }
    }
    
    if(Pstream::parRun())
    {
        reduce(nParcels, sumOp<label>());
    }
    
    Info << "nParcels (measured) = " << nParcels 
        << ", measured number density = " << nParcels*cloud_.nParticle()/volume_
        << ", target number density = " << density_
        << endl;

    //2. compute no of parcels to insert
    
    label nParcelsNew = label(density_*volume_/cloud_.nParticle()) - nParcels;
    bool insert = true;


    // deleting instead?
    if(nParcelsNew < 0)
    {
        nParcelsNew = -nParcelsNew;
        insert = false;
        
        Info << "nParcels to delete: " << nParcelsNew << endl;        
    }
    else
    {
        Info << "nParcels to insert: " << nParcelsNew << endl;    
    }
    
    //3. pick random points to insert new set of parcels
    
    vectorField parcelPositionToInsert(nParcelsNew, vector::zero);

    
    if(Pstream::master())
    {
        //randomly generate points in square
        forAll(parcelPositionToInsert, p)
        {
            parcelPositionToInsert[p] = randomSquare_.randomPoint();
            parcelPositionToInsert[p].z() = 0.0;
        }
    }
    
    if(Pstream::parRun())
    {
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << parcelPositionToInsert;
                }
            }
        }
   
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                vectorField parcelPositionToInsertProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> parcelPositionToInsertProc;
                }

                parcelPositionToInsert += parcelPositionToInsertProc;
            }
        }
    }
    
    // insert parcels
    if(insert)
    {
        label nParcelsAdded = 0;
        
        forAll(parcelPositionToInsert, p)
        {
            const vector& position = parcelPositionToInsert[p];
            
            label cell = -1;
            label tetFace = -1;
            label tetPt = -1;

            mesh_.findCellFacePt
            (
                position,
                cell,
                tetFace,
                tetPt
            );

            if(cell != -1)
            {
                //insert parcel
                const label& typeId = typeIds_[0];
                
                // new velocity picked from Maxwellian distribution
                scalar mass = cloud_.constProps(typeId).mass();

                scalar sigma = sqrt(physicoChemical::k.value()*temperature_/mass);

        
                vector U
                (
                    sigma*rndGen_.GaussNormal(),
                    sigma*rndGen_.GaussNormal(),
                    sigma*rndGen_.GaussNormal()
                );
                
                U += velocity_;
                
                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    temperature_,
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                );
                
                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    temperature_,
                    cloud_.constProps(typeId).vibrationalDegreesOfFreedom(),
                    typeId
                );
                
                const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);
                
                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    temperature_,
                    cP.degeneracyList(),
                    cP.electronicEnergyList(),
                    typeId
                );
                
                const scalar& RWF = cloud_.getRWF_cell(cell);
              
                cloud_.addNewParcel
                (
                    position,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cell,
                    tetFace,
                    tetPt,
                    typeId,
                    -1,
                    0,
                    vibLevel
                );
                
                nParcelsAdded++;
            }
        }

    
        if(Pstream::parRun())
        {
            reduce(nParcelsAdded, sumOp<label>());
        }

        cloud_.reBuildCellOccupancy(); 
        
        Info << "no. of parcels added = " << nParcelsAdded << endl;
    }
    else //delete parcels
    {
        label nParcelsDeleted = 0;
        
        forAll(parcelPositionToInsert, p)
        {
            const vector& position = parcelPositionToInsert[p];
            
            label cell = -1;
            label tetFace = -1;
            label tetPt = -1;

            mesh_.findCellFacePt
            (
                position,
                cell,
                tetFace,
                tetPt
            );

            if(cell != -1)
            {
                const List<dsmcParcel*>& molsInCell = cloud_.cellOccupancy()[cell];
                label cellMolRemoveId = -1;

                scalar rClosest = GREAT;
                
                forAll(molsInCell, mIC)
                {
                    dsmcParcel* p = molsInCell[mIC];

                    const label& typeId = p->typeId();
                    
                    if(findIndex(typeIds_, typeId) != -1) 
                    {
                        scalar rD = mag(p->position() - position);
                        
                        if(rD < rClosest)
                        {
                            rClosest = rD;
                            cellMolRemoveId = mIC;
                        }
                    }
                }

                if(cellMolRemoveId != -1)
                {
                    dsmcParcel* delParcel = molsInCell[cellMolRemoveId];
                    
                    //- remove from cellOccupancy before deleting it from cloud
                    cloud_.removeParcelFromCellOccupancy(cellMolRemoveId, cell);
                    cloud_.deleteParticle(*delParcel);
                    
                    nParcelsDeleted++;
                }
            }
        }
        
        if(Pstream::parRun())
        {
            reduce(nParcelsDeleted, sumOp<label>());
        }

        Info << "no. of parcels deleted = " << nParcelsDeleted << endl;        
    }
    
    Info << "new size of cloud = " << cloud_.size() << endl;
}


void densityControlZone::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

}

void densityControlZone::setProperties()
{

}


}
// End namespace Foam

// ************************************************************************* //
