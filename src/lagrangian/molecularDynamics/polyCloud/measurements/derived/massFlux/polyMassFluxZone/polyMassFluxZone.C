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

#include "polyMassFluxZone.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMassFluxZone, 0);

addToRunTimeSelectionTable(polyField, polyMassFluxZone, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void polyMassFluxZone::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");
    bb.resetBoundedBox(startPoint, endPoint);
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMassFluxZone::polyMassFluxZone
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    totalVolume_(0.0),
    length_(readScalar(propsDict_.lookup("length"))),
    unitVector_(propsDict_.lookup("unitVector")),
    useBoundBox_(false),
    molIds_(),
    massFlux_()
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyMassFluxZone::polyMassFluxZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();


    //-set the total volume
    const labelList& cells = cellZones[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        totalVolume_ += mesh_.cellVolumes()[cellI];
    }

    if (Pstream::parRun())
    {
        reduce(totalVolume_, sumOp<scalar>());
    }
    
    if (propsDict_.found("useBoundBox"))
    {
        useBoundBox_ = Switch(propsDict_.lookup("useBoundBox"));
        
        if(useBoundBox_)
        {
            setBoundBox(propsDict_, bb_, "samplingRegion");
        }
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMassFluxZone::~polyMassFluxZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMassFluxZone::createField()
{}

void polyMassFluxZone::calculateField()
{
    
    const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();
    
    const labelList& cells = mesh_.cellZones()[regionId_];
    
    vector mom = vector::zero;
    
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];
        
        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];
            
            if(findIndex(molIds_, molI->id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(molI->id());
                
                if(useBoundBox_)
                {
                    if(bb_.contains(molI->position()))
                    {
                        mom += molI->v()*massI;
                    }
                }
                else
                {
                    mom += molI->v()*massI;
                }
            }
        }
    }
    
    if(Pstream::parRun())
    {
        reduce(mom, sumOp<vector>());
    }


    scalar massFlux =(mom & unitVector_)/(length_);
    massFlux_.append(massFlux);

}

void polyMassFluxZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

        if(Pstream::master())
        {

            massFlux_.shrink();
            scalarField timeField (massFlux_.size(), 0.0);
            scalarField massFlux (massFlux_.size(), 0.0);
            
            massFlux.transfer(massFlux_);
            massFlux_.clear();
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "fluxZone_"+regionName_+"_"+fieldName_+"_M.xy",
                timeField,
                massFlux,
                true
            );


            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    casePath_,
                    "fluxZone_"+regionName_+"_"+fieldName_+"_M_SI.xy",
                    timeField*rU.refTime(),
                    massFlux*rU.refMassFlux(),
                    true
                );
            }
        }
    }
}

void polyMassFluxZone::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void polyMassFluxZone::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& polyMassFluxZone::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
