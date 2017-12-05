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

#include "poly1DMassFluxZone.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(poly1DMassFluxZone, 0);
addToRunTimeSelectionTable(polyField, poly1DMassFluxZone, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
poly1DMassFluxZone::poly1DMassFluxZone
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
    binModel_(),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("fieldName")),
    length_(readScalar(propsDict_.lookup("length"))),
    unitVector_(propsDict_.lookup("unitVector")),
    molIds_()

{
 
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("poly1DMassFluxZone::poly1DMassFluxZone()")
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

    
    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );

    const label& nBins = binModel_->nBins();

    nBins_ = nBins;
    
    massFlowRate_.setSize(nBins);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

poly1DMassFluxZone::~poly1DMassFluxZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void poly1DMassFluxZone::createField()
{}


void poly1DMassFluxZone::calculateField()
{
    vectorField mom(nBins_, vector::zero);
    
    forAll(mesh_.cellZones()[regionId_], c)
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());
                    mom[n] += massI*molI->v();
                }
            }
        }
    }

    //- parallel communication
    if(Pstream::parRun())
    {
        forAll(mom, i)
        {
            reduce(mom[i], sumOp<vector>());
        }        
    }
    


    // collect and compute properties
    
    forAll(mom, n)
    {
        scalar massFlux = (mom[n] & unitVector_)/(length_);
        massFlowRate_[n].append(massFlux);
    }
    
}


void poly1DMassFluxZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            const reducedUnits& rU = molCloud_.redUnits();            
            scalarField bins = binModel_->binPositions();
            vectorField vectorBins = binModel_->bins();

            label nBins = nBins_;
            label nTimeSteps = massFlowRate_[0].size();
            
            for (int j = 0; j < nTimeSteps; j++)
            {
                scalarField massFlux(nBins, 0.0);
                
                forAll(massFlowRate_, i)
                {
                    massFlux[i] = massFlowRate_[i][j];
                }
                
                writeTimeData
                (
                    casePath_,
                    "bins_instant_oneDim_"+fieldName_+"_massFlowRate_SI.xy",
                    massFlux*rU.refMassFlux(),
                    "sidewaysAppend",
                    true
                );
            }
        }
        
            
        // clear fields 
        forAll(massFlowRate_, i)
        {
            massFlowRate_[i].clear();
        }
    }
}

void poly1DMassFluxZone::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void poly1DMassFluxZone::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& poly1DMassFluxZone::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
