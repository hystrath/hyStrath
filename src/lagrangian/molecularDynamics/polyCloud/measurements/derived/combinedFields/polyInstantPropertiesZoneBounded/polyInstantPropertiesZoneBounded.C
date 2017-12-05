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

#include "polyInstantPropertiesZoneBounded.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyInstantPropertiesZoneBounded, 0);

addToRunTimeSelectionTable(polyField, polyInstantPropertiesZoneBounded, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyInstantPropertiesZoneBounded::setBoundBoxes()
{
 
    PtrList<entry> boxList(propsDict_.lookup("boxes"));

    boxes_.setSize(boxList.size());

    forAll(boxList, b)
    {
        const entry& boxI = boxList[b];
        const dictionary& dict = boxI.dict();

        vector startPoint = dict.lookup("startPoint");
        vector endPoint = dict.lookup("endPoint");
        boxes_[b].resetBoundedBox(startPoint, endPoint);
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyInstantPropertiesZoneBounded::polyInstantPropertiesZoneBounded
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
    boxes_(),
    molIds_()
{

        // build bound boxes
    setBoundBoxes();

    // choose molecule ids to sample
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    velocityField_.clear();
//     forceField_.clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyInstantPropertiesZoneBounded::~polyInstantPropertiesZoneBounded()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyInstantPropertiesZoneBounded::createField()
{
 
}

void polyInstantPropertiesZoneBounded::calculateField()
{
    scalar mols = 0.0;
    vector vel = vector::zero;
//     vector frc = vector::zero;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            forAll(boxes_, b)
            {
                if(boxes_[b].contains(mol().position()))
                {
//                     const polyMolecule::constantProperties& constProp = molCloud_.constProps(mol().id());
                    mols += 1.0;
                    vel += mol().v();
//                     frc += molCloud_.cP().mass(mol().id())*mol().a();
                }
            }
        }
    }    
    
    // - parallel processing
    if(Pstream::parRun())
    {
        reduce(mols, sumOp<scalar>());
        reduce(vel, sumOp<vector>());
//         reduce(frc, sumOp<vector>());
    }
    
    vector velocity = vector::zero;
//     vector force = vector::zero;
    
    if(mols > 0.0)
    {
        velocity = vel/mols;
//         force = frc/mols;
    }
    
    velocityField_.append(velocity);
//     forceField_.append(force);
}


void polyInstantPropertiesZoneBounded::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            velocityField_.shrink();
//             forceField_.shrink();

            scalarField timeField (velocityField_.size(), 0.0);
            vectorField velocity (velocityField_.size(), vector::zero);
//            vectorField force (forceField_.size(), vector::zero);
            
            velocity.transfer(velocityField_);
            velocityField_.clear();
//            force.transfer(forceField_);
//            forceField_.clear();
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "zone_instant_"+fieldName_+"_velocity.xy",
                timeField,
                velocity,
                true
            );
/*            
            writeTimeData
            (
                casePath_,
                "zone_instant_"+fieldName_+"_force.xy",
                timeField,
                force,
                true
            );
*/
        }
    }
}
void polyInstantPropertiesZoneBounded::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void polyInstantPropertiesZoneBounded::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

const propertyField& polyInstantPropertiesZoneBounded::fields() const
{
    return fields_;
}

// void polyInstantPropertiesZoneBounded::updateProperties(const dictionary& newDict)
// {
//     //- the main properties should be updated first
//     updateBasicFieldProperties(newDict);

// }

} // End namespace Foam

// ************************************************************************* //
