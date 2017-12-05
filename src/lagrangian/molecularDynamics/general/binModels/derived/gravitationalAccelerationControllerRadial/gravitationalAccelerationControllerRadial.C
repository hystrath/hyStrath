/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gravitationalAccelerationControllerRadial.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(gravitationalAccelerationControllerRadial, 0);

addToRunTimeSelectionTable(dsmcStateController, gravitationalAccelerationControllerRadial, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gravitationalAccelerationControllerRadial::gravitationalAccelerationControllerRadial
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    rotationalAxis_(propsDict_.lookup("rotationalAxis")),
    m_(propsDict_.lookup("rotationalPoint")),
    rVx_(propsDict_.lookup("referenceVectorX")),
    rVy_(propsDict_.lookup("referenceVectorY")),
//     radius_(readScalar(propsDict_.lookup("radius"))),
    thetaStart_(readScalar(propsDict_.lookup("angleStart"))),
    thetaEnd_(readScalar(propsDict_.lookup("angleEnd"))),
    acc_(readScalar(propsDict_.lookup("acceleration")))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    singleValueController() = true;
    
    Info << "rotationalAxis = " << rotationalAxis_ << endl;
    
    //- check
    
    rVx_ /= mag(rVx_);
    rVy_ /= mag(rVy_);
    
    // check - the x and y reference vectors should give the cross product = to the rotational axis
    
    
    vector tryA = rVx_ ^ rVy_;
    vector tryB = rVy_ ^ rVx_;
    
    Info << "check on rVx x rVy: " << tryA << nl 
         << ", check on rVy x rVx: " << tryB << endl;

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gravitationalAccelerationControllerRadial::~gravitationalAccelerationControllerRadial()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gravitationalAccelerationControllerRadial::initialConfiguration()
{
}

void gravitationalAccelerationControllerRadial::calculateProperties()
{
}

void gravitationalAccelerationControllerRadial::controlParcelsBeforeMove()
{
    Info<< "DEBUG BeforeMove" << nl
        << "    Check acceleration              = "
        << acc_ << nl
        << "    Delta time                      = "
        << mesh_.time().deltaTValue() << nl
        << endl;
        
    forAll(controlZone(), c)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            const vector& rI = p->position();

            vector rIm = rI - m_;
            vector rImMod = (rVy_ & rIm)*rVy_ + (rVx_ & rIm)*rVx_;
            scalar theta = acos((rVy_ & rImMod)/mag(rImMod));
            scalar sign = rVx_ & rImMod;
            
            if(sign < 0.0)
            {
                theta = 2.0*constant::mathematical::pi - theta;
            }
            
            if((theta >= thetaStart_ ) && (theta <= thetaEnd_))
            {
                vector normal = rotationalAxis_ ^ (rImMod /mag(rImMod));
                vector acceleration = acc_*normal;
                
                p->U() += 0.5*acceleration*mesh_.time().deltaTValue();
            }
        }
    }
}


void gravitationalAccelerationControllerRadial::output
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

void gravitationalAccelerationControllerRadial::controlParcelsBeforeCollisions()
{
    Info<< "DEBUG BeforeCollisions" << nl
        << "    Check acceleration              = "
        << acc_ << nl
        << "    Delta time                      = "
        << mesh_.time().deltaTValue() << nl
        << endl;
        
    forAll(controlZone(), c)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
        
        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            const vector& rI = p->position();
            
            vector rIm = rI - m_;
            vector rImMod = (rVy_ & rIm)*rVy_ + (rVx_ & rIm)*rVx_;
            scalar theta = acos((rVy_ & rImMod)/mag(rImMod));
            scalar sign = rVx_ & rImMod;
            
            if(sign < 0.0)
            {
                theta = 2.0*constant::mathematical::pi - theta;
            }
            
            if((theta >= thetaStart_ ) && (theta <= thetaEnd_))
            {
                vector normal = rotationalAxis_ ^ (rImMod /mag(rImMod));
                vector acceleration = acc_*normal;
                
                Info << "particle: " << p << ", acceleration: " << acceleration << endl;
                
                p->U() -= 0.5*acceleration*mesh_.time().deltaTValue();
            }
        }
    }
}

void gravitationalAccelerationControllerRadial::controlParcelsAfterCollisions()
{
    Info<< "DEBUG AfterCollisions" << nl
        << "    Check acceleration              = "
        << acc_ << nl
        << "    Delta time                      = "
        << mesh_.time().deltaTValue() << nl
        << endl;
        
    forAll(controlZone(), c)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            const vector& rI = p->position();
            
            vector rIm = rI - m_;
            vector rImMod = (rVy_ & rIm)*rVy_ + (rVx_ & rIm)*rVx_;
            scalar theta = acos((rVy_ & rImMod)/mag(rImMod));
            scalar sign = rVx_ & rImMod;
            
            if(sign < 0.0)
            {
                theta = 2.0*constant::mathematical::pi - theta;
            }
            
            if((theta >= thetaStart_ ) && (theta <= thetaEnd_))
            {
                vector normal = rotationalAxis_ ^ (rImMod /mag(rImMod));
                vector acceleration = acc_*normal;            
                p->U() += acceleration*mesh_.time().deltaTValue();
            }
        }
    }
}

void gravitationalAccelerationControllerRadial::updateProperties(const dictionary& newDict)
{
    updateStateControllerProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    setProperties();
}

void gravitationalAccelerationControllerRadial::setProperties()
{
    acc_ = readScalar(propsDict_.lookup("acceleration"));
}


}

// ************************************************************************* //
