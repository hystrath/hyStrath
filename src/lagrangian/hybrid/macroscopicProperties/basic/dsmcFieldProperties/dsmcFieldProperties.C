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

#include "dsmcFieldProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dsmcFieldProperties::dsmcFieldProperties
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    dsmcFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(),
    fieldNames_(),
    fieldIds_(),
    fields_()
{}


dsmcFieldProperties::dsmcFieldProperties
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    dsmcFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    fieldList_(dsmcFieldPropertiesDict_.lookup("dsmcFields")),
    fieldNames_(fieldList_.size()),
    fieldIds_(fieldList_.size()),
    fields_(fieldList_.size())
{
    if(fields_.size() > 0 )
    {
        fileName fieldPath(time_.path()/"fieldMeasurements");

        if(isDir(fieldPath) )
        {
            rmDir(fieldPath);
        }

        mkDir(fieldPath);
        
        Info << "Creating fields: " << nl << endl;
    
        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();
    
            fields_[f] = autoPtr<dsmcField>
            (
                dsmcField::New(time_, mesh, cloud, fieldIDict)
            );
    
            fieldNames_[f] = fields_[f]->type();
            fieldIds_[f] = f;
            fields_[f]->casePath() = fieldPath;
        }
    }
}


dsmcFieldProperties::~dsmcFieldProperties()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void dsmcFieldProperties::createFields()
{
    Info << nl << "Initialising the measurement fields" << nl << endl;

    forAll(fields_, f)
    {
        fields_[f]->createField();
    }
}


//- this function is to be called at the beginning of the MD time-step,
//  in order to update the time scheme used by each measurement property.
void dsmcFieldProperties::updateTimeInfo()
{

    forAll(fields_, f)
    {
        fields_[f]->updateTime();
    }
}


void dsmcFieldProperties::calculateFields()
{
	Info << "Calculate fields" << endl;
 
    forAll(fields_, f)
    {
        fields_[f]->calculateField();
    }
}


// Set the reset variable for the combined fields
void dsmcFieldProperties::setReset
(
    bool& stopReset
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridVolFields")
        {
            fields_[f]->setReset(stopReset);
        }
    }
}

void dsmcFieldProperties::setStartTime
(
    fileName& startTimeName
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax")
        {
            fields_[f]->setStartTime(startTimeName);
        }
    }
}

void dsmcFieldProperties::setSubReset
(
    bool& subreset
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax")
        {
            fields_[f]->setSubReset(subreset);
        }
    }
}


void dsmcFieldProperties::setSubrelax
(
    bool& subrelax
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax")
        {
            fields_[f]->setSubrelax(subrelax);
        }
    }
}

void dsmcFieldProperties::setCellZone
(
    word& zoneToModify
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax")
        {
            fields_[f]->setCellZone(zoneToModify);
        }
    }
}

void dsmcFieldProperties::getFields
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
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax" or fieldNames_[f] == "hybridSubrelaxTest")
        {
            fields_[f]->getFields
            (
                nIter,
                rho,
                Ttr,
                U,
                Tt,
                Tr,
                Tv,
                Y
            );
        }
    }
}

void dsmcFieldProperties::set0Fields
(
    volScalarField& rho,
    volScalarField& Ttr,
    volVectorField& U,
    PtrList<volScalarField>& Tv,
    PtrList<volScalarField>& Y
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax")
        {
            fields_[f]->set0Fields
            (
                rho,
                Ttr,
                U,
                Tv,
                Y
            );
        }
    }
}

void dsmcFieldProperties::setItNumber
(
    label& itNumber
)
{
    forAll(fields_, f)
    {
        if(fieldNames_[f] == "hybridSubrelax")
        {
            fields_[f]->setItNumber
            (
                itNumber
            );
        }
    }
}

//- Note, not all fields automatically write out to hard disc. 
void dsmcFieldProperties::writeFields()
{
    const Time& runTime = time_;

    fileName timePath(runTime.path()/runTime.timeName()/"uniform");

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            if (!isDir(timePath))
            {
                mkDir(timePath);
            }
        }
    }

    forAll(fields_, f)
    {
        fields_[f]->timePath() = timePath;
        fields_[f]->writeField();
    }
   
    if(runTime.outputTime())
    {
        //- Checking for modifications in the IOdictionary
        //  this allows for run-time tuning of any parameters.

        // NOTES: 
        // At the moment, the dictionary is forced to be re-read every write-interval, 
        // and properties within the abstract and models are re-set to. 
        // The "ideal" case is to have the code identify when the dictionary has been
        // modified, before re-reading it in again. Unfortunately the .modified() function
        //  is not working properly. 

        fieldList_.clear();
    
        fieldList_ = dsmcFieldPropertiesDict_.lookup("dsmcFields");
    
        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();

            fields_[f]->updateProperties(fieldIDict);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
