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

#include "polyFieldProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyFieldProperties::polyFieldProperties
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    polyFieldPropertiesDict_
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
    fields_(),
    measurementsDuringForceComp_(),
    measurementsDuringForceCompSite_()
{}

// constructor
polyFieldProperties::polyFieldProperties
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    time_(t),
    polyFieldPropertiesDict_
    (
        IOobject
        (
            "fieldPropertiesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    fieldList_(polyFieldPropertiesDict_.lookup("polyFields")),
    fieldNames_(fieldList_.size()),
    fieldIds_(fieldList_.size()),
    fields_(fieldList_.size()),
    measurementsDuringForceComp_(),
    measurementsDuringForceCompSite_()
{
    DynamicList<label> measurementsDuringForceComp(0);
    DynamicList<label> measurementsDuringForceCompSite(0);

    if(fields_.size() > 0 )
    {
        fileName fieldPath(time_.path()/"fieldMeasurements");
    
        if (!isDir(fieldPath))
        {
            mkDir(fieldPath);
        }

        // rule: clear fieldMeasurements every time the simulation is run

        // create directory: case/fieldMeasurements/poly
        fileName polyFieldsPath(fieldPath/"poly");

        if(isDir(polyFieldsPath) )
        {
            rmDir(polyFieldsPath);
        }

        mkDir(polyFieldsPath);

        Info << "Creating fields: " << nl << endl;
    
        forAll(fields_, f)
        {
            const entry& fieldI = fieldList_[f];
            const dictionary& fieldIDict = fieldI.dict();
    
            fields_[f] = autoPtr<polyField>
            (
                polyField::New(time_, mesh, molCloud, fieldIDict)
            );
    
            fieldNames_[f] = fields_[f]->type();
            fieldIds_[f] = f;

            fields_[f]->casePath() = polyFieldsPath;

            if(fields_[f]->measureInterForces())
            {
                measurementsDuringForceComp.append(f);
            }
            if(fields_[f]->measureInterForcesSites())
            {
                measurementsDuringForceCompSite.append(f);
            }
        }
    }

    //measurementsDuringForceComp_.transfer(measurementsDuringForceComp.shrink());
    //measurementsDuringForceCompSite_.transfer(measurementsDuringForceCompSite.shrink());

    measurementsDuringForceComp_.transfer(measurementsDuringForceComp);
    measurementsDuringForceCompSite_.transfer(measurementsDuringForceCompSite);
}


polyFieldProperties::~polyFieldProperties()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void polyFieldProperties::createFields()
{
    Info << nl << "Initialising the measurement fields" << nl << endl;

    forAll(fields_, f)
    {
        fields_[f]->createField();
    }
}

void polyFieldProperties::calculateFields()
{
	Info << "Calculate fields" << endl;
 
    forAll(fields_, f)
    {
        fields_[f]->calculateField();
    }
}

void polyFieldProperties::measurementsDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{
    forAll(measurementsDuringForceComp_, n)
    {
        const label& f = measurementsDuringForceComp_[n];
        fields_[f]->measureDuringForceComputation(molI, molJ);
    }
}

void polyFieldProperties::measurementsDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{
    forAll(measurementsDuringForceCompSite_, n)
    {
        const label& f = measurementsDuringForceCompSite_[n];
        fields_[f]->measureDuringForceComputationSite(molI, molJ, sI, sJ);
    }
}

//- Note, not all fields automatically write out to disk.
void polyFieldProperties::writeFields()
{
    const Time& runTime = time_;

    fileName timePath(runTime.path()/runTime.timeName()/"uniform"/"poly");

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
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
