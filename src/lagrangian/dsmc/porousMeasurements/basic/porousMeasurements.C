/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    porousMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "porousMeasurements.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(porousMeasurements, 0);
    defineRunTimeSelectionTable(porousMeasurements, fvMesh);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
Foam::porousMeasurements::porousMeasurements
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::porousMeasurements>
Foam::porousMeasurements::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
{
    word porousMeasurementModel =
        cloud.particleProperties().lookupOrDefault<word>
        (
            "porousMeasurement",
            "no"
        );
        
    porousMeasurementModel = "dsmc" 
        + static_cast<word>(std::toupper(porousMeasurementModel[0]))
        + porousMeasurementModel.substr(1) + "PorousMediumMeasurements";    
    
    Info<< "Selecting the porous measurement model:" << tab 
        << porousMeasurementModel << "\n" << endl;
    
    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(porousMeasurementModel);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "porousMeasurements::New(Time&, const polyMesh& mesh,"
            "porousMeasurements& cloud)"
        )   << "Unknown porous measurements type "
            << porousMeasurementModel << endl << endl
            << "Valid porous measurements types are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<porousMeasurements>(cstrIter()(t, mesh, cloud));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porousMeasurements::~porousMeasurements()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
