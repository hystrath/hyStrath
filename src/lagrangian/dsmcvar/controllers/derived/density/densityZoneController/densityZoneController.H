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

Class
    densityZoneController

Description


SourceFiles
    densityZoneController.C


\*---------------------------------------------------------------------------*/

#ifndef densityZoneController_H
#define densityZoneController_H

#include "dsmcStateController.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class densityZoneController Declaration
\*---------------------------------------------------------------------------*/

class densityZoneController
:
   public dsmcStateController
{

private:
    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
        label typeId_;

        // target mass density
        scalar targetParcelDensity_;

        // list of cells in the zone (no of parcels to control)
        labelField nMols_;
        labelField nMolsActualSum_;

        // sampling zone (type ids)
        scalar avParcelDensity_;
        label measuredParcels_;

        scalar controlVolume_;
        scalarField volumeWeighting_;

    // Private Member Functions

        void setProperties();

        void deleteParcels(const label& nMols, const label& c);

        void insertParcels(const label& nMols, const label& c);

public:

    //- Runtime type information
    TypeName("densityZoneController");

    // Constructors

        //- Construct from components
        densityZoneController
        (
	        Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
        // Destructor
        ~densityZoneController();
    
    // Member Functions
    
        void initialConfiguration();

        void controlParcelsBeforeMove();
        void controlParcelsBeforeCollisions();
        void controlParcelsAfterCollisions();
    
        //- calculate any required properties
        void calculateProperties();

        //- output data
        void output
        (
	        const fileName& fixedPathName, 
	        const fileName& timePath
        );
    
        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);
    
        //- set no of molecules to control
        void nMolsToControl();
    


        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
