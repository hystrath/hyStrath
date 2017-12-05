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
    temperatureController

Description


SourceFiles
    temperatureController.C


\*---------------------------------------------------------------------------*/

#ifndef temperatureController_H
#define temperatureController_H

#include "dsmcStateController.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class temperatureController Declaration
\*---------------------------------------------------------------------------*/

class temperatureController
:
   public dsmcStateController
{

private:
    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
        
        dictionary timeDictVel_;
        timeDataMeas timeVel_;

        //- time constant
        scalar tauT_;

        // component control

        bool componentControl_;
        bool X_;
        bool Y_;
        bool Z_;

        labelList typeIds_;

        //- measurements
        scalarField massV_;
        vectorField momV_;
        vectorField UMean_;
        scalarField mcc_;
        scalarField m_;
        scalarField nParcels_;

        scalarField measuredTranslationalTemperature_;

        //- velocity rescaling term
        scalarField chi_;


    // Private Member Functions

		void setProperties();

public:

    //- Runtime type information
    TypeName("temperatureController");

    // Constructors

        //- Construct from components
        temperatureController
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
        // Destructor
        ~temperatureController();
    
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

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
