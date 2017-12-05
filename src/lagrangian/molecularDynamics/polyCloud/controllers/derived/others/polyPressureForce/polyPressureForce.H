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
    polyPressureForce

Description
    A specific model for applying a pressure forcing on a graphene sheet,
    to act as a piston.
    Pressure force could be ramped up gradually too.

SourceFiles
    polyPressureForceI.H
    polyPressureForce.C
    polyPressureForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyPressureForce_H
#define polyPressureForce_H

#include "polyStateController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "boundBox.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyPressureForce Declaration
\*---------------------------------------------------------------------------*/

class polyPressureForce
:
	public polyStateController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        word fileName_;        
        
        //- ids of mols
        List<label> molIds_;

        labelList trackingNumbers_;
        
        List<vector> sites_;
        
        vector d_;
        
        scalar area_;
        scalar mass_;
        
        scalar targetPressure_;
        scalar endTime_;
        scalar deltaT_;
        scalar currentTime_;
        vector force_; 
        scalar dt_;

        vector velocity_;
        vector position_;
        vector netForce_;
        vector netAcc_;
        
    // Private Member Functions

        void readProperties();
        
        void centreOfMass();
        void setVelocities();
        
public:

    //- Runtime type information
    TypeName("polyPressureForce");

    // Constructors

        //- Construct from components
        polyPressureForce
		(
            Time& t,
//             const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyPressureForce();


    // Member Functions
        
        //- create an initial configuration
        void initialConfiguration();

        //- control molecules at various stages of time-step
        void controlBeforeVelocityI();
        void controlBeforeMove();
        void controlBeforeForces();
        void controlDuringForces
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );        
        void controlAfterForces();
        void controlAfterVelocityII();
        
        //- calculate any required properties (end of time-step)
        void calculateProperties(); 
        
        //- output data (optional)
        void output
        (
            const fileName& fixedPathName, 
            const fileName& timePath
        );
        
        //- update properties from a modified dictionary (optional)
        void updateProperties(const dictionary&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
