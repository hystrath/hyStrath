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
    polyMovingWallRitos

Description
    Model that applies a gravitational force on all molecules.

SourceFiles
    polyMovingWallRitosI.H
    polyMovingWallRitos.C
    polyMovingWallRitosIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyMovingWallRitos_H
#define polyMovingWallRitos_H

#include "polyStateController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "wallMotion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyMovingWallRitos Declaration
\*---------------------------------------------------------------------------*/

class polyMovingWallRitos
:
	public polyStateController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        autoPtr<wallMotion> wallMotionModel_;

        //- ids of mols
        List<label> molIds_;

        vector velocity_;
        
    // Private Member Functions

        void readProperties();

public:

    //- Runtime type information
    TypeName("polyMovingWallRitos");

    // Constructors

        //- Construct from components
        polyMovingWallRitos
		(
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyMovingWallRitos();


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
