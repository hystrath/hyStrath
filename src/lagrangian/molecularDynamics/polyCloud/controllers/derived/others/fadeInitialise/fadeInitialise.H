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
    fadeInitialise

Description
    Use FADE to initalise systems that produce overlaps.
    

SourceFiles
    fadeInitialiseI.H
    fadeInitialise.C
    fadeInitialiseIO.C

\*---------------------------------------------------------------------------*/

#ifndef fadeInitialise_H
#define fadeInitialise_H

#include "polyStateController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fadeInitialise Declaration
\*---------------------------------------------------------------------------*/

class fadeInitialise
:
	public polyStateController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;


        //- ids of mols
        List<label> molIds_;

//         scalar nSteps_;
//         scalar dF_;
//         scalar runningSum_;
        scalar deltaT_;
        scalar t_;
        scalar tauT_;
        label n_;
        
    // Private Member Functions



public:

    //- Runtime type information
    TypeName("fadeInitialise");

    // Constructors

        //- Construct from components
        fadeInitialise
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~fadeInitialise();


    // Member Functions
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

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
