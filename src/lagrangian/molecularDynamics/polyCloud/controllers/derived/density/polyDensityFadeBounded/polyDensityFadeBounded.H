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
    polyDensityFadeBounded

Description
    - insertion / deletion of molecules using the FADE algorithm
    - NEW



SourceFiles
    polyDensityFadeBounded.C


\*---------------------------------------------------------------------------*/

#ifndef polyDensityFadeBounded_H
#define polyDensityFadeBounded_H

#include "polyStateController.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "polyFadeII.H"
#include "distributePoints.H"
#include "boundedBox.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDensityFadeBounded Declaration
\*---------------------------------------------------------------------------*/

class polyDensityFadeBounded
:
	public polyStateController
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        //- sampling zone
        boundedBox bbsampling_;
        scalar samplingVolume_;

        //- control zone
        boundedBox bbcontrol_;        
        scalar controlVolume_;
        
        distributePoints pointsBox_;
        
        //- id of injected polyMolecule
        label molId_;
        
        // max number of molecules to insert/delete
        label maxN_;

        polyFadeII insertionScheme_;
        
        scalar molMass_;
        
        scalar density_;
        vector velocity_;
        scalar temperature_;        
        
        bool deltaN_;
        label molsToControlN_;
        label molsControlled_;
        
        // - measurements
        scalar mols_;
        scalar nMeas_;
        
        scalar rhoN_;
        
        DynamicList<vector> molPoints_;
        
        label controlTimeIndex_;
        label nControlSteps_;
        
        // distribution
        bool random_;
        bool uniform_;
        
        //temporary
        
        bool output_;
        
    // Private Member Functions

        void setBoundBox
        (
            const dictionary& propsDict,
            boundedBox& bb,
            const word& name 
        );

        void measureProperties();
        
        void averageProperties();
        
        label noOfMolsToControl();
        void nMolsToInsert(label& molsToControl);
        
public:

    //- Runtime type information
    TypeName("polyDensityFadeBounded");

    // Constructors

        //- Construct from components
        polyDensityFadeBounded
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~polyDensityFadeBounded();


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
