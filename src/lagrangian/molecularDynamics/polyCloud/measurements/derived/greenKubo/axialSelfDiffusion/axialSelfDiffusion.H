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
    axialSelfDiffusion


Description
    Used to determine self diffusion coefficient inside a carbon nanotube. 


SourceFiles
    axialSelfDiffusionI.H
    axialSelfDiffusion.C
    axialSelfDiffusionIO.C

\*---------------------------------------------------------------------------*/

#ifndef axialSelfDiffusion_H
#define axialSelfDiffusion_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "boundedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class axialSelfDiffusion Declaration
\*---------------------------------------------------------------------------*/

class axialSelfDiffusion
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        List<label> molIds_;        

        vector unitVector_;

        // counter for one acf
        label nS_;        
        
        // number of MD time-steps for generating one acf
        label nSteps_;

        // outer counter to average over all acfs
        scalar nBatch_;          
        
        // just the liquid mols
        label nLiquidMols_;
        
        label nTNs_;
        
        List<label> tNaddress_;
        
//         label nMols_;
        
        // inner list = number of liquid mols
        // outer list = time (nSteps_)        
        List<List<scalar> > velocities_;
        
        // outer list = number of liquid mols
        // inner list = nSteps_
        // values = accumulation sum
        List< List<scalar> > acf_;
        
//         List<label> mols_;
        
        // self diffusion
        scalar D_;

        // mean ACF - for outputting (list = number of steps)
        List<scalar> ACF_;
        
    // Private Member Functions
        
        void setVelocities();
        
//         scalar getIntegral();

        
public:

    //- Runtime type information
    TypeName("axialSelfDiffusion");

    // Constructors

        //- Construct from components
        axialSelfDiffusion
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~axialSelfDiffusion();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        const propertyField& fields() const;

        //- write field
        void writeField();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
