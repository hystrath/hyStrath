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
    pdDiffuseNeutralisingWallPatch

Description
    Wall interaction setting microscopd velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    pdDiffuseNeutralisingWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef pdDiffuseNeutralisingWallPatch_H
#define pdDiffuseNeutralisingWallPatch_H

#include "pdPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "pdCloud.H"
#include "distribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdDiffuseNeutralisingWallPatch Declaration
\*---------------------------------------------------------------------------*/

class pdDiffuseNeutralisingWallPatch
:
	public pdPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;

        //label nTotReactions_;
        //label nReactionsPerTimeStep_;
        scalar              typeElec_;
        List<label>         ionIds_;
        List<List <label> > productsToNeut_;

        //bool writeRatesToTerminal_;

        //labelList catalysisTypeIds_;
        //labelList catalysedTypeIds_;
        //bool oxygenSwitch_;
        //bool nitrogenSwitch_;

        void setProperties();



public:

    //- Runtime type information
    TypeName("pdDiffuseNeutralisingWallPatch");

    // Constructors

        //- Construct from components
        pdDiffuseNeutralisingWallPatch
		(
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~pdDiffuseNeutralisingWallPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the parcel
        void controlParticle(pdParcel& p, pdParcel::trackingData& td);

        //- diffusely reflect parcel
        void diffuseReflect(pdParcel& p, pdParcel::trackingData& td);

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

		void writeField();

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);



        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
