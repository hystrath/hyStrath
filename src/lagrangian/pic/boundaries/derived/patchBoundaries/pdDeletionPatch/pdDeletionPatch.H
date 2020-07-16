/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    pdDeletionPatch

Description

    Fully deletion boundary patch: pdParcels get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    pdDeletionPatch.C


\*---------------------------------------------------------------------------*/

#ifndef pdDeletionPatch_H
#define pdDeletionPatch_H

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
                           Class pdDeletionPatch Declaration
\*---------------------------------------------------------------------------*/

class pdDeletionPatch
:
	public pdPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;

        //label timeIndex_;

        //scalar writeInterval_;

        //label writeIntSteps_;

        //scalar startTime_;

        //scalarField deletedMassFlux_;

        bool allSpecies_;

        labelList typeIds_;


        //scalar massFlux_;


        void setProperties();

public:

    //- Runtime type information
    TypeName("pdDeletionPatch");

    // Constructors

        //- Construct from components
        pdDeletionPatch
		(
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud,
            const dictionary& dict
		);

    // Destructor

	~pdDeletionPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(pdParcel& p, pdParcel::trackingData& td);

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
