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
    dsmcDeletionPatch

Description

    Fully deletion boundary patch: dsmcParcels get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    dsmcDeletionPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcDeletionPatch_H
#define dsmcDeletionPatch_H

#include "dsmcPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dsmcDeletionPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcDeletionPatch
:
	  public dsmcPatchBoundary
{

private:

    // Private data

        dictionary propsDict_;


protected:

    // Protected data

        //label timeIndex_;

        //scalar writeInterval_;

        //label writeIntSteps_;

        //scalar startTime_;

        //scalarField deletedMassFlux_;

        bool allSpecies_;

        labelList typeIds_;


        //scalar massFlux_;


    // Protected member functions

        void setProperties();



public:

    //- Runtime type information
    TypeName("dsmcDeletionPatch");

    // Constructors

        //- Construct from components
        dsmcDeletionPatch
		    (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		    );

    // Destructor

	    ~dsmcDeletionPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td);

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
