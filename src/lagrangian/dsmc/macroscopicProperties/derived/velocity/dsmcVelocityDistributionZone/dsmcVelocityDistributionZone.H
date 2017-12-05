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
    dsmcVelocityDistributionZone

Description


SourceFiles
    dsmcVelocityDistributionZoneI.H
    dsmcVelocityDistributionZone.C
    dsmcVelocityDistributionZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcVelocityDistributionZone_H
#define dsmcVelocityDistributionZone_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcVelocityDistributionZone Declaration
\*---------------------------------------------------------------------------*/

class dsmcVelocityDistributionZone
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;

//         propertyField fields_;

        word fieldName_;

        labelField typeIds_;
//         bool oneSpecie_;

        dictionary timeDictVel_;
        timeDataMeas timeVel_;

        word regionName_;
        label regionId_;

        vector UMean_;
        vector Ucollected_;
        label nParcels_;

        scalar binWidth_;
        distribution distrX_;
        distribution distrY_;
        distribution distrZ_;



    // Private Member Functions


public:

    //- Runtime type information
    TypeName("dsmcVelocityDistributionZone");

    // Constructors

        //- Construct from components
        dsmcVelocityDistributionZone
		(
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~dsmcVelocityDistributionZone();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        //- return field
//         const volScalarField& densityField() const;

        //- write field
        void writeField();

//         const propertyField& fields() const;

        void updateProperties(const dictionary& dict);

}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
