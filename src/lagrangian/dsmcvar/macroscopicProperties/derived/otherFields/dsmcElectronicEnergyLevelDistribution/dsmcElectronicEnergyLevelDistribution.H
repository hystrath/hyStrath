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
    dsmcElectronicEnergyLevelDistribution

Description


SourceFiles
    dsmcElectronicEnergyLevelDistributionI.H
    dsmcElectronicEnergyLevelDistribution.C
    dsmcElectronicEnergyLevelDistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcElectronicEnergyLevelDistribution_H
#define dsmcElectronicEnergyLevelDistribution_H

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
                           Class dsmcElectronicEnergyLevelDistribution Declaration
\*---------------------------------------------------------------------------*/

class dsmcElectronicEnergyLevelDistribution
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;

        word fieldName_;

        label typeId_;

        word regionName_;
        label regionId_;
	
        scalar binWidth_;
        distribution distr_;



    // Private Member Functions


public:

    //- Runtime type information
    TypeName("dsmcElectronicEnergyLevelDistribution");

    // Constructors

        //- Construct from components
        dsmcElectronicEnergyLevelDistribution
	(
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
	);

    // Destructor

	~dsmcElectronicEnergyLevelDistribution();


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
