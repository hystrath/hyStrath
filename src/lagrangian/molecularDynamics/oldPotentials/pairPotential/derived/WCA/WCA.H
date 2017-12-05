/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::pairPotentials::WCA

Description


SourceFiles
    WCA.C

\*---------------------------------------------------------------------------*/

#ifndef WCA_H
#define WCA_H

#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

/*---------------------------------------------------------------------------*\
                         Class WCA Declaration
\*---------------------------------------------------------------------------*/

class WCA
:
    public pairPotential
{
    // Private data

        dictionary WCACoeffs_;

        scalar sigma_;
        scalar epsilon_;


public:

    //- Runtime type information
    TypeName("WCA");


    // Constructors

        //- Construct from components
    	WCA
		(
			const word& name,
			const reducedUnits& rU,
			const dictionary& pairPotentialProperties
		);


    //- Destructor
    ~WCA()
    {}


    // Member Functions

        scalar unscaledEnergy(const scalar r) const;

        //- Read dictionary
		bool read(const dictionary& pairPotentialProperties, const reducedUnits& rU);

		const dictionary& dict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
