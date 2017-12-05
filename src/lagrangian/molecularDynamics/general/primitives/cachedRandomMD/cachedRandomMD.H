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
    Foam::cachedRandomMD

Description
    Cached random number generator for use in molecular dynamics code.

    Pre-computes and caches samples on construction, so that when sample01()
    is called, the function simply returns the next (pre-computed) sample. On
    reaching the last sample, the sample sequence is repeated.

    Constructed using a seed and sample count. If the supplied count is
    negative, no caching is performed, and a new sample is generated on each
    call to sample01().

    Note: the copy constructor cannot be used if count = -1.

    Functionality extended beyond the basic Foam::cachedRandom by allowing the
    GaussNormal() and vector01() functions to be called on the cached values in
    the same manner as Foam::Random.

SourceFiles
    cachedRandomMDI.H
    cachedRandomMD.C
    cachedRandomMDTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cachedRandomMD_H
#define cachedRandomMD_H

#include "vector.H"
#include "scalarList.H"
#include <fpu_control.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cachedRandomMD;

/*---------------------------------------------------------------------------*\
                       Class cachedRandomMD Declaration
\*---------------------------------------------------------------------------*/

class cachedRandomMD
{
private:

    // Private data

        //- Initial random number seed
        label seed_;

        //- List of scalar samples
        scalarList samples_;

        //- Current sample marker
        label sampleI_;

        //- Current cache marker
		label cacheI_;

        //- Cache size multiplier
		label cacheSizeMult_;


    // Private Member Functions

        //- Returns the current sample
        inline scalar scalar01();

        //- vector with every component scalar01
        inline vector vector01();

        //- Return a normal Gaussian randon number
		//  with zero mean and unity variance N(0, 1)
		inline scalar GaussNormal();

	    //- Repopulate the value pool based on the initial seed
		inline void repopulate();

public:


    // Constructors

        //- Construct given seed and sample count
        cachedRandomMD(const label seed, const label cacheSizeMult);

        //- Copy constructor with optional reset of sampleI
        cachedRandomMD(const cachedRandomMD& cr, const bool reset = false);

        cachedRandomMD
        (
            const label seed,
            const label cacheSizeMult,
            const label numMols 
        );
        
    // Destructor
    ~cachedRandomMD();


    // Member functions

        // Access

            //- Return const access to the initial random number seed
            inline label seed() const;

            //- Return const access to the sample size
			inline label sampleSize() const;

            //- Return const access to the list of samples
            inline const scalarList& samples() const;

            //- Return the current sample marker
            inline label sampleI() const;

            //- Return the cacheSizeMult_
            inline label cacheSizeMult() const;            
            
        // Manipulation

            //- Return non-const access to the sample marker
            inline label& sampleI();

            //- Initialise the value pool
            void initialise(const label numMols);

        // Evaluation

            //- label [lower..upper]
			label integer(const label lower, const label upper);

            //- Return a sample whose components lie in the range 0-1
            template<class Type>
            Type sample01();

            //- Return a sample between start and end
            template<class Type>
            Type position(const Type& start, const Type& end);

            //- Randomise value in the range 0-1
            template<class Type>
            void randomise01(Type& value);

            //- Return a sample vector whose components lie in the range 0-1
			template<class Type>
			Type sampleVectorMD();

			//- Return a normal Gaussian random number
			//  with zero mean and unity variance N(0, 1)
			template<class Type>
			Type GaussNormalMD();


        // Operators

            //- Assignment operator
            void operator=(const cachedRandomMD& cr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Template specialisations

template<>
Foam::label cachedRandomMD::sample01<label>();

template<>
Foam::scalar cachedRandomMD::sample01<scalar>();

template<>
Foam::label cachedRandomMD::position<label>(const label& start, const label& end);

template<>
Foam::scalar cachedRandomMD::position<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
Foam::vector cachedRandomMD::sampleVectorMD<vector>();

template<>
Foam::scalar cachedRandomMD::GaussNormalMD<scalar>();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cachedRandomMDI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cachedRandomMDTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
