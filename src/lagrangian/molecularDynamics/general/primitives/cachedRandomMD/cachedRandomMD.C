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

\*---------------------------------------------------------------------------*/

#include "cachedRandomMD.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar Foam::cachedRandomMD::scalar01()
{
    if (sampleI_ < 0) //The sample iterator is uninitialised
    {
        return osRandomDouble();
    }

    if (sampleI_ == samples_.size() - 1) //The end of the cache has been reached
    {
    	scalar s = samples_[sampleI_]; //Grab the final value
    	cacheI_++; //Increment the cache counter
    	repopulate(); //Call function to repopulate the cache
        sampleI_ = 0; //Reset the sample iterator
        return s;
    }
    else //Normal operation
    {
        scalar s = samples_[sampleI_]; //Grab the value
        sampleI_++; //Increment the sample iterator
        return s;
    }
}

inline Foam::vector Foam::cachedRandomMD::vector01()
{
    vector rndVec;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = scalar01();
    }

    return rndVec;
}

// return a normal Gaussian random number
// with zero mean and unity variance N(0, 1)
inline Foam::scalar Foam::cachedRandomMD::GaussNormal()
{
	static bool iset = false;
	static scalar gset;

    if (!iset)
    {
    	volatile long double v1 = 0.0, v2 = 0.0, rsq = 0.0;

    	while(rsq >= 1.0 || rsq == 0.0)
    	{
    		v1 = (2.0*scalar01())-1.0;
			v2 = (2.0*scalar01())-1.0;
			rsq = (v1*v1)+(v2*v2);
    	}

    	//Although log(rsq) should always return a negative number, this is wrapped in an abs() to ensure the value passed to sqrt() is positive
    	scalar absVal = abs(static_cast<scalar>(-2.0*(Foam::log(static_cast<scalar>(rsq))/static_cast<scalar>(rsq))));
		scalar fac = Foam::sqrt(absVal);
        gset = static_cast<scalar>(v1)*fac;
        iset = true;

        return static_cast<scalar>(v2)*fac;
    }
    else
    {
        iset = false;
        return gset;
    }
}

inline void Foam::cachedRandomMD::repopulate()
{
	// Initialise samples based on initial seed
	osRandomSeed(seed_); //Call system seed based on original seed value from application initialisation

	//Create cache equal to the parameterised size but only store the n_th chunk based on the current value of cacheI_
	label sampleSize = samples_.size()+(samples_.size()*cacheI_); //Total sample size to be iterated through
	label sampleStart = samples_.size()*cacheI_; //The start point to begin storing samples
	label sampleIt = 0; //Internal iterator
	scalar rndValue;

	for(label i=0; i<sampleSize; i++)
	{
		//Generate a random value
		rndValue = osRandomDouble();

		//Only start to store once past the threshold
		if(i>=sampleStart){
			samples_[sampleIt] = rndValue;
			sampleIt++; //Increment internal iterator
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cachedRandomMD::cachedRandomMD(const label seed, const label cacheSizeMult)
:
    seed_(1),
    samples_(0),
    sampleI_(-1),
	cacheI_(0),
	cacheSizeMult_(cacheSizeMult)
{
    if (seed > 1)
    {
        seed_ = seed;
    }
}


Foam::cachedRandomMD::cachedRandomMD(const cachedRandomMD& cr, const bool reset)
:
    seed_(cr.seed_),
    samples_(cr.samples_),
    sampleI_(cr.sampleI_),
	cacheI_(cr.cacheI_),
	cacheSizeMult_(cr.cacheSizeMult_)
{
    if (sampleI_ == -1)
    {
        WarningIn
        (
            "Foam::cachedRandomMD::cachedRandomMD(const cachedRandomMD& cr)"
        )   << "Copy constructor called, but samples not being cached. "
            << "This may lead to non-repeatable behaviour" << endl;

        osRandomSeed(seed_);
    }
    else if (reset)
    {
        sampleI_ = 0;
    }
}

Foam::cachedRandomMD::cachedRandomMD
(
    const label seed,
    const label cacheSizeMult,
    const label numMols 
)
:
    seed_(1),
    samples_(0),
    sampleI_(-1),
    cacheI_(0),
    cacheSizeMult_(cacheSizeMult)
{
    if (seed > 1)
    {
        seed_ = seed;
    }
    
    initialise(numMols);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cachedRandomMD::~cachedRandomMD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cachedRandomMD::initialise(const label numMols)
{
	// Samples will be cached if mol count > 0
	if (numMols > 0)
	{
		label cacheSize = numMols * cacheSizeMult_; //Calculate cache size
		samples_.setSize(cacheSize); //Resize sample list
		sampleI_ = 0; //Set sample pointer to zero

		// Initialise samples
		osRandomSeed(seed_); //Call system seed based on original seed value from application initialisation

		//Create cache equal to the parameterised size
		forAll(samples_, i)
		{
			samples_[i] = osRandomDouble();
		}
	}
}

Foam::label Foam::cachedRandomMD::integer(const label lower, const label upper)
{
	if (sampleI_ < 0) //The sample iterator is uninitialised
	{
		return osRandomInteger();
	}

	label rndInt = 0;

    if (sampleI_ == samples_.size() - 1) //The end of the cache has been reached
	{
    	rndInt = static_cast<label>(samples_[sampleI_]);
    	cacheI_++; //Increment the cache counter
		repopulate(); //Call function to repopulate the cache
		sampleI_ = 0;
	}
	else
	{
		rndInt = static_cast<label>(samples_[sampleI_]);
		sampleI_++;
	}

	return lower + (rndInt % (upper+1-lower));
}

template<>
Foam::label Foam::cachedRandomMD::sample01()
{
    return round(scalar01());
}


template<>
Foam::scalar Foam::cachedRandomMD::sample01()
{
    return scalar01();
}


template<>
Foam::label Foam::cachedRandomMD::position(const label& start, const label& end)
{
    return start + round(scalar01()*(end - start));
}


template<>
Foam::scalar Foam::cachedRandomMD::position
(
    const scalar& start,
    const scalar& end
)
{
    return start + scalar01()*(end - start);
}


template<>
Foam::vector Foam::cachedRandomMD::sampleVectorMD()
{
    return vector01();
}


template<>
Foam::scalar Foam::cachedRandomMD::GaussNormalMD()
{
    return GaussNormal();
}


void Foam::cachedRandomMD::operator=(const cachedRandomMD& cr)
{
    seed_ = cr.seed_;
    samples_ = cr.samples_;
    sampleI_ = cr.sampleI_;
    cacheSizeMult_ = cr.cacheSizeMult_;
}


// ************************************************************************* //
