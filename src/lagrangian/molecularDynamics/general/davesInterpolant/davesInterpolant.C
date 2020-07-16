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
    davesInterpolant

Description

\*----------------------------------------------------------------------------*/

#include "davesInterpolant.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

davesInterpolant::davesInterpolant
(
    const IOdictionary& dict
)
:
    propsDict_(dict)
{
    PtrList<entry> infoList(propsDict_.lookup("data"));

    nDataPoints_ = infoList.size();

    density_.setSize(nDataPoints_, 0.0);
    strain_.setSize(nDataPoints_);
    viscosity_.setSize(nDataPoints_);
    slipLength_.setSize(nDataPoints_);

    forAll(infoList, i)
    {
        const entry& setI = infoList[i];
        const dictionary& dictI = setI.dict();

        scalar rhoN = readScalar(dictI.lookup("numberDensity"));

        density_[i] = rhoN;

        List<scalar> strain = List<scalar>(dictI.lookup("strainRate"));
        List<scalar> viscosity = List<scalar>(dictI.lookup("viscosity"));
        List<scalar> slipLength = List<scalar>(dictI.lookup("slipLength"));

        label sizeOfList = strain.size();

        if( sizeOfList != viscosity.size() )
        {
            FatalErrorIn("davesInterpolant::davesInterpolant()")
                << "You are not consistent with the sizes of the list."
                << ", Check the viscosity at rhoN = " << rhoN
                << exit(FatalError);
        }

        if( sizeOfList != slipLength.size() )
        {
            FatalErrorIn("davesInterpolant::davesInterpolant()")
                << "You are not consistent with the sizes of the list."
                << ", Check the slip-length at rhoN = " << rhoN
                << exit(FatalError);
        }

        strain_[i].setSize(sizeOfList, 0.0);
        viscosity_[i].setSize(sizeOfList, 0.0);
        slipLength_[i].setSize(sizeOfList, 0.0);

        forAll(strain, j)
        {
            strain_[i][j] = strain[j];
            viscosity_[i][j] = viscosity[j];
            slipLength_[i][j] = slipLength[j];
        }
    }

    Info << "densities = " <<  density_ << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

davesInterpolant::~davesInterpolant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
