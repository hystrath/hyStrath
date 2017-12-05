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

#include "maitlandSmith.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(maitlandSmith, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    maitlandSmith,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maitlandSmith::maitlandSmith
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& maitlandSmith
)
:
    pairPotential(name, rU, maitlandSmith),
    maitlandSmithCoeffs_(maitlandSmith.subDict(typeName + "Coeffs")),
    m_(readScalar(maitlandSmithCoeffs_.lookup("m"))),
    gamma_(readScalar(maitlandSmithCoeffs_.lookup("gamma"))),
    rm_(readScalar(maitlandSmithCoeffs_.lookup("rm"))),
    epsilon_(readScalar(maitlandSmithCoeffs_.lookup("epsilon")))
{
    // temporary error
    FatalErrorIn("maitlandSmith::maitlandSmith()")
        << "Check and modify the code for the maitlandSmith model (i.e. you need to make sure that the coefficents are changed to reduced units if you are using reduced units)"
        << nl << "in: " << "potentialDict"
        << exit(FatalError);

    // see below too
    if(rU.runReducedUnits())
    {
        rm_ /= rU.refLength();
        epsilon_ /= rU.refEnergy();
    }

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar maitlandSmith::unscaledEnergy(const scalar r) const
{
    scalar nr = (m_ + gamma_*(r/rm_ - 1.0));

    return epsilon_
       *(
            (6.0 / (nr - 6.0))*Foam::pow(r/rm_, -nr)
          - (nr / (nr - 6.0))*Foam::pow(r/rm_, -6)
        );
}


bool maitlandSmith::read(const dictionary& pairPotentialProperties, const reducedUnits& rU)
{
    pairPotential::read(pairPotentialProperties, rU);

    maitlandSmithCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    maitlandSmithCoeffs_.lookup("m") >> m_;
    maitlandSmithCoeffs_.lookup("gamma") >> gamma_;
    maitlandSmithCoeffs_.lookup("rm") >> rm_;
    maitlandSmithCoeffs_.lookup("epsilon") >> epsilon_;

    if(rU.runReducedUnits())
    {
        rm_ /= rU.refLength();
        epsilon_ /= rU.refEnergy();
    }

    return true;
}

const dictionary& maitlandSmith::dict() const
{
    return pairPotentialProperties_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
