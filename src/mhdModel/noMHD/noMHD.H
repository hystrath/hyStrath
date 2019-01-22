/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::radiation::noRadiation

Description
    No MHD - does nothing to energy and momentum equation source terms
    (returns zeros)

SourceFiles
    noMHD.C

\*---------------------------------------------------------------------------*/

#ifndef noMHD_H
#define noMHD_H

#include "mhdModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                         Class noMHD Declaration
\*---------------------------------------------------------------------------*/

class noMHD
:
    public mhdModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noMHD(const noMHD&);

        //- Disallow default bitwise assignment
        void operator=(const noMHD&);


public:

    //- Runtime type information
    TypeName("noMHD");


    // Constructors

        //- Construct from components
        noMHD(const rho2ReactionThermo& thermo);

        //- Construct from components
        noMHD(const dictionary& dict, const rho2ReactionThermo& thermo);


    //- Destructor
    virtual ~noMHD();


    // Member functions

        // Edit

            //- Read radiationProperties dictionary
            bool read();
       //Access
     //- Source term component for energy equation
     tmp<volScalarField> Q(const volVectorField& U) const;
     // - Source term component for momentum equation
     tmp<volVectorField> F(const volVectorField& U) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
