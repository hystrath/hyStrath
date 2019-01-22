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
    Standart mhd model for flows with low magnetic Reynolds number

SourceFiles
    lowReMag.C

\*---------------------------------------------------------------------------*/

#ifndef lowReMag_H
#define lowReMag_H

#include "mhdModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                         Class noMHD Declaration
\*---------------------------------------------------------------------------*/

class lowReMag
:
    public mhdModel
{
    // Private Member Functions
        //magnetic field
        volVectorField B_;
        //conductivity value
        volScalarField sigma_;

        //- Disallow default bitwise copy construct
        lowReMag(const lowReMag&);

        //- Disallow default bitwise assignment
        void operator=(const lowReMag&);


public:

    //- Runtime type information
    TypeName("lowReMag");


    // Constructors

        //- Construct from components
        lowReMag(const rho2ReactionThermo& thermo);

        //- Construct from components
        lowReMag(const dictionary& dict, const rho2ReactionThermo& thermo);


    //- Destructor
    virtual ~lowReMag();


    // Member functions

        // Edit

            //- Read radiationProperties dictionary
            bool read();
       //Access
     void calculateSigma();
     // - Current density
     tmp<volVectorField> j(const volVectorField& U) const;
     //- Source term component for energy equation
     tmp<volScalarField> Q(const volVectorField& U) const;
     // - Source term component for momentum equation
     tmp<volVectorField> F(const volVectorField& U) const;
     // - Source term component for momentum equation
     tmp<volScalarField> Stuart(const volVectorField& U) const;
     //- Calculate Hall Effect tensor
     tmp<volTensorField> hallCorrection() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
