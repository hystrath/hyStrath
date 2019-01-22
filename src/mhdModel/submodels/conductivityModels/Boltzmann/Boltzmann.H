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
    Foam::mhd::Boltzmann

Description
    Boltzmann conductivity model for MHD terms

SourceFiles
    SpitzerHarm.C

\*---------------------------------------------------------------------------*/

#ifndef mhdBoltzmann_H
#define mhdBoltzmann_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                  Class constantAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class Boltzmann
:
    public conductivityModel
{
    // Private data
    word crossSections_;

public:

    //- Runtime type information
    TypeName("Boltzmann");


    // Constructors

        //- Construct from components
        Boltzmann(const mhdModel& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~Boltzmann();


    // Member Functions

        // Access


                //- conductivity field
                volScalarField sigma() const;

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
