/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::heRho2Thermo

Description
    Energy for a mixture based on density

SourceFiles
    heRho2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRho2Thermo_H
#define heRho2Thermo_H

#include "rho2Thermo.H"
#include "he2Thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRho2Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsi2Thermo, class MixtureType>
class heRho2Thermo
:
    public he2Thermo<BasicPsi2Thermo, MixtureType>
{
    // Private Member Functions

      //- Initialiase the thermo variables
      void init2();
        
      //- Calculate the thermo variables (2T activated)
      void calculate2();

      //- Construct as copy (not implemented)
      heRho2Thermo(const heRho2Thermo<BasicPsi2Thermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("heRho2Thermo");


    // Constructors

        //- Construct from mesh and phase name
        heRho2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~heRho2Thermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "heRho2Thermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
