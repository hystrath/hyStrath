/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::rho2ChemistryHTC

Description
    Density-based chemistry model wrapper for hTC2 models

SourceFiles
    rho2ChemistryHTC.C

\*---------------------------------------------------------------------------*/

#ifndef rho2ChemistryHTC_H
#define rho2ChemistryHTC_H

#include "autoPtr.H"
#include "rho2HTCModel.H"
#include "rho2ChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hTC2Models
{

/*---------------------------------------------------------------------------*\
                    class rho2ChemistryHTC Declaration
\*---------------------------------------------------------------------------*/

class rho2ChemistryHTC
:
    public rho2HTCModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rho2ChemistryHTC(const rho2ChemistryHTC&);

        //- Disallow default bitwise assignment
        void operator=(const rho2ChemistryHTC&);


protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<rho2ChemistryModel> chemistryPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        rho2ChemistryHTC(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~rho2ChemistryHTC();


    // Member Functions

        //- Return access to the thermo package
        virtual rho2ReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const rho2ReactionThermo& thermo() const;

        //- Return const access to the density field
        virtual tmp<volScalarField> rho() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hTC2Models
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
