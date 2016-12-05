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
    Foam::psiThermoHTC

Description
    Compressibility-based thermo model wrapper for hTC models

SourceFiles
    psiThermoHTC.C

\*---------------------------------------------------------------------------*/

#ifndef psiThermoHTC_H
#define psiThermoHTC_H

#include "autoPtr.H"
#include "psiHTCModel.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hTCModels
{

/*---------------------------------------------------------------------------*\
                    class psiThermoHTC Declaration
\*---------------------------------------------------------------------------*/

class psiThermoHTC
:
    public psiHTCModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        psiThermoHTC(const psiThermoHTC&);

        //- Disallow default bitwise assignment
        void operator=(const psiThermoHTC&);


protected:

    // Protected data

        //- Pointer to chemistry model
        autoPtr<psiReactionThermo> thermoPtr_;


public:

    // Constructors

        //- Construct from components and thermo
        psiThermoHTC(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~psiThermoHTC();


    // Member Functions

        //- Return access to the thermo package
        virtual psiReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const psiReactionThermo& thermo() const;

        //- Return const access to the density field
        virtual tmp<volScalarField> rho() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hTCModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
