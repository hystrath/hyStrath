/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::NoElectron

Description
    Boltzman electron fluid model

\*---------------------------------------------------------------------------*/

#ifndef Boltz_H
#define Boltz_H

#include "ElectronModel.H"
#include "pdCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class NoElectron Declaration
\*---------------------------------------------------------------------------*/

class Boltz
:
    public ElectronModel
{

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Electron reference density
        const dimensionedScalar n0_;

        //- Electron reference temperature (eV)
        dimensionedScalar TeV_;

        //- Electron reference temperature (K)
        const dimensionedScalar Te_;

        //- Reference potential [V]
        const dimensionedScalar phi0_;

        //- Non-linear solver convergence criteria
        const scalar convTol_;

        //- Number of iterations done by non-linear solver
        const scalar nCorr_;

        volScalarField phiN_;

        volScalarField kappa_;

        volScalarField dphi_;



public:

    //- Runtime type information
    TypeName("Boltz");


    // Constructors

        //- Construct from dictionary
        Boltz
        (
            const dictionary& dict,
            pdCloud& cloud
        );


    //- Destructor
        ~Boltz();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        bool active() const;

        //Calculate model's contribution to potential
        void calculatePotential();

        //Calculate model's electron density
        void calculateRhoe();

        const dictionary& coeffDict() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "Boltz.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
