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
    Foam::binaryDiffusivityModels::collisionDataD

Description
    collisionDataD binary diffusivity model

SourceFiles
    collisionDataD.C

\*---------------------------------------------------------------------------*/

#ifndef collisionDataD_H
#define collisionDataD_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binaryDiffusivityModels
{

/*---------------------------------------------------------------------------*\
                       Class collisionDataD Declaration
\*---------------------------------------------------------------------------*/

class collisionDataD
:
    public binaryDiffusivityModel
{
    // Private data

        //- Model's coefficients
        FixedList<scalar, 4> Dbar_;


    // Private member functions

        //- Model's coefficients given in cm2*atm/sec
        //  Changed to m^2*Pa/sec in the constructor
        inline scalar DijBar(const scalar T, const scalar pe = 0) const
        {
            scalar Dij = Dbar_[3]
                *pow
                (
                    T,
                    Dbar_[0]*sqr(log(T))
                  + Dbar_[1]*log(T) + Dbar_[2]
                );

            // ion-ion, ion-electron or electron-electron
            if(this->collisionType_ > 1 and pe > SMALL)
            {
               Dij *= nonNominalpeCorrection(T, pe/1.01325e5);
            }

            return Dij;
        }

        //- Correction factor for non-nominal electron pressure
        //  (see Eq. 40.b [1989])
        inline scalar nonNominalpeCorrection
        (
            const scalar T,
            const scalar pe
        ) const
        {
            return log
                   (
                      2.09e-14*pow(T, 4.0)
                    + 1.52*pow(T/1000.0, 8.0/3.0)
                   )
                   /log
                   (
                      2.09e-14*pow(T/pe, 4.0)
                    + 1.52*pow(T/(1000.0*pe), 8.0/3.0)
                   );
        }


public:

    //- Runtime type information
    TypeName("collisionDataD");


    // Constructors

        //- construct from components
        collisionDataD
        (
            const word& name1,
            const word& name2,
            const dictionary& dictThermo,
            const dictionary& dictTransport,
            const volScalarField& p,
            const volScalarField& pe,
            const volScalarField& T
        );


    // Destructor

        ~collisionDataD()
        {}


    // Member Functions

        //- Binary diffusion coefficient [m^2/s]
        virtual tmp<volScalarField> D() const;

        //- Binary diffusion coefficient for patch [m^2/s]
        virtual tmp<scalarField> D
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;

        //- Binary diffusion coefficient for patch [m^2/s]
        virtual tmp<scalarField> D
        (
            const scalarField& p,
            const scalarField& pe,
            const scalarField& T,
            const label patchi
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryDiffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
