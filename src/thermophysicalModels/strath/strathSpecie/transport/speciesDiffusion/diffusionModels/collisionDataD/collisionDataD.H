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
    Foam::binaryDiffusionModels::collisionDataD

Description
    Curve-fitted binary diffusion coefficients with correction term for non-
    nominal electron pressure (i.e., p_e != 1 atm)
    _
    D_ij = exp(Dbar_3) * T^( Dbar_0*(lnT)^2 + Dbar_1*lnT + Dbar_2 )

           _
    D_ij = D_ij / p

SourceFiles
    collisionDataD.C

\*---------------------------------------------------------------------------*/

#ifndef collisionDataD_H
#define collisionDataD_H

#include "binaryDiffusionModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binaryDiffusionModels
{

/*---------------------------------------------------------------------------*\
                        Class collisionDataD Declaration
\*---------------------------------------------------------------------------*/

class collisionDataD
:
    public binaryDiffusionModel
{
    // Private data

        //- Model's coefficients
        FixedList<scalar, 4> Dbar_;


    // Private member functions

        //- Model's coefficients given in cm^2*atm/sec
        //  Changed to m^2*Pa/sec in the constructor (see Dbar_[3])
        inline scalar DijBar(const scalar T, const scalar pe) const
        {
            // Binary diffusion coefficient for neutral-neutral interactions
            scalar logT = log(T);
            scalar pDij = Dbar_[3]
              * pow
                (
                    T,
                    Dbar_[0]*sqr(logT) + Dbar_[1]*logT + Dbar_[2]
                );

            // Correction for ion-ion, ion-electron and electron-electron
            // interactions
            if (this->collisionType_ > 1)
            {
                pDij *= nonNominalPeCorrection(T, pe/Pstd);
            }

            return pDij;
        }

        //- Correction factor for non-nominal electron pressure
        //  See Eq. 40.b in:
        //  NASA-TM-101528 (Gupta, Yos, Thompson: Feb. 1989)
        //  Document ID: 19890011822
        //  A review of reaction rates and thermodynamic and transport 
        //  properties for the 11-species air model for chemical and thermal
        //  nonequilibrium calculations to 30000 K
        inline scalar nonNominalPeCorrection
        (
            const scalar T,
            const scalar pe
        ) const
        {
            scalar tprime = 1.0e-12*pow(T, 4.0);
            scalar factor = tprime/(pe+VSMALL);
            
            // the function output is in the range [6e-4, 1]
            return
                log
                (
                    2.09e-02*tprime
                  + 1.52*pow(tprime, 2.0/3.0)
                )
              / log
                (
                    2.09e-02*factor
                  + 1.52*pow(factor, 2.0/3.0)
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
            const scalarField& pe,
            const scalarField& T,
            const label patchi
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryDiffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
