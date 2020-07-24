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
    Foam::binaryDiffusionModels::collisionDataO

Description
    Binary diffusion coefficients obtained from curve-fitted collision integrals
    with correction term for non-nominal electron pressure (i.e., p_e != 1 atm)

SourceFiles
    collisionDataO.C

\*---------------------------------------------------------------------------*/

#ifndef collisionDataO_H
#define collisionDataO_H

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
                       Class collisionDataO Declaration
\*---------------------------------------------------------------------------*/

class collisionDataO
:
    public binaryDiffusionModel
{
    // Private data

        //- Local copy of the Boltzmann constant
        const scalar localkB_;
        
        //- Ratio of fundamental constants in the Qc expression
        const scalar e4OverkB2_;

        //- Constant factor in the tabulated collision term expression
        scalar constantFactor_;
        
        //- Model's coefficients
        FixedList<scalar, 4> piOmega1_, piOmega2_;


    // Private member functions

        //- First collision term
        inline scalar collisionTerm1(const scalar T, const scalar pe) const
        {
            scalar Deltaij1 = constantFactor_/sqrt(T)*piOmega11(T);
            
            if (this->collisionType_ > 1)
            {
                Deltaij1 *= nonNominalPeCorrection(T, pe/Pstd);
            }

            return Deltaij1;
        }

        //- Second collision term
        inline scalar collisionTerm2(const scalar T, const scalar pe) const
        {
            scalar Deltaij2 = 6.0*constantFactor_/(5.0*sqrt(T))*piOmega22(T);
            
            if (this->collisionType_ > 1)
            {
                Deltaij2 *= nonNominalPeCorrection(T, pe/Pstd);
            }

            return Deltaij2;
        }

        //- Collision integrals pi*Omega_ij^(1,1) for tabulated colliding pairs
        //  in Angstrom^2
        inline scalar piOmega11(const scalar T) const
        {
            scalar logT = log(T);
            return piOmega1_[3]
              * pow
                (
                    T,
                    piOmega1_[0]*sqr(logT) + piOmega1_[1]*logT + piOmega1_[2]
                );
        }

        //- Collision integrals pi*Omega_ij^(2,2) for tabulated colliding pairs
        //  in Angstrom^2
        inline scalar piOmega22(const scalar T) const
        {
            scalar logT = log(T);
            return piOmega2_[3]
              * pow
                (
                    T,
                    piOmega2_[0]*sqr(logT) + piOmega2_[1]*logT + piOmega2_[2]
                );
        }

        //- Collision integrals pi*Omega_ij^(1,1) for non-tabulated colliding
        //  pairs: charged / charged interactions (in Angstrom^2)
        inline
        scalar piOmega11NonNeutral(const scalar T, const scalar pe) const
        {
            // It is assumed that Z = 1, i.e., that the species is
            // singly-ionised (see Eq. 23.c in Ref Gupta 1989)
/*            if (this->name1_ == "e-" and this->name2_ == "e-") {*/
/*            for (scalar tt = 500; tt <=30000; tt+=500 )*/
/*            {*/
/*                Info<< tt << tab 
                      << 0.795*e4OverkB2_*log(shieldingParameter(tt, 1))/sqr(tt)
                      << endl;*/
/*            }*/
/*            FatalErrorIn*/
/*            (*/
/*                "DiffusionModel - my error"*/
/*            )   << tab << exit(FatalError);}*/
            
            return 0.795*e4OverkB2_/sqr(T)*log(shieldingParameter(T, pe));
        }

        //- Collision integrals pi*Omega_ij^(2,2) for non-tabulated colliding
        //  pairs: e- / e- interaction (in Angstrom^2)
        inline
        scalar piOmega22NonNeutralE(const scalar T, const scalar pe) const
        {
            return 1.29*e4OverkB2_/sqr(T)*log(shieldingParameter(T, pe));
        }

        //- Collision integrals pi*Omega_ij^(2,2) for non-tabulated colliding
        //  pairs: ion / ion interactions (in Angstrom^2)
        inline
        scalar piOmega22NonNeutralI(const scalar T, const scalar pe) const
        {
            // It is assumed that Z = 1, i.e., that the species is
            // singly-ionised (see Eq. 23.c in Ref Gupta 1989)
            return 1.36*e4OverkB2_/sqr(T)*log(shieldingParameter(T, pe));
        }

        //- Shielding parameter
        inline scalar shieldingParameter(const scalar T, const scalar pe) const
        {
            scalar factor = 1.0e-12*pow(T, 4.0)/pe;
            return sqrt(2.09e-02*factor + 1.52*pow(factor, 2.0/3.0));
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
    TypeName("collisionDataO");


    // Constructors

        //- construct from components
        collisionDataO
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

        ~collisionDataO()
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
