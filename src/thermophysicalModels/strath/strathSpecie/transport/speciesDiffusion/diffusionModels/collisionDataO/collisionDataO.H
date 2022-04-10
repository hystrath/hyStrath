/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
    with correction term for non-nominal electron pressure (i.e., p_e != p_em)

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

        //- Switch to decide whether to use the tabulated collision pair data
        //  (when existing) or the default provisions
        bool isTabulatedData_;
        
        //- Label to decide which limiting electron pressure model to use
        //  0: 1 atm, as in NASA-TM-101528 (Gupta, Yos, Thompson: Feb. 1989)
        //  1: 0.0975*(T/1e3)^4, as in NASA-RP-1232 (Gupta, Yos, Thompson, Lee:
        //                       Aug. 1990)
        label limitingElectronPressureModel_;
        
        //- Local copy of the Boltzmann constant
        const scalar localkB_;
        
        //- Local copy of pi
        const scalar localpi_;
        
        //- Ratio of fundamental constants in the Qc expression
        //  NB: in this ratio, the units of the electric charge ar egiven in esu
        //  and the units of the Boltzmann constant are given in erg/K
        const scalar e2TokB_;

        //- Constant factor in the tabulated collision term expression
        scalar constantFactor_;
        
        //- Model's coefficients
        FixedList<scalar, 4> piOmega1_, piOmega2_;
        
        //- Coefficients cn, Cn, Dn for the computation of ion-ion, ion-electron
        //  and electron-electron collision integrals (when required)
        //  NB: cgs units
        FixedList<scalar, 3> cCDCoeffs1_, cCDCoeffs2_;
        

    // Private member functions

        //- First collision term
        inline scalar collisionTerm1
        (
            const scalar T,
            const scalar pe,
            const scalar Te
        ) const
        {
            scalar Deltaij1 = 0.0;
            
            if (isTabulatedData_)
            {
                Deltaij1 = constantFactor_/sqrt(T)*piOmega11(T);
                
                if (this->collisionType_ > 1)
                {
                    Deltaij1 /= nonNominalPeCorrection(T, pe/(1.01325e5));
                }
            }
            else
            {
                //- Compute the electron number density in cm^3
                const scalar nDe = 1e-6*pe/(localkB_*Te);
                
                Deltaij1 = constantFactor_/sqrt(T)*localpi_
                    *Omega11_charged(T, nDe);
            }
            
            return Deltaij1;
        }

        //- Second collision term
        inline scalar collisionTerm2
        (
            const scalar T,
            const scalar pe,
            const scalar Te
        ) const
        {
            scalar Deltaij2 = 0.0;
            
            if (isTabulatedData_)
            {
                Deltaij2 = 6.0*constantFactor_/(5.0*sqrt(T))*piOmega22(T);
                
                if (this->collisionType_ > 1)
                {
                    Deltaij2 /= nonNominalPeCorrection(T, pe/(1.01325e5));
                }
            }
            else
            {
                //- Compute the electron number density in cm^3
                const scalar nDe = 1e-6*pe/(localkB_*Te);
                
                Deltaij2 = constantFactor_/sqrt(T)*localpi_
                    *Omega22_charged(T, nDe);
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
                
            if (limitingElectronPressureModel_ == 0)
            {
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
            else
            {
                // in this case, the limiting electron pressure is equal to
                // 0.0975*(T/1e3)^4
                return
                    2.0
                  / log
                    (
                        2.09e-02*factor
                      + 1.52*pow(factor, 2.0/3.0)
                    );
            }
        }
        
      //- Additional functions to computate of ion-ion, ion-electron and
      //  electron-electron collision integrals when the input data is missing 
        
        //- Debye shielding length (in cm)
        //  See Eq. 4 in:
        //  M.J. Wright, D. Bose, G.E. Palmer, and E. Levin, "Recommended
        //  Collision Integrals for Transport Property Computations, Part 1: 
        //  Air Species", AIAA Journal, Vol. 43, No. 12, 2005
        inline scalar
        debyeShieldingLength(const scalar T, const scalar nDe) const
        {
            return sqrt(T/(4.0*localpi_*max(nDe,1e-6)*e2TokB_));
        }
        
        //- Dimensionless reduced temperature, T*
        //  See Eq. 5 in:
        //  M.J. Wright, D. Bose, G.E. Palmer, and E. Levin, "Recommended
        //  Collision Integrals for Transport Property Computations, Part 1: 
        //  Air Species", AIAA Journal, Vol. 43, No. 12, 2005
        inline scalar tStar(const scalar T, const scalar nDe) const
        {
            return debyeShieldingLength(T, nDe)*T/e2TokB_;
        }
        
        //- Collision integral Omega^(1,1) (in Angstrom^2)
        //  See Eq. 6 in:
        //  M.J. Wright, D. Bose, G.E. Palmer, and E. Levin, "Recommended
        //  Collision Integrals for Transport Property Computations, Part 1: 
        //  Air Species", AIAA Journal, Vol. 43, No. 12, 2005
        inline scalar Omega11_charged
        (
            const scalar T,
            const scalar nDe
        ) const
        {
            const scalar c1 = cCDCoeffs1_[0];
            const scalar C1 = cCDCoeffs1_[1];
            const scalar D1 = cCDCoeffs1_[2];
            
            const scalar dsl = debyeShieldingLength(T, nDe);
            const scalar tZ = tStar(T, nDe);
            const scalar factor = 5.0e15*sqr(dsl/tZ);
            
            return factor*log( 1.0 + D1*tZ*( 1.0 - C1*exp(-c1*tZ) ) );
        }
        
        //- Collision integral Omega^(2,2) (in Angstrom^2)
        //  See Eq. 6 in:
        //  M.J. Wright, D. Bose, G.E. Palmer, and E. Levin, "Recommended
        //  Collision Integrals for Transport Property Computations, Part 1: 
        //  Air Species", AIAA Journal, Vol. 43, No. 12, 2005
        inline scalar Omega22_charged
        (
            const scalar T,
            const scalar nDe
        ) const
        {
            const scalar c2 = cCDCoeffs2_[0];
            const scalar C2 = cCDCoeffs2_[1];
            const scalar D2 = cCDCoeffs2_[2];
            
            const scalar dsl = debyeShieldingLength(T, nDe);
            const scalar tZ = tStar(T, nDe);
            const scalar factor = 5.0e15*sqr(dsl/tZ);
            
            return factor*log( 1.0 + D2*tZ*( 1.0 - C2*exp(-c2*tZ) ) );
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
            const volScalarField& T,
            const volScalarField& Te
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
            const scalarField& Te,
            const label patchi
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryDiffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
