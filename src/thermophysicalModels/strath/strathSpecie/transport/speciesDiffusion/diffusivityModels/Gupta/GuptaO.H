/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::binaryDiffusivityModels::GuptaO

Description
    GuptaO binary diffusivity model.
    
    See:
    R. N. Gupta, J. M. Yos, R. A. Thompson, and A. Richard 1989
    A Review of Reaction Rates and Thermodynamic and Transport Properties
    for the 11-Species Air Model for Chemical and Thermal Nonequilbrium 
    Calculations to 30000 K. 
    NASA TM-101528
    
    and
    
    R. N. Gupta, J. M. Yos, R. A. Thompson, and K.-P. Lee. 1990
    A Review of Reaction Rates and Thermodynamic and Transport Properties for
    an 11-Species Air Model for Chemical and Thermal Nonequilbrium Calculations
    to 30000 K.
    NASA RP-1232
    
SourceFiles
    GuptaO.C

\*---------------------------------------------------------------------------*/

#ifndef GuptaO_H
#define GuptaO_H

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
                           Class GuptaO Declaration
\*---------------------------------------------------------------------------*/

class GuptaO
:
    public binaryDiffusivityModel
{
    // Private data

        //- Model coefficients
        FixedList<scalar, 4> piOmega1_;
        
        //- Molecular weights
        const scalar W1_, W2_;
        
        //- Fundamental constants (more efficient)
        const scalar pi, kB, Runi;
        
        //- Constant in the neutral collision term
        const scalar constantFactorInCollisionTerm_;
        
        // Ratio of fundamental constants in the expression of Qc
        const scalar e4OverkB2_;

        
    // Private member functions
        
        //- Collision term 1
        inline scalar collisionTerm1(const scalar T, const scalar pe = 0) const
        {
            if(T != 0.0)
            {
                if(this->collisionType_ < 2)
                {
                    return constantFactorInCollisionTerm_/(3.0*sqrt(T))
                        *piOmega11Neutral(T);
                    
                }
                else if(pe > SMALL)
                {
                    return constantFactorInCollisionTerm_
                        *piOmega11NonNeutral(T, pe);    
                }
            }
            
            return VSMALL;
        }
        
        //- Collision term 2 (unused)
        /*inline scalar collisionTerm2Neutral(const scalar T) const 
        {
            if(T != 0.0)
            {
                if(this->collisionType_ < 2)
                {
                    return 2.0*constantFactorInCollisionTerms_/(5.0*sqrt(T))
                        *piOmega22Neutral(T); // to implement should it be needed, would require piOmega2_
                    
                }
                else if(pe > SMALL)
                {
                    if(this->collisionType_ != 4)
                    {
                        return constantFactorInCollisionTerm_
                            *piOmega22NonNeutralI(T, pe);    
                    }
                    else
                    {
                        return constantFactorInCollisionTerm_
                            *piOmega22NonNeutralE(T, pe);
                    }
                }
            }

            return VSMALL;
        }*/
        
        //- Collision cross-section (1,1) for neutral interactions
        //  in Angstrom^2
        inline scalar piOmega11Neutral(const scalar T) const 
        {
            return piOmega1_[3]
                *pow
                (
                    T, 
                    piOmega1_[0]*sqr(log(T)) 
                  + piOmega1_[1]*log(T) + piOmega1_[2]
                );
        } 
        
        //- Collision cross-section (2,2) for neutral interactions
        //  in Angstrom^2 (unused)
        /*inline scalar piOmega22Neutral(const scalar T) const 
        {
            return piOmega2_[3]
                *pow
                (
                    T, 
                    piOmega2_[0]*sqr(log(T)) 
                  + piOmega2_[1]*log(T) + piOmega2_[2]
                );
        }*/
        
        //- Collision cross-section (1,1) for non-neutral interactions
        //  in Angstrom^2
        inline scalar piOmega11NonNeutral(const scalar T, const scalar pe) const 
        {
            // It is assuned that Z = 1, name that the species is 
            // singly-ionised (see Eq. 23.c in Ref Gupta 1989)
            return 0.795*e4OverkB2_*log(shieldingParameter(T, pe))/sqr(T);
        }
        
        //- Collision cross-section (2,2) for e- / e- interactions
        //  in Angstrom^2 (unused)
        inline scalar piOmega22NonNeutralE(const scalar T, const scalar pe) const 
        {
            return 1.29*e4OverkB2_*log(shieldingParameter(T, pe))/sqr(T);
        }
        
        //- Collision cross-section (2,2) for ion / ion interactions
        //  in Angstrom^2 (unused)
        inline scalar piOmega22NonNeutralI(const scalar T, const scalar pe) const 
        {
            // It is assuned that Z = 1, name that the species is 
            // singly-ionised (see Eq. 23.c in Ref Gupta 1989)
            return 1.36*e4OverkB2_*log(shieldingParameter(T, pe))/sqr(T);
        }
        
        //- Shielding parameter function of the temperature and electron 
        //  pressure.
        inline scalar shieldingParameter(const scalar T, const scalar pe) const
        {
            return sqrt(2.09e-14*pow(T, 4.0)/pe 
                + 1.52*pow(pow(T, 4.0)/pe, 2.0/3.0));
        }   


public:

    //- Runtime type information
    TypeName("GuptaO");


    // Constructors

        //- construct from components
        GuptaO
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

        ~GuptaO()
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
