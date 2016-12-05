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
    Foam::VTRelaxationModels::SchwartzSlawskyHerzfeld

Description
    SchwartzSlawskyHerzfeld (SSH) V-T Relaxation model.
    
SourceFiles
    SchwartzSlawskyHerzfeld.C

\*---------------------------------------------------------------------------*/

#ifndef SchwartzSlawskyHerzfeld_H
#define SchwartzSlawskyHerzfeld_H

#include "VTRelaxationModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "fundamentalConstants.H"
#include "physicoChemicalConstants.H"

#include "FixedList.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace VTRelaxationModels
{

/*---------------------------------------------------------------------------*\
                           Class SchwartzSlawskyHerzfeld Declaration
\*---------------------------------------------------------------------------*/

class SchwartzSlawskyHerzfeld
:
    public VTRelaxationModel
{
    // Private data

        // Are species 1 and 2 one and only species?
        bool matchingSpeciesIndices_;
        
        // Species positions in the species() list
        label species1_, species2_;
        
        // SSH OR MWP VT energy exchange?
        bool SHHon_;
        
        // Temperatures for the linear interpolation
        scalar Tlow_, Thigh_;
        
        // Weight of a single particle, species 1 [kg]
        scalar m1_; 
        
        // Reduced mass [kg]
        scalar mu12_; 
        
        // Characteristic vibrational temperature, species 1 [K]
        scalar TH1_; 
        
        // Interatomic distnce of the molecule 1, at the ground vibrational level
        scalar rp1_; 
        
        // factor in ZVs
        scalar fms_; 
        
        // Number density of the mixture
        scalarField nD12_;   
        
        // Lennard and Jones parameters
        FixedList<scalar,2> sigma12_, epsilon12_; 
        
        
        //- Interaction with atoms: switch to MWP
        
        // Coefficients
        scalar A12_, B12_; 
        
        // Molecular weight, species i [g/mol]
        scalar W1_;   
        
        // Offset term in the exponential
        scalar offset_;   
        
        // Coefficients related to Park's correction
        scalar sigma1_, sigma2_; 
        
        // Fundamental constants
        scalar pi, kB, hPlanck;


public:

    //- Runtime type information
    TypeName("SSH"); //- SSH for Schwartz, Slawsky, and Herzfeld


    // Constructors

        //- construct from components
        SchwartzSlawskyHerzfeld
        (
            const word& name1,
            const word& name2,
            const label& lname1,
            const label& lname2,
            const dictionary& dict2T,
            const dictionary& dictThermoPhy,
            const volScalarField& p,
            const volScalarField& Tt,
            const PtrList<volScalarField>& Tv,
            const PtrList<volScalarField>& nD
        );


    // Destructor

        ~SchwartzSlawskyHerzfeld()
        {}


    // Member Functions

        //- V-T relaxation time [s]
        virtual tmp<volScalarField> tauVT() const;
         
        //- V-T relaxation time for patch [s]
        virtual tmp<scalarField> tauVT
        (
            const label patchi,
            const scalarField& p,
            const scalarField& Tt,
            const PtrList<scalarField>& Tv,
            const PtrList<scalarField>& nD 
        ) const;
        
        scalar linearInterpolation //- bounded
        (
            const scalar xpt,
            const scalar x1, 
            const scalar x2,
            const scalar y1, 
            const scalar y2
        ) const
        {
            /*if (xpt<x1)
            {
                return y1;
            }
            else if (xpt>x2)
            {
                return y2;
            }
            else
            {*/
                return max(min(y1 + (xpt-x1)*(y2-y1)/(x2-x1), y2), y1);
            //}
        }
        
        scalar sigmaLJmixture 
        (
            const word name1,
            const scalar molW1,
            const scalar sigma1,
            const scalar epsilon1,
            const word name2,
            const scalar sigma2,
            const scalar epsilon2,
            const scalar polarizability2
        ) const
        {
            scalar ksi = 0.0;

            if(similarPolarization(name1, name2))
            {
                return 0.5*(sigma1+sigma2);
            }
            else
            {
                ksi = 1.0; //+ 0.25*polarizability2*sqr(molW1/2.0/sqrt(epsilon1*pow3(sigma1)))*sqrt(epsilon1/epsilon2)/pow3(sigma2);
                return 0.5*(sigma1+sigma2)*pow(ksi, -1.0/6.0);
            }
        }
        
        scalar epsilonLJmixture 
        (
            const scalar molW1,
            const scalar sigma1,
            const scalar epsilon1,
            const scalar sigma2,
            const scalar epsilon2,
            const scalar polarizability2
        ) const
        {
            scalar ksi = 1.0; //+ 0.25*polarizability2*sqr(molW1/2.0/sqrt(epsilon1*pow3(sigma1)))*sqrt(epsilon1/epsilon2)/pow3(sigma2);
            return sqrt(epsilon1*epsilon2)*sqr(ksi);
        }
        
        
        bool similarPolarization
        (
            const word name1,
            const word name2
        ) const
        {
            bool pol1(false), pol2(false);
            if (name1.find("+") != std::string::npos or name1.find("-") != std::string::npos) 
            {
                pol1 = true;
            }
            if (name2.find("+") != std::string::npos or name2.find("-") != std::string::npos) 
            {
                pol2 = true;
            }
            return (not(pol1 xor pol2));
        } 
        
        
        inline bool similarSpecies
        (
            const label name1,
            const label name2
        ) const
        {
            if (name1 == name2) 
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
        
        template <class T>
        T nDmixture 
        (
            const T nD1,
            const T nD2
        ) const
        {
            if(matchingSpeciesIndices_)
            {
                return nD1;
            }
            else
            {
                return nD1+nD2;
            }
        }
        
        inline scalar deltaE() const
        {
            return kB*TH1_;
        }
        
        inline scalar beta
        (
            const scalar Tt
        ) const
        {
            scalar coeff1 = pow(2.0*epsilon12interpolation(Tt)/mu12_, 9.0);
            scalar coeff2 = pow6((3.0*hPlanck*mu12_)/(sqr(pi)*sigma12interpolation(Tt)
                *kB*Tt*deltaE()));
            return pow(0.5*coeff1*coeff2, 1.0/19.0);
        }
        
        inline scalar deltaStarmixture
        (
            const scalar Tt
        ) const
        {
            return 1.0/sigma12interpolation(Tt)*(12.0/pow(2.0*beta(Tt), 1.0/6.0))*(1.0+21.0/19.0*beta(Tt));
        }
        
        inline scalar alphasr
        (
            const scalar Tt
        ) const
        {
            return (16.0*pow4(pi)*mu12_*sqr(deltaE()))
                /(sqr(hPlanck*deltaStarmixture(Tt))*kB);
        }
        
        inline scalar chiStarsr
        (
            const scalar Tt
        ) const
        {
            return 0.5*pow(alphasr(Tt)/Tt, 1.0/3.0);
        }
        
        inline scalar rcStarmixture
        (
            const scalar Tt
        ) const
        {
            return sigma12interpolation(Tt)/(pow(2.0*beta(Tt), 1.0/6.0)*(1.0+2.0/19.0*beta(Tt)));
        }
        
        inline scalar collisionalArea
        (
            const scalar Tt
        ) const
        {
            return sqr(rcStarmixture(Tt)/sigma12interpolation(Tt));
        }
        
        inline scalar Z0s
        (
            const scalar Tt
        ) const
        {
            return deltaStarmixture(Tt)*rp1_ + 5.0/2.0/sqr(deltaStarmixture(Tt)*rp1_);
        }
        
        inline scalar ZPlussr
        (
            const scalar Tt
        ) const
        {
            scalar inter = Tt*kB;
            return exp(-4.0/pi*sqrt(epsilon12interpolation(Tt)*chiStarsr(Tt)/inter) -16.0/(3.0*sqr(pi))*epsilon12interpolation(Tt)/inter);
        }
        
        
        inline scalar ZTsr
        (
            const scalar Tt
        ) const
        {
            scalar inter1 = exp(1.5*pow(alphasr(Tt) /Tt, 1.0/3.0) - deltaE()/(2.0*Tt*kB));
            scalar inter2 = sqr(pi)*sqrt(1.5/pi)*sqr(deltaE()/(kB*alphasr(Tt))) * pow(Tt/alphasr(Tt), 1.0/6.0);
            return inter1*inter2;
        }
        
        inline scalar ZVs
        (
            const label i,
            const scalar Tt
        ) const
       {
            scalar inter1 = 0;
            for(label vibState=0; vibState < 35; vibState++)
            {
                inter1 += 1.0/((vibState+1)*sqr(pi)); //*alphasr(Tt)/TH1_*sqr(kB*TH1_/deltaE());
            }
            //inter1 += 0.5/sqr(pi)*alphasr(Tt)/TH1_*sqr(kB*TH1_/deltaE());
            return inter1*reducedMass(m1_,m1_)/mu12_*alphasr(Tt)/TH1_*sqr(kB*TH1_/deltaE());
        }
        
        inline scalar fms
        (
            const word name
        ) const
        {
            if (name.find("2") != std::string::npos) 
            {
                return 1.0;
            }
            else
            {
                // I am assuming the species is NO
                scalar mN = 23.25e-27;
                scalar mO = 26.56e-27;
                scalar mNO = mN + mO;
                scalar Cpbarsquare = (sqr(mN)+sqr(mO))/(2.0*mN*mO*(mN+mO));
                return 1.0/(mNO*Cpbarsquare);
            }
        }
        
        inline scalar P10sr
        (
            const scalar Tt
        ) const
        {
            /*Info << "collisionalArea(Tt): " << collisionalArea(Tt)<< endl;
            Info << "Z0s(Tt): " << Z0s(Tt)<< endl;
            Info << "ZVs(1,Tt): " << ZVs(1,Tt)<< endl;
            Info << "ZTsr(Tt): " << ZTsr(Tt)<< endl;
            Info << "ZPlussr(Tt): " << ZPlussr(Tt)<< endl;*/
            return collisionalArea(Tt)/(Z0s(Tt)*ZVs(1,Tt)*ZTsr(Tt)*ZPlussr(Tt));
        }
        
        inline scalar Zcollsr
        (
            const scalar Tt,
            const scalar nD1,
            const scalar nD2
        ) const
        {
            return 2.0*sqr(sigma12interpolation(Tt))*max(nDmixture<scalar>(nD1, nD2), Foam::SMALL)
              * sqrt(2.0*pi*kB*Tt/mu12_);
        }
        
        inline scalar sigma12interpolation
        (
            const scalar Tt
        ) const
        {
            return linearInterpolation(Tt, Tlow_, Thigh_, sigma12_[0], sigma12_[1]);
        }
        
        inline scalar epsilon12interpolation
        (
            const scalar Tt
        ) const
        {
            return linearInterpolation(Tt, Tlow_, Thigh_, epsilon12_[0], epsilon12_[1]);
        }
        
        inline scalar reducedMass
        (
            const scalar m1,
            const scalar m2
        ) const
        {
            return m1*m2/(m1+m2);
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VTRelaxationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
