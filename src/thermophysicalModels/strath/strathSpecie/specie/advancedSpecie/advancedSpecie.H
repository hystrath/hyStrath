/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::advancedSpecie

Description
    Base class of the thermophysical property types.

SourceFiles
    advancedSpecieI.H
    advancedSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef advancedSpecie_H
#define advancedSpecie_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"

#include "DynamicList.H" // NEW VINCENT 24/02/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class advancedSpecie Declaration
\*---------------------------------------------------------------------------*/

class advancedSpecie
{
    // Private data
      
        //- Name of specie
        word name_;

        //- Number of moles of this component in the mixture
        scalar nMoles_;

        //- Molecular weight of specie [kg/kmol]
        scalar molWeight_;
        
        //- Type of particle
        label particleType_;
        
        //- Charge of the particle
        label particleCharge_;
        
        //- Diameter of the particle
        scalar diameter_;
        
        //- Temperature-viscosity coefficient
        scalar omega_;
        
        //- Vibrational list with degeneracy and characteristic vibrational temperature [K] for each vib. mode
        DynamicList<scalar> vibrationalList_;
        
        //- Dissociation potential of the particle [J/kg]
        scalar dissociationPotential_;
        
        //- Number of characteristic vibrational temperatures
        label noVibrationalTemp_;
        
        //- Number of electronic levels
        label noElectronicLevels_;
        
        //- Fraction of the first ionisation energy of a given species [J/kg]
        scalar iHat_;


public:

    //- Runtime type information
    ClassName("advancedSpecie");


    // Public constants

        // Thermodynamic constants

            //- Universal gas constant [J/(kmol K)]
            static const scalar RR;

            //- Standard pressure [Pa]
            static const scalar Pstd;

            //- Standard temperature [K]
            static const scalar Tstd;


    // Constructors


        //- Construct from components without name
        inline advancedSpecie
        (
            const scalar nMoles, 
            const scalar molWeight,
            const label particleType,
            const label particleCharge,
            const scalar diameter,
            const scalar omega,
            const DynamicList<scalar> vibrationalList,
            const scalar dissociationPotential,
            const label noVibrationalTemp,
            const label noElectronicLevels,
            const scalar iHat
        );

        //- Construct from components with name
        inline advancedSpecie
        (
            const word& name,
            const scalar nMoles,
            const scalar molWeight,
            const label particleType,
            const label particleCharge,
            const scalar diameter,
            const scalar omega,
            const DynamicList<scalar> vibrationalList,
            const scalar dissociationPotential,
            const label noVibrationalTemp,
            const label noElectronicLevels,
            const scalar iHat
        );

        //- Construct as copy
        inline advancedSpecie(const advancedSpecie&);

        //- Construct as named copy
        inline advancedSpecie(const word& name, const advancedSpecie&);

        //- Construct from Istream
        advancedSpecie(Istream&);

        //- Construct from dictionary
        advancedSpecie(const dictionary& dict);


    // Member Functions

        // Access
          
            //- Name
            inline const word& nameSpecies() const;

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- No of moles of this species in mixture
            inline scalar nMoles() const;

            //- Gas constant [J/(kg K)]
            inline scalar R() const;
            
            //- Return const access to the type of particle
            inline label particleType() const;
            
            //- Return const access to the charge of the particle
            inline label particleCharge() const;
            
            //- Return const access to the diameter of the particle
            inline scalar diameter() const;
            
            //- Return const access to the temperature-viscosity coefficient
            inline scalar omega() const;
            
            //- Return const access to the vibrational list
            inline DynamicList<scalar> vibrationalList() const;
            
            //- Return const access to the dissociation potential of the particle [K]
            inline scalar dissociationPotential() const;
            
            //- Return const access to the number of characteristic vibrational temperatures
            inline label noVibrationalTemp() const;
            
            //- Return const access to the number of electronic energy levels
            inline label noElectronicLevels() const;
            
            //- Return const access to the fraction of the first ionisation energy
            inline scalar iHat() const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const advancedSpecie&);

        inline void operator+=(const advancedSpecie&);
        inline void operator-=(const advancedSpecie&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend advancedSpecie operator+(const advancedSpecie&, const advancedSpecie&);
        inline friend advancedSpecie operator-(const advancedSpecie&, const advancedSpecie&);

        inline friend advancedSpecie operator*(const scalar, const advancedSpecie&);

        inline friend advancedSpecie operator==(const advancedSpecie&, const advancedSpecie&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const advancedSpecie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "advancedSpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
