/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::ImpactIonisationArrhenius2ReactionRate

Description
    Arrhenius reaction rate given by:

        k = A * T^beta * exp(-Ta/T)

SourceFiles
    ImpactIonisationArrhenius2ReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef ImpactIonisationArrhenius2ReactionRate_H
#define ImpactIonisationArrhenius2ReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

#include "DynamicList.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ImpactIonisationArrhenius2ReactionRate Declaration
\*---------------------------------------------------------------------------*/

class ImpactIonisationArrhenius2ReactionRate
{
    // Private data

        // Bounds delimiting the different temperature sub-ranges
        DynamicList<scalar> Tbound_;
        
        // List of the standard Arrhenius coefficients
        DynamicList<scalar> A_;
        DynamicList<scalar> beta_;
        DynamicList<scalar> Ta_;
        DynamicList<scalar> C_;
        
        // Smoothing between the sub-ranges
        scalar smoothingTemperatureRange_;
        
    // Private member function
        
        label findBounds
        (
            const scalar xpt
        ) const
        {
            label bMax = Tbound_.size();
            label bCurrent = 0;
            
            if (bMax == 0)
            {
                return bCurrent;
            }
            else
            {
                while (bCurrent < bMax and xpt > Tbound_[bCurrent])
                {
                    bCurrent += 1;   
                }
                return bCurrent;
            }
        }

public:

    // Constructors

        //- Construct from components
        inline ImpactIonisationArrhenius2ReactionRate
        (
            const DynamicList<scalar> Tbound,
            const DynamicList<scalar> A,
            const DynamicList<scalar> beta,
            const DynamicList<scalar> Ta,
            const DynamicList<scalar> C,
            const scalar smoothingTemperatureRange
        );

        //- Construct from Istream
        inline ImpactIonisationArrhenius2ReactionRate
        (
            const speciesTable& species,
            Istream& is
        );

        //- Construct from dictionary
        inline ImpactIonisationArrhenius2ReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "ElectronImpactIonisationArrhenius";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;
        

    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const ImpactIonisationArrhenius2ReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ImpactIonisationArrhenius2ReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
