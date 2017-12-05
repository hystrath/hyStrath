/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::reducedUnits

Description

SourceFiles
    reducedUnitsI.H
    reducedUnits.C
    reducedUnitsIO.C

\*---------------------------------------------------------------------------*/

#ifndef reducedUnits_H
#define reducedUnits_H

#include "scalar.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class reducedUnits Declaration
\*---------------------------------------------------------------------------*/

class reducedUnits
{
    // Private data

        // Reduced units

                bool reducedUnits_;

                bool outputSI_;

            // Fundamental values

                scalar refLength_;

                scalar refTime_;

                scalar refMass_;

                scalar refCharge_;


            // Derived values

                scalar refEnergy_;

                scalar refTemp_;

                scalar refForce_;

                scalar refVelocity_;

                scalar refVolume_;

                scalar refPressure_;

                scalar refMassDensity_;

                scalar refNumberDensity_;

                scalar refHeatFlux_;

                scalar refAmpere_;




    // Private Member Functions

        void calcRefValues();

        //- Disallow default bitwise copy construct
        reducedUnits(const reducedUnits&);

        //- Disallow default bitwise assignment
        void operator=(const reducedUnits&);


public:

    // Static data members

        //- Static data 
        static const scalar kb;
        static const scalar elementaryCharge;
        static const scalar vacuumPermittivity;

    // Constructors

        //- Construct with no argument
        reducedUnits();

        //- Construct from components
        reducedUnits
        (
            scalar refLength,
            scalar refTime,
            scalar refMass,
            scalar refCharge
        );

        reducedUnits
        (
            label unity
        );        
        
        //- Construct from dictionary
        reducedUnits(const IOdictionary& reducedUnitsDict);
        reducedUnits(Time& runTime, const polyMesh& mesh);


    // Destructor

        ~reducedUnits();


    // Member Functions

        void setRefValues
        (
            scalar refLength,
            scalar refTime,
            scalar refMass,
            scalar refCharge
        );

        void setRefValues(const IOdictionary& reducedUnitsDict);


        // Access

            inline bool runReducedUnits() const;

            inline bool outputSIUnits() const;

            inline scalar refLength() const;
//             inline scalar& refLength();
            
            inline scalar refTime() const;
//             inline scalar& refTime();
            
            inline scalar refMass() const;
//             inline scalar& refMass();
            
            inline scalar refCharge() const;
//             inline scalar& refCharge();
            
            inline scalar refTemp() const;

            inline scalar refEnergy() const;

            inline scalar refForce() const;

            inline scalar refVelocity() const;

            inline scalar refVolume() const;

            inline scalar refPressure() const;

            inline scalar refMassDensity() const;

            inline scalar refNumberDensity() const;


            inline scalar refMolFlux() const;

            inline scalar refMassFlux() const;

            inline scalar refHeatFlux() const;

            inline scalar refAmpere() const;


            //- calculated properties   

            inline scalar kB() const;

            inline scalar epsilonPermittivity() const;

            inline scalar epsilonCharge() const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const reducedUnits&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reducedUnitsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
