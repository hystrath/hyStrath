/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::eVRelaxationModels::BourdonVervisch

Description
    BourdonVervisch e-V Relaxation model.
    
SourceFiles
    BourdonVervisch.C

\*---------------------------------------------------------------------------*/

#ifndef BourdonVervisch_H
#define BourdonVervisch_H

#include "eVRelaxationModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "fundamentalConstants.H"
#include "physicoChemicalConstants.H"

#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace eVRelaxationModels
{

/*---------------------------------------------------------------------------*\
                           Class BourdonVervisch Declaration
\*---------------------------------------------------------------------------*/

class BourdonVervisch
:
    public eVRelaxationModel
{
    // Private data

        // Labels of the species
        label species1_;
        
        // Name of the species
        word name1_;
        
        // Temperatures for the linear interpolation
        scalar Tlow_, Tmed_, Thigh_;
        
        // Coefficients
        FixedList<scalar,2> BV1_, BV2_, BV3_;
        
        
    // Private member functions
    
        FixedList<scalar,2> findTemperatureRange
        (
            const scalar T
        ) const
        {
            FixedList<scalar, 2> tRange(Tmed_);
            if (T<Tmed_)
            {
                tRange[0] = Tlow_;
            }
            else
            {
                tRange[1] = Thigh_;
            }
            return tRange;
        } 
        
        scalar boundedLinearInterpolation
        (
            const scalar xpt,
            const scalar x1, 
            const scalar x2,
            const scalar y1, 
            const scalar y2
        ) const
        {
            if (xpt<x1)
            {
                return y1;
            }
            else if (xpt>x2)
            {
                return y2;
            }
            else
            {
                return y1 + (xpt-x1)*(y2-y1)/(x2-x1);
            }
        }     
        

public:

    //- Runtime type information
    TypeName("BourdonVervisch");


    // Constructors

        //- construct from components
        BourdonVervisch
        (
            const word& name1,
            const label& lname1,
            const dictionary& dict2T,
            const dictionary& dictThermoPhy,
            const volScalarField& p,
            const volScalarField& Tv
        );


    // Destructor

        ~BourdonVervisch()
        {}


    // Member Functions

        //- V-T relaxation time [s]
        virtual tmp<volScalarField> taueV() const;
         
        //- V-T relaxation time for patch [s]
        virtual tmp<scalarField> taueV
        (
            const label patchi,
            const scalarField& p,
            const scalarField& Tv
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eVRelaxationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
