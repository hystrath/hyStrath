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
    Foam::basicMfpModels::variableHardSphere

Description
    variableHardSphere model.
    
SourceFiles
    variableHardSphere.C

\*---------------------------------------------------------------------------*/

#ifndef variableHardSphere_H
#define variableHardSphere_H

#include "basicMfpModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace basicMfpModels
{

/*---------------------------------------------------------------------------*\
                           Class variableHardSphere Declaration
\*---------------------------------------------------------------------------*/

class variableHardSphere
:
    public basicMfpModel
{
    // Private data

        // Molecular weight of the species [g/mol]
        scalar molW_;
        
        // Temperature coefficient of viscosity for the specie
        scalar omega_;  


public:

    //- Runtime type information
    TypeName("variableHardSphere");


    // Constructors

        //- construct from components
        variableHardSphere
        (
            const word& name,
            const label& speciesIndex,
            const dictionary& dict,
            const dictionary& dictThermoPhy,
            const volScalarField& p,
            const volScalarField& Tt
        );


    // Destructor

        ~variableHardSphere()
        {}


    // Member Functions

        //- Mean free path [m]
        virtual tmp<volScalarField> mfp() const;
         
        //- Mean free path for patch [m]
        virtual tmp<scalarField> mfp
        (
            const label patchi,
            const scalarField& p,
            const scalarField& Tt
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace basicMfpModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
