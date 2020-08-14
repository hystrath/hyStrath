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
    Foam::mhdModel::lowReMag

Description
    Low magnetic Reynolds number MHD model

SourceFiles
    lowReMag.C

\*---------------------------------------------------------------------------*/

#ifndef lowReMag_H
#define lowReMag_H

#include "mhdModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                          Class lowReMag Declaration
\*---------------------------------------------------------------------------*/

class lowReMag
:
    public mhdModel
{
    // Private Member Functions
    
        //- Const reference to the trans-rotational temperature field
        const volScalarField& T_;
        
        //- Const reference to the electron pressure field
        const volScalarField& pe_;
        
        //- Local copy of the Boltzmann constant
        const scalar localkB_;
        
        //- Local copy of the electric charge
        const scalar localElecCharge_;
        
        //- Magnetic field
        volVectorField B_;
        
        //- Electric conductivity field
        volScalarField sigma_;

        //- Disallow default bitwise copy construct
        lowReMag(const lowReMag&);

        //- Disallow default bitwise assignment
        void operator=(const lowReMag&);


public:

    //- Runtime type information
    TypeName("lowReMag");


    // Constructors

        //- Construct from components
        lowReMag(const rho2ReactionThermo& thermo);

        //- Construct from components
        lowReMag(const dictionary& dict, const rho2ReactionThermo& thermo);


    //- Destructor
    virtual ~lowReMag();


    // Member functions

        // Edit

          //- Read mhdProperties dictionary
          virtual bool read();
          
          //- Update MHD model
          virtual void update();
          
       // Access
       
         //- Calculate the current density
         tmp<volVectorField> j(const volVectorField& U) const;
         
         //- Source term component for energy equation TODO
         tmp<volScalarField> Q(const volVectorField& U) const;
         
         //- Source term component for momentum equations TODO
         tmp<volVectorField> F(const volVectorField& U) const;
         
         //- TODO
         tmp<volScalarField> Stuart(const volVectorField& U) const;
         
         //- Calculate the Hall Effect tensor TODO
         tmp<volTensorField> hallCorrection() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
