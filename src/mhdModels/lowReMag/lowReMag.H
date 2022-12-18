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

protected:

    // Protected data
    
        //- Const reference to the trans-rotational temperature field
        const volScalarField& T_;
        
        //- Const reference to the electron pressure field
        const volScalarField& pe_;
        
        //- Local copy of the Boltzmann constant
        const scalar localkB_;
        
        //- Local copy of the electric charge
        const scalar localElecCharge_;
        
        //- Total magnetic field [Tesla = kg/s^2/A]
        volVectorField B_;
        
        //- Electric field [V/m = kg*m/s^-3/A]
        volVectorField E_;
        
        //- Electric potential [V = kg*m^2/s^-3/A]
        volScalarField elecPot_;
        
        //- Current density field [A/m^2]
        volVectorField j_;

        //- Disallow default bitwise copy construct
        lowReMag(const lowReMag&);

        //- Disallow default bitwise assignment
        void operator=(const lowReMag&);
        
    // Private Member Functions

        //- Calculate the electric current density
        void updateElectricCurrentDensity(const volVectorField& U);
         
        //- Calculate the Hall tensor
        tmp<volTensorField> hallTensor() const;


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
          virtual void update(const volVectorField& U);
          
          //- Return non-const access to the electric field
          virtual volVectorField& E()
          {
              return E_;
          }
          
          //- Return non-const access to the electric potential field
          virtual volScalarField& elecPot()
          {
              return elecPot_;
          }
          
          //- Return non-const access to the total magnetic field
          virtual volVectorField& B()
          {
              return B_;
          }
          
        // Access
       
          //- Joule heating: source term in the total energy equation
          virtual tmp<volScalarField> jouleHeating
          (
              const volVectorField& U
          ) const;
         
          //- Lorentz force: source term in the momentum equations
          virtual tmp<volVectorField> lorentzForce() const;
         
          //- Return the electrical conductivity field
          virtual volTensorField sigma() const;
         
          //- Return const access to the electric field
          virtual const volVectorField& E() const
          {
              return E_;
          }
      
          //- Return const access to the electric potential field
          virtual const volScalarField& elecPot() const
          {
              return elecPot_;
          }
            
          //- Return const access to the total magnetic field
          virtual const volVectorField& B() const
          {
              return B_;
          }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
