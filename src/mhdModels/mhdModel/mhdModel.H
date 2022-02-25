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

Namespace
    Foam::mhd

Description
    Namespace for mhd modelling

Class
    Foam::mhd::mhdModel

Description
    Base abtract class for MHD modelling

SourceFiles
    mhdModel.C

\*---------------------------------------------------------------------------*/

#ifndef mhdModel_H
#define mhdModel_H

#include "rho2ReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class rho2ReactionThermo;
class fvMesh;

namespace mhd
{

// Forward declaration of classes
class electricalConductivityModel;


/*---------------------------------------------------------------------------*\
                       Class mhdModel Declaration
\*---------------------------------------------------------------------------*/

class mhdModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the thermophysical model
        const rho2ReactionThermo& thermo_;

        //- MHD model on/off flag
        Switch active_;
        
        //- Hall effect on/off flag
        Switch hallEffect_;
        
        //- Constant Hall parameter;
        scalar constBeta_;

        //- MHD model dictionary
        dictionary coeffs_;


        // References to the mhd sub-models

            //- Electrical conductivity model
            autoPtr<electricalConductivityModel> electricalConductivity_;
            
            
    // Protected Member Functions
    
        //- Initialise models
        void initialise();
        

private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Disallow default bitwise copy construct
        mhdModel(const mhdModel&);

        //- Disallow default bitwise assignment
        void operator=(const mhdModel&);


public:

    //- Runtime type information
    TypeName("mhdModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mhdModel,
            thermo,
            (
                const rho2ReactionThermo& thermo
            ),
            (thermo)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            mhdModel,
            dictionary,
            (
                const dictionary& dict,
                const rho2ReactionThermo& thermo
            ),
            (dict, thermo)
        );


    // Constructors

        //- Null constructor
        mhdModel(const rho2ReactionThermo& thermo);

        //- Construct from components
        mhdModel(const word& type, const rho2ReactionThermo& thermo);

        //- Construct from components
        mhdModel
        (
            const word& type,
            const dictionary& dict,
            const rho2ReactionThermo& thermo
        );


    // Selectors

        //- Return a reference to the selected mhd model
        static autoPtr<mhdModel> New(const rho2ReactionThermo& thermo);

        //- Return a reference to the selected mhd model
        static autoPtr<mhdModel> New
        (
            const dictionary& dict,
            const rho2ReactionThermo& thermo
        );


    //- Destructor
    virtual ~mhdModel();


    // Member Functions

        // Edit

            //- Read mhdProperties dictionary
            virtual bool read() = 0;
            
            //- Update MHD model
            virtual void update(const volVectorField& U) = 0;
            
            //- Return non-const access to the electric field
            virtual volVectorField& E() = 0;
            
            //- Return non-const access to the electric potential field
            virtual volScalarField& elecPot() = 0;
            
            //- Return non-const access to the total magnetic field
            virtual volVectorField& B() = 0;

        // Access
        
            //- Joule heating: source term in the total energy equation
            virtual tmp<volScalarField> jouleHeating
            (
                const volVectorField& U
            ) const = 0;
            
            //- Lorentz force: source term in the momentum equations
            virtual tmp<volVectorField> lorentzForce() const = 0;

            //- Return const access to the state of the MHD on/off switch
            const Switch active() const
            {
                return active_;
            }
            
            //- Return const access to the thermophysical model
            const rho2ReactionThermo& thermo() const
            {
                return thermo_;
            }
            
            //- Return const access to the electrical conductivity field
            virtual const volTensorField& sigma() const = 0;
            
            //- Return const access to the electric field
            virtual const volVectorField& E() const = 0;
            
            //- Return const access to the electric potential field
            virtual const volScalarField& elecPot() const = 0;
            
            //- Return const access to the total magnetic field
            virtual const volVectorField& B() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToMhdRunTimeSelectionTables(model)                                 \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        mhdModel,                                                             \
        model,                                                                \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        mhdModel,                                                             \
        model,                                                                \
        thermo                                                                \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
