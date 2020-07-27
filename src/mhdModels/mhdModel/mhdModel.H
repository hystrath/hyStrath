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

        //- MHD model dictionary
        dictionary coeffs_;


        // References to the mhd sub-models

            //- Electrical conductivity model
            autoPtr<electricalConductivityModel> electricalConductivity_;

private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Initialise
        void initialise();

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
            virtual void update() = 0;

        // Access
        
            //- Source term component for the energy equation
            virtual tmp<volScalarField> Q(const volVectorField& U) const = 0;
            
            //- Source term component for the momentum equations
            virtual tmp<volVectorField> F(const volVectorField& U) const = 0;
            
            //- TODO
            virtual tmp<volScalarField> Stuart
            (
                const volVectorField& U
            ) const = 0;

            //- Return the state of the MHD on/off switch
            const Switch active() const
            {
                return active_;
            }
            
            //- Access function to the thermophysical model
            const rho2ReactionThermo& thermo() const
            {
                return thermo_;
            }
            
            //- Access function to the electrical conductivity model
            const electricalConductivityModel& electricalConductivity() const;
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
