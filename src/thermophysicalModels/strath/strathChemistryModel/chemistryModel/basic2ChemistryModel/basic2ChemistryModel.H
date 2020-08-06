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
    Foam::basic2ChemistryModel

Description
    Base class for chemistry models

SourceFiles
    basic2ChemistryModelI.H
    basic2ChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef basic2ChemistryModel_H
#define basic2ChemistryModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "volMesh.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   class basic2ChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class basic2ChemistryModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        basic2ChemistryModel(const basic2ChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const basic2ChemistryModel&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Chemistry activation switch
        Switch chemistry_;

        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Latest estimation of integration step
        DimensionedField<scalar, volMesh> deltaTChem_;

        //- Modified temperature switch when bounding is required 
        //  (cf Scalabrin's PhD thesis (2007), Eq.3.13, p37)
        Switch modifiedTemperature_;

        //- Coefficients for the modified temperature calculation
        PtrList<scalar> modTCoeffs_;

        //- Model for chemistry-vibrational coupling
        word CVModel_;

        // Park TTV model
        
            //- Park's exponent for heavy-particle temperature
            scalar exponentPark_;

            //- Model for the chemistry-vibrational coupling source term
            word ScvModel_;

            //- Model for the preferential chemistry-vibrational coupling
            //  source term, if chosen
            word preferentialModel_;


    // Protected Member Functions

        //- Return non-const access to the latest estimation of integration
        //  step, e.g. for multi-chemistry model
        inline DimensionedField<scalar, volMesh>& deltaTChem();

        //- Correct function - updates due to mesh changes
        void correct();


public:

    //- Runtime type information
    TypeName("basic2ChemistryModel");


    // Constructors

        //- Construct from mesh
        basic2ChemistryModel(const fvMesh& mesh);


    // Selectors

        //- Generic New for each of the related chemistry model
        template<class Thermo>
        static autoPtr<Thermo> New(const fvMesh&);


    //- Destructor
    virtual ~basic2ChemistryModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Chemistry activation switch
        inline Switch chemistry() const;

        //- The number of species
        virtual label nSpecie() const = 0;

        //- The number of reactions
        virtual label nReaction() const = 0;

        //- Return the latest estimation of integration step
        inline const DimensionedField<scalar, volMesh>& deltaTChem() const;

        //- Modified temperature activation switch
        inline Switch modifiedTemperature() const;


        // Functions to be derived in derived classes

            // Fields

                //- Return const access to chemical source terms [kg/m3/s]
                virtual const DimensionedField<scalar, volMesh>& RR
                (
                    const label i
                ) const = 0;

                //- Return access to chemical source terms [kg/m3/s]
                virtual DimensionedField<scalar, volMesh>& RR
                (
                    const label i
                ) = 0;

                //- Return const access to forward chemical source terms
                //  [kg/m3/s]
                virtual const DimensionedField<scalar, volMesh>& RRf
                (
                    const label i
                ) const = 0;

                //- Return access to forward chemical source terms [kg/m3/s]
                virtual DimensionedField<scalar, volMesh>& RRf
                (
                    const label i
                ) = 0;

                //- Return reaction rate of the specieI in reactionI
                virtual tmp<DimensionedField<scalar, volMesh> > calculateRR
                (
                    const label reactionI,
                    const label specieI
                ) const = 0;


            // Chemistry solution

                //- Calculates the reaction rates
                virtual void calculate() = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalar deltaT) = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalarField& deltaT) = 0;

                //- Return the chemical time scale
                virtual tmp<volScalarField> tc() const = 0;

                //- Return source for enthalpy equation [kg/m/s3]
                virtual tmp<volScalarField> Sh() const = 0;

                //- Return chemistry-vibrational coupling source term
                //  in the vibrational equation [kg/m/s3]
                virtual tmp<volScalarField> Scv() const = 0;
                
                //- Return chemistry-vibrational coupling source term
                //  in the vibrational equation for specie i [kg/m/s3]
                virtual tmp<volScalarField> Scv(const label i) const = 0;

                //- Return energy removed from free-electrons during impact
                //  ionisation reactions [kg/m/s3]
                virtual tmp<volScalarField> Seiir() const = 0;
                
                //- Return energy removed from free-electrons during impact
                //  ionisation reactions for specie i [kg/m/s3]
                virtual tmp<volScalarField> Seiir(const label i) const = 0;

                //- Return the heat release, i.e. enthalpy/sec [m2/s3]
                virtual tmp<volScalarField> dQ() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basic2ChemistryModelI.H"

#ifdef NoRepository
#   include "basic2ChemistryModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
