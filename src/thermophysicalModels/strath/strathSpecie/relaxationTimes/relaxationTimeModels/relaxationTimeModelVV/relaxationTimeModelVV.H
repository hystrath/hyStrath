/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Foam::relaxationTimeModelVV

Description
    Abstract base class for V-V relaxation models.

SourceFiles
   relaxationTimeModelVVI.H
    relaxationTimeModelVV.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationTimeModelVV_H
#define relaxationTimeModelVV_H

#include "fvc.H"

#include "rho2ChemistryModel.H"
#include "turbulentFluidThermoModel.H"
#include "VVModel.H"
#include "multi2ComponentMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class relaxationTimeModelVV Declaration
\*---------------------------------------------------------------------------*/

class relaxationTimeModelVV
:

    public IOdictionary
{
    
protected:  
  
    // Protected data

        //- Mesh
        const fvMesh& mesh_;

        //- Reference to the thermodynamic package
        rho2ReactionThermo& thermo_;

        //- Reference to the turbulence model
        const compressibleTurbulenceModel& turbulence_;

        //- Interspecies V-V models
        autoPtr<VVModel> tauVVijModel_;

        //- V-V relaxation time [s]
        PtrList<volScalarField> tauVV_;
        
        //- V-V source term in the vib. energy equations [kg/m^3/s]
        PtrList<volScalarField> QVV_;


private:
  
    // Private member functions
      
        //- Disallow copy construct
        relaxationTimeModelVV(const relaxationTimeModelVV&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationTimeModelVV&);


public:
 
    //- Runtime type information
    TypeName("relaxationTimeModelVV");

      
    // Declare run-time constructor selection table
       
#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        relaxationTimeModelVV,
        fvMesh,
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        ),
        (thermo, turbulence)
    );
#endif
  

    // Constructors	  

        //- Construct from components
        relaxationTimeModelVV
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


    // Selectors	  
  
        //- Return a reference to the selected multispecies transport model

        static autoPtr<relaxationTimeModelVV> New
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


    //- Destructor

        virtual ~relaxationTimeModelVV()
        {}


    // Member Functions

        // Access   
           
            //- Return the const interspecies V-V relaxation time field
            inline const volScalarField& tauVVij(label i, label j) const;  
            
            //- Return the const source term in vib. energy equations 
            //  due to V-V energy transfer
            inline const volScalarField& QVV(label speciei) const;
            
            //- Return the table of species
            inline const speciesTable& species() const;
            
            //- Return the table of species having a vib. eq. being solved
            inline const speciesTable& solvedVibEqSpecies() const;

            //- Molar weight for a specie given by index
            virtual scalar W(label i) const = 0;


        //- Source term in vib. energy equations due to V-V energy transfer
        //tmp<volScalarField> VVRelaxationSource();
        
        //- Correct the V-V relaxation times and V-V energy transfer
        virtual void correct() = 0;

        //- Read thermo2T dictionary
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "relaxationTimeModelVVI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
