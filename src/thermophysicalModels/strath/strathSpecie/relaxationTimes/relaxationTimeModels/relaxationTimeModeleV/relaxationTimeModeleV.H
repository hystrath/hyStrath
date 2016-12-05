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
    Foam::relaxationTimeModeleV

Description
    Abstract base class for e-V relaxation models.

SourceFiles
    relaxationTimeModeleVI.H
    relaxationTimeModeleV.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationTimeModeleV_H
#define relaxationTimeModeleV_H

#include "fvc.H"

#include "rho2ChemistryModel.H"
#include "turbulenceModel2.H"
#include "eVModel.H"
#include "multi2ComponentMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class relaxationTimeModeleV Declaration
\*---------------------------------------------------------------------------*/

class relaxationTimeModeleV
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
        const compressible::turbulenceModel2& turbulence_;

        //- e-V models [s]
        autoPtr<eVModel> taueViModel_;

        //- e-V relaxation time [s]
        PtrList<volScalarField> taueV_;
        
        //- e-V source term in the vib. energy equations [kg/m^3/s]
        PtrList<volScalarField> QeV_;


private:
  
    // Private member functions
      
        //- Disallow copy construct
        relaxationTimeModeleV(const relaxationTimeModeleV&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationTimeModeleV&);


public:
 
    //- Runtime type information
    TypeName("relaxationTimeModeleV");

      
    // Declare run-time constructor selection table
       
#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        relaxationTimeModeleV,
        fvMesh,
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel2& turbulence
        ),
        (thermo, turbulence)
    );
#endif
  

    // Constructors	  

        //- Construct from components
        relaxationTimeModeleV
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel2& turbulence
        );


    // Selectors	  
  
        //- Return a reference to the selected multispecies transport model

        static autoPtr<relaxationTimeModeleV> New
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel2& turbulence
        );


    //- Destructor

        virtual ~relaxationTimeModeleV()
        {}


    // Member Functions

        // Access   
           
            //- Return the const e-V relaxation time field
            inline const volScalarField& taueVi(label i) const;  
            
            //- Return the const source term in vib. energy equations 
            //  due to e-V energy transfer
            inline const volScalarField& QeV(label speciei) const;
            
            //- Return the table of species
            inline const speciesTable& species() const;
            
            //- Return the table of species having a vib. eq. being solved
            inline const speciesTable& solvedVibEqSpecies() const;

            //- Molar weight for a specie given by index
            virtual scalar W(label i) const = 0;


        //- Source term in vib. energy equations due to e-V energy transfer
        //tmp<volScalarField> eVRelaxationSource();
        
        //- Correct the e-V relaxation times and e-V energy transfer
        virtual void correct() = 0;

        //- Read thermo2T dictionary
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "relaxationTimeModeleVI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
