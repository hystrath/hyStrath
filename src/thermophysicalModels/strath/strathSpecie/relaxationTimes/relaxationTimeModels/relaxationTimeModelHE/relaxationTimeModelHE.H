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
    Foam::relaxationTimeModelHE

Description
    Abstract base class for heavy-particle - electron (HE) energy transfer.

SourceFiles
   relaxationTimeModelHEI.H
    relaxationTimeModelHE.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationTimeModelHE_H
#define relaxationTimeModelHE_H

#include "fvc.H"

#include "rho2ChemistryModel.H"
#include "turbulenceModel2.H"
#include "multi2ComponentMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                    Class relaxationTimeModelHE Declaration
\*---------------------------------------------------------------------------*/

class relaxationTimeModelHE
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
        
        //- H-E source term in the vib. energy equations [kg/m^3/s]
        volScalarField QHE_;


private:
  
    // Private member functions
      
        //- Disallow copy construct
        relaxationTimeModelHE(const relaxationTimeModelHE&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationTimeModelHE&);


public:
 
    //- Runtime type information
    TypeName("relaxationTimeModelHE");

      
    // Declare run-time constructor selection table
       
#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        relaxationTimeModelHE,
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
        relaxationTimeModelHE
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel2& turbulence
        );


    // Selectors	  
  
        //- Return a reference to the selected multispecies transport model

        static autoPtr<relaxationTimeModelHE> New
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel2& turbulence
        );


    //- Destructor

        virtual ~relaxationTimeModelHE()
        {}


    // Member Functions

        // Access   
           
            //- Return the const source term in vib. energy equations 
            //  due to H-E energy transfer
            inline const volScalarField& QHE() const;
            
            //- Return the table of species
            inline const speciesTable& species() const;
            
            //- Molar weight for a specie given by index
            virtual scalar W(label i) const = 0;


        //- Source term in vib. energy equations due to H-E energy transfer
        //tmp<volScalarField> HERelaxationSource();
        
        //- Correct the source term
        virtual void correct() = 0;

        //- Read thermo2T dictionary
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "relaxationTimeModelHEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
