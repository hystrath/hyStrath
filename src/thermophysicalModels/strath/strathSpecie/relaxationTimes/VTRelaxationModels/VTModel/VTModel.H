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
    Foam::VTModel

Description
    V-T relaxation model for multispecies binary diffusion coefficients 

SourceFiles
    VTModel.C

\*---------------------------------------------------------------------------*/

#ifndef VTModel_H
#define VTModel_H

#include "VTRelaxationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class VTModel Declaration
\*---------------------------------------------------------------------------*/

class VTModel
{

protected:

    // Protected data
    
        //- Two-temperature model dictionary
        const IOdictionary dict2T_;
        
        //- Thermophysical properties dictionary
        const IOdictionary dictThermoPhy_;

        //- Species having a vib. eq. being solved
        const wordList& solvedVibEqSpecies_;
        
        //- Species
        const wordList& species_;
        
        //- Pressure [Pa]
        const volScalarField& p_;

        //- Translational temperature [K]
        const volScalarField& Tt_;
        
        //- Vibrational temperatures [K]
        const PtrList<volScalarField>& Tv_;
        
        //- Number densities [m^-3]
        const PtrList<volScalarField>& nD_;
        

        //- V-T relaxation models 
        PtrList<VTRelaxationModel> tauVTijModels_;   

        //- Interspecies V-T relaxation time
        PtrList<volScalarField> tauVTij_;         


    // Private Member Functions

        //- Disallow copy construct
        VTModel(const VTModel&);

        //- Disallow default bitwise assignment
        void operator=(const VTModel&);


public:

    //- Runtime type information
    TypeName("VTModel");

    // Constructors

        //- Construct from components
        VTModel
        (
            const word& dict2T,
            const word& dictThermoPhy,
            const wordList& solvedVibEqSpecies,
            const wordList& species,
            const volScalarField& p,
            const volScalarField& Tt,
            const PtrList<volScalarField>& Tv,
            const PtrList<volScalarField>& nD 
        );


    // Destructor

        virtual ~VTModel()
        {}


    // Member Functions

        inline const volScalarField& tauVT(label i, label j) const
        {
            label k = solvedVibEqSpecies_.size()*i+j;
            return tauVTij_[k];
        }

        //- Update the inter-species V-T relaxation times
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
