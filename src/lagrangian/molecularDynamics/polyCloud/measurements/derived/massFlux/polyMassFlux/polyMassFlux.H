/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    polyMassFlux

Description
    Instantaneous measurements of mass flow rate (here instantaneous is defined 
    at every MD timestep)
    
    Apply to an existing surface 
    
    
    
SourceFiles
    polyMassFlux.C


\*---------------------------------------------------------------------------*/

#ifndef polyMassFlux_H
#define polyMassFlux_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
// #include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyMassFlux Declaration
\*---------------------------------------------------------------------------*/

class polyMassFlux
:
    public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
    
        word fieldName_;

        fileName pathName_;
        fileName nameFile_;

        propertyField fields_;

        label regionId_;
        word faceZoneName_;

        //- this term is not used here
        scalar zoneSurfaceArea_;

        //- ids of molecules
        List<label> molIds_;

        vector fluxDirection_;
    
        scalar cumulativeFlux_;

        DynamicList<scalar> massFlux_;


    // Private Member Functions


 public:

    //- Runtime type information
    TypeName("polyMassFlux");

    // Constructors

        //- Construct from components
        polyMassFlux
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~polyMassFlux();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();
        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        //- write field
        void writeField();

        const propertyField& fields() const;

        void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
