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
    polyMassFluxSurface

Description
    measure the average mass flux at a user-defined surface


NOTES
    To understand this measurement model, one needs to look at how the cloudTracker works, 
    as well as how processor patches are defined on the mesh.

    - cloud tracker fields are generated every time-step (used by measurement models)
      then fields are cleaned.

    - in parallel processing runs a facezone on a processor patch does not get
      divided uniquely. A pair of coupled processor patches will both identify that
      a face belonging to the face occurs. care must therefore be taken to get the flux
      right. For example note how the surface area of the total face zone is calculated.
      Faces on processor patches are given a weighting of 0.5 so that when communication
      occurs the correct zonal surface area is computed. The same is said for polyMolecule
      crossings. In the cloudTracker class only one face on the processor patch gets 
      gets assigned the flux. Thereby synchronisation between processors gives the 
      correct total flux.        

SourceFiles
    polyMassFluxSurface.C


\*---------------------------------------------------------------------------*/

#ifndef polyMassFluxSurface_H
#define polyMassFluxSurface_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyMassFluxSurface Declaration
\*---------------------------------------------------------------------------*/

class polyMassFluxSurface
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

        scalar molsZone_;
        scalar massZone_;
        scalar absMomZone_;
        vector momZone_;
        scalar averagingTime_;
        scalar writeInterval_;

        bool averagingAcrossManyRuns_;
        bool resetAtOutput_;

        scalarField molFluxZone_;
        scalarField massFluxZone_;
        scalarField absMomFluxZone_;
        scalarField momFluxZone_;

        // computing error bars
        bool computeErrorBars_;

        // time-steps between independent samples
        label nSteps_;
        label stepCounter_;
        scalar averagingCounter_;
        scalar meanMassFlux_;

        scalar stdTerm_;


    // Private Member Functions

            void writeToStorage();

            bool readFromStorage();

 public:

    //- Runtime type information
    TypeName("polyMassFluxSurface");

    // Constructors

        //- Construct from components
        polyMassFluxSurface
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyMassFluxSurface();

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
