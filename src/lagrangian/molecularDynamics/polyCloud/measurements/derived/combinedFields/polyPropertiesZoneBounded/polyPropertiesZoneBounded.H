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
    polyPropertiesZoneBounded


Description
    computing properties in a selected zone

    parallelised

    so far the following timeSchemes may be used:

    for instantaneous measurements:

        timeOption          general;
        nSamples        1;
        nAverages       1;

    semi-instantaneous measurements
        timeOption          general;
        nSamples        20;
        nAverages       20;

    time-averaged measurements (outputed every 100 timesteps)
        timeOption      general;
        nSamples        1;
        nAverages       100;
        resetAtOutput   off;
    or
        timeOption          write;
        resetAtOutput       off;

SourceFiles
    polyPropertiesZoneBoundedI.H
    polyPropertiesZoneBounded.C
    polyPropertiesZoneBoundedIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyPropertiesZoneBounded_H
#define polyPropertiesZoneBounded_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "boundedBox.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyPropertiesZoneBounded Declaration
\*---------------------------------------------------------------------------*/

class polyPropertiesZoneBounded
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
        
        fileName pathName_;
        fileName nameFile_;

        propertyField fields_;

        word fieldName_;

        List<boundedBox> boxes_;

        scalar totalVolume_;
        
        scalar mols_;
        scalar mass_;
        vector mom_;
        scalar kE_;
        vector velocityB_;
        tensor virialTensor_;
        tensor kineticTensor_;
        tensor kineticTensor2_;
        scalar dof_;
        scalar keSum_;
        scalar peSum_;
        scalar angularKeSum_;

        List<label> molIds_;

        label timeIndex_;

        scalarField molField_;
        scalarField densityField_;
        scalarField massDensityField_;
        vectorField velocityFieldA_;
        vectorField velocityFieldB_;
        vectorField momentumField_;
        scalarField temperatureField_;
        scalarField temperatureField2_;
        tensorField stressField_;
        scalarField pressureField_;
        tensorField stressField2_;
        scalarField pressureField2_;
        scalarField pKinField_;
        scalarField pVirField_;
        scalarField kEField_;
        scalarField pEField_;
        scalarField angularkEField_;
        scalarField tEField_;

        List<bool> outputField_;

        scalar nAvTimeSteps_;
        bool resetAtOutput_;        
        
    // Private Member Functions

        void setBoundBoxes();

        void writeToStorage();

        bool readFromStorage();

public:

    //- Runtime type information
    TypeName("polyPropertiesZoneBounded");

    // Constructors

        //- Construct from components
        polyPropertiesZoneBounded
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyPropertiesZoneBounded();

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
        
        const propertyField& fields() const;

        //- write field
        void writeField();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
