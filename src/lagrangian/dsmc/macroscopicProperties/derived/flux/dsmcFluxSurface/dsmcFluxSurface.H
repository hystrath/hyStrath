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
    dsmcFluxSurface

Description
    Measure average fluxes / flow rates through user-defined faceZones.

    This is based on dsmcMassFluxSurface (which has not been updated and does
    not work in the v1706 version of dsmcFoam+).

    This works in conjunction with the faceTracker object. The faceTracker is
    responsible for recording the correct fluxes on the faces. This class is
    used to calculate the total fluxes / flow rates on a given faceZone (i.e.
    a collection of faces defined by topoSet / setSet).

SourceFiles
    dsmcFluxSurface.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcFluxSurface_H
#define dsmcFluxSurface_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dsmcFluxSurface Declaration
\*---------------------------------------------------------------------------*/

class dsmcFluxSurface
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;

        word fieldName_;

        label regionId_;

        word faceZoneName_;
        scalar faceZoneArea_;

        labelList typeIds_;

        vector fluxDirection_;

        scalar particlesAccumulated_;
        scalar massAccumulated_;
        scalar averagingCounter_;
        label timeIndex_;

        // particle flux in #/m^2/s
        scalarField particleFlux_;
        // particle flow rate in #/s
        scalarField particleFlowRate_;
        // mass flux in kg/m^2/s
        scalarField massFlux_;
        // mass flow rate in kg/s
        scalarField massFlowRate_;

        bool averagingAcrossManyRuns_;

    // Private Member Functions

        // for support of averaging across many runs
        void readIn();
        void writeOut();

        // (re-)calculate the area of all faces in the faceZone
        scalar calculateFaceZoneArea(const labelList& faces);

 public:

    // Runtime type information
    TypeName("dsmcFluxSurface");

    // Constructors

        dsmcFluxSurface
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcFluxSurface();


    // Member Functions

        // Create initial field
        virtual void createField();

        // Calculate fields
        virtual void calculateField();

        // Reset fields (relevant if mesh has changed)
        virtual void resetField();

        // Write out fields
        virtual void writeField();

        virtual void updateProperties(const dictionary& dict);

        // Necessary functions (currently defined as virtual functions in base
        // class, therefore these also have to be defined here)
        virtual scalar translationalT(const label cellI);
        virtual scalar overallT(const label cellI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
