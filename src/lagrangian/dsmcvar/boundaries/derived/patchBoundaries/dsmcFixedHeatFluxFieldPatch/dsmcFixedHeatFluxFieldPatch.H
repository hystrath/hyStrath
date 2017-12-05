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
    dsmcFixedHeatFluxFieldPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcFixedHeatFluxFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcFixedHeatFluxFieldPatch_H
#define dsmcFixedHeatFluxFieldPatch_H

#include "dsmcPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"
#include "timeDataMeas.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcFixedHeatFluxFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcFixedHeatFluxFieldPatch
:
    public dsmcPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;
        
        volScalarField wallTemperature_;

        scalarField EcTot_;

        scalarField EcTotSum_;
        
        scalarField newWallTemperature_;
        
        // user-defined wall heat flux
        scalar desiredHeatFlux_;

        // alpha parameter
        scalar relaxationFactor_;

        label stepCounter_;
        
        label nSamples_;
        
        label nSamplesAverage_;
        
        scalar referenceCp_;
        
        scalar referenceRho_;
        
        scalar referenceU_;
        
        scalar referenceTemperature_;
        
        bool resetAtOutput_;

        void setProperties();
        
        void initialTemperatureEstimate();

public:

    //- Runtime type information
    TypeName("dsmcFixedHeatFluxFieldPatch");

    // Constructors

        //- Construct from components
        dsmcFixedHeatFluxFieldPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcFixedHeatFluxFieldPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td);

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);



        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
