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
    dsmcWangPressureInlet

Description

    Fully deletion boundary patch: atomisticMolecules get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    dsmcWangPressureInlet.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcWangPressureInlet_H
#define dsmcWangPressureInlet_H

#include "dsmcGeneralBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcWangPressureInlet Declaration
\*---------------------------------------------------------------------------*/

class dsmcWangPressureInlet
:
    public dsmcGeneralBoundary
{

private:

    // Private data
        dictionary propsDict_;

    labelList typeIds_;  // I Added.
    scalarField moleFractions_; // I Added.
    scalar inletPressure_;
    scalar inletTemperature_;
    scalar n_; // inlet number density
    scalarField cellVolume_;
    List <scalarField> accumulatedParcelsToInsert_;   // I Added.
    vectorField inletVelocity_;
    vectorField totalMomentum_;
    scalarField totalMass_;
//     scalarField totalMassDensity_;
//     scalarField totalNumberDensity_;
//     scalarField totalRotationalEnergy_;
//     scalarField totalRotationalDof_;
//     List <scalarField> totalVibrationalEnergy_;
//     List <scalarField> vibT_;
//     List <scalarField> vDof_;
//     scalarField velSqrMeanX_;
//     scalarField velSqrMeanY_;
//     scalarField velSqrMeanZ_;
//     scalarField velMeanSqrX_;
//     scalarField velMeanSqrY_;
//     scalarField velMeanSqrZ_;
//     scalarField totalVelSqrMeanX_;
//     scalarField totalVelSqrMeanY_;
//     scalarField totalVelSqrMeanZ_;
//     scalarField totalVelMeanSqrX_;
//     scalarField totalVelMeanSqrY_;
//     scalarField totalVelMeanSqrZ_;
    scalarField nTotalParcels_;
//     scalarField nTotalParcelsInt_;
//     List <scalarField> nTotalParcelsSpecies_;
    scalar nTimeSteps_;
    scalarField mcc_;
    vectorField UMean_;
    vectorField UCollected_;

    void setProperties();

public:

    //- Runtime type information
    TypeName("dsmcWangPressureInlet");

    // Constructors

    //- Construct from components
    dsmcWangPressureInlet
    (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
    );

    // Destructor

    ~dsmcWangPressureInlet();


    // Member Functions

    void initialConfiguration();

    //- calculate any required properties
    void calculateProperties();

    //- control the mols
    void controlParcelsBeforeMove();
    void controlParcelsBeforeCollisions();
    void controlParcelsAfterCollisions();

    //- output data
    void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    //- update properties from a modified dictionary
    void updateProperties(const dictionary&);

    //reset fields if mesh has changed
    void setNewBoundaryFields();



    // Access

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
