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
    dsmcIsothermalPressureOutletSpecifiedMolarFraction

Description

    Fully deletion boundary patch: atomisticMolecules get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    dsmcIsothermalPressureOutletSpecifiedMolarFraction.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcIsothermalPressureOutletSpecifiedMolarFraction_H
#define dsmcIsothermalPressureOutletSpecifiedMolarFraction_H

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
                           Class dsmcIsothermalPressureOutletSpecifiedMolarFraction Declaration
\*---------------------------------------------------------------------------*/

class dsmcIsothermalPressureOutletSpecifiedMolarFraction
:
    public dsmcGeneralBoundary
{

private:

    // Private data
    dictionary propsDict_;

    scalar outletPressure_;
    scalar outletTemperature_;
    scalar outletNumberDensity_;
    scalar nTimeSteps_;
    labelList typeIds_;
    vectorField UMean_;
    vectorField outletVelocity_;
    vectorField UCollected_;
    List <scalarField> accumulatedParcelsToInsert_;
    List <scalarField> totalVibrationalEnergy_;
    List <scalarField> vibT_;
    List <scalarField> vDof_;
    List <scalarField> nTotalParcelsSpecies_;
    scalarField moleFractions_;
    scalarField cellVolume_;
//     scalarField outletNumberDensity_;
//     scalarField outletMassDensity_;
    vectorField totalMomentum_;
    scalarField totalMass_;
    scalarField totalMassDensity_;
    scalarField totalNumberDensity_;
    scalarField totalRotationalEnergy_;
    scalarField totalRotationalDof_;
    scalarField nTotalParcels_;
    scalarField nTotalParcelsInt_;
    scalarField mcc_;

    void setProperties();

public:

    //- Runtime type information
    TypeName("dsmcIsothermalPressureOutletSpecifiedMolarFraction");

    // Constructors

        //- Construct from components
    dsmcIsothermalPressureOutletSpecifiedMolarFraction
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    ~dsmcIsothermalPressureOutletSpecifiedMolarFraction();


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
