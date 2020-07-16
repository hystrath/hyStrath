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
    dsmcCLLWallFieldPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcCLLWallFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcCLLWallFieldPatch_H
#define dsmcCLLWallFieldPatch_H

#include "dsmcPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcCLLWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcCLLWallFieldPatch
:
	public dsmcPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;


        void setProperties();



public:

    scalar normalAccommodationCoefficient_;

    scalar tangentialAccommodationCoefficient_;

	scalar rotationalEnergyAccommodationCoefficient_;

    //- boundary temperature
    volScalarField boundaryT_;

    //- boundary velocity
    volVectorField boundaryU_;

    //- Runtime type information
    TypeName("dsmcCLLWallFieldPatch");

    // Constructors

    //- Construct from components
    dsmcCLLWallFieldPatch
    (
	Time& t,
	const polyMesh& mesh,
	dsmcCloud& cloud,
	const dictionary& dict
    );

    // Destructor

    ~dsmcCLLWallFieldPatch();


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
