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
    dsmcAxisymmetric

Description

SourceFiles
    dsmcAxisymmetricI.H
    dsmcAxisymmetric.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcAxisymmetric_H
#define dsmcAxisymmetric_H

#include "dsmcCoordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcAxisymmetric Declaration
\*---------------------------------------------------------------------------*/


class dsmcAxisymmetric
:
    public dsmcCoordinateSystem
{

    // Private data

        //- Reference to the DSMC cloud
        dsmcCloud& cloud_;

        //- Axis of revolution
        //  Default: x-axis, i.e., component 0
        label revolutionAxis_;

        //- Polar axis
        //  Default: (revolutionAxis_ + 1)%3
        label polarAxis_;

        //- Angular coordinate
        //  Default: (revolutionAxis_ + 2)%3
        label angularCoordinate_;

        //- Radial extent of the domain
        scalar radialExtent_;

        //- Maximum radial scaling factor
        scalar maxRWF_;

        //- Radial weighting method employed
        //  Default: "cell", other options: "particleAverage"
        word rWMethod_;

        //- Radial weighting factor field
        volScalarField RWF_;


    // Private Member Functions

        //- Clone or delete particles based on their old and new RWF
        void axisymmetricWeighting();

        //- Update the radial weighting factor field
        void updateRWF();

        //- Recalculate and return the radial weighting factor for a cell
        scalar recalculateRWF(const label cellI) const;

        //- Write axisymmetric info
        void writeAxisymmetricInfo() const;

        //- Return the radial extent of the domain
        inline scalar radialExtent() const;

        //- Return the maximum value for RWF
        inline scalar maxRWF() const;



public:

    //- Runtime type information
    TypeName("dsmcAxisymmetric");

    // Constructors

        //- Constructor
        dsmcAxisymmetric
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcAxisymmetric();


    // Member Functions

      // Access

        //- Return const access to the number of real particles represented by
        //  a DSMC parcel
        virtual inline const volScalarField& nParticles() const;

        //- Return const access to the number of real particles represented by
        //  a DSMC parcel for a given cell
        virtual inline scalar nParticles(const label celli) const;

        //- Return const access to the number of real particles represented by
        //  a parcel for a given patch/face couple
        virtual inline scalar nParticles
        (
            const label patchi,
            const label facei
        ) const;

        //- Return the local RWF for a given cell
        virtual inline scalar RWF(const label celli) const;

        //- Return the local RWF for a given face
        virtual inline scalar pRWF
        (
            const label patchi,
            const label facei
        ) const;


      // Check

        //- Check the validity of the coordinate system inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct. Pass init as true if used in the constructor of
        //  the dsmcInitialise application
        virtual void checkCoordinateSystemInputs(const bool init = false);


      // Edit

        //- Evolve function
        virtual void evolve();


      // Write

        //- Write coordinate system info
        virtual void writeCoordinateSystemInfo() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcAxisymmetricI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
