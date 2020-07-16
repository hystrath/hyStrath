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
    dsmcReflectiveParticleMembranePatch

Description

    RPM method.

    Li J, Liao DY, Yip S (1998) Coupling continuum to molecular-dynamics
    simulation: reflecting particle method and the field estimator.
    Phys Rev E 57(6):7259


SourceFiles
    dsmcReflectiveParticleMembranePatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcReflectiveParticleMembranePatch_H
#define dsmcReflectiveParticleMembranePatch_H

#include "dsmcCyclicBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dsmcReflectiveParticleMembranePatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcReflectiveParticleMembranePatch
:
	public dsmcCyclicBoundary
{

private:

    // Private data

        dictionary propsDict_;

        //- typeIds to be adsorbed
        labelList typeIds_;

        //- Probability a particle has of being specularly reflected when
        //  stricking the membrane
        scalar specularReflectionProb_;

        //- Probability a particle of a given species has of being specularly
        //  reflected when stricking the membrane
        scalarList specularReflectionProbs_;

        //- Cumulative number of specular reflections
        label nReflections_;

        //- Cumulative number of particle passing through the interface
        label nRejections_;


    // Private Member Functions

        void setProperties();

        void readProperties();

        //- Find the position on which the particle initially striked
        vector findOriginalPosition(const dsmcParcel& p, const label fI);


public:

    //- Runtime type information
    TypeName("dsmcReflectiveParticleMembranePatch");

    // Constructors

        //- Construct from components
        dsmcReflectiveParticleMembranePatch
		    (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		    );

    // Destructor

 		    ~dsmcReflectiveParticleMembranePatch();


    // Member Functions

        void calculateProperties();

        //- create an initial configuration
        void initialConfiguration();

        //- control the molecules
        void controlMol
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        );

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
