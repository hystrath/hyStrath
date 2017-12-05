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
    dsmcReflectiveParticleMembranePatch

Description

    RPM technique Li et al.


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

        scalar p_;

        label nReflections_;
        label nRejections_;



    // Private Member Functions

        void readProperties();


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
            dsmcParcel& mol,
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
