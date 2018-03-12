/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::functionObjects::wallBoundedStreamLine

Group
    grpFieldFunctionObjects

Description
    Generates streamline data by sampling a set of user-specified fields along a
    particle track, transported by a user-specified velocity field, constrained
    to a patch.

Usage
    Example of function object specification:
    \verbatim
    wallBoundedStreamLine1
    {
        type            wallBoundedStreamLine;
        libs            ("libfieldFunctionObjects.so");

        writeControl    writeTime;

        setFormat       vtk;
        U               UNear;
        trackForward    yes;

        fields
        (
            UNear
            p
        );

        lifeTime        10000;
        trackLength     1e-3;
        bounds          (0.2 -10 -10)(0.22 10 10);
        cloudName       particleTracks;

        seedSampleSet   patchSeed;
        patchSeedCoeffs
        {
            type        patchSeed;
            patches     (wall);
            axis        x;
            maxPoints   20000;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: wallBoundedStreamLine| yes |
        setFormat    | Output data type        | yes         |
        U            | Tracking velocity field name | yes    |
        fields       | Fields to sample        | yes         |
        lifetime     | Maximum number of particle tracking steps | yes |
        trackLength  | Tracking segment length | no          |
        nSubCycle    | Number of tracking steps per cell | no|
        cloudName    | Cloud name to use       | yes         |
        bounds       | Bounding box to trim tracks | no | greatBox
        seedSampleSet| Seeding method (see below)| yes       |
    \endtable

    Where \c seedSampleSet is typically one of
    \plaintable
        uniform | uniform particle seeding
        cloud   | cloud of points
        patchSeed | seeding via patch faces
        triSurfaceMeshPointSet | points according to a tri-surface mesh
    \endplaintable

Note
    When specifying the track resolution, the \c trackLength OR \c nSubCycle
    option should be used

See also
    Foam::functionObjects::streamLineBase

SourceFiles
    wallBoundedStreamLine.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallBoundedStreamLine_H
#define functionObjects_wallBoundedStreamLine_H

#include "streamLineBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class wallBoundedStreamLine Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLine
:
    public streamLineBase
{
    // Private Member Functions

        //- Find wall tet on cell
        tetIndices findNearestTet
        (
            const PackedBoolList& isWallPatch,
            const point& seedPt,
            const label celli
        ) const;

        //- Disallow default bitwise copy construct
        wallBoundedStreamLine(const wallBoundedStreamLine&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const wallBoundedStreamLine&) = delete;


public:

    //- Runtime type information
    TypeName("wallBoundedStreamLine");


    // Constructors

        //- Construct from Time and dictionary
        wallBoundedStreamLine
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~wallBoundedStreamLine();


    // Member Functions

        //- Read settings
        virtual bool read(const dictionary&);

        //- Do the actual tracking to fill the track data
        virtual void track();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
