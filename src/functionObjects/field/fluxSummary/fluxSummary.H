/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjects::fluxSummary

Group
    grpFieldFunctionObjects

Description
    This function object calculates the flux across selections of faces.

    Output comprises, per set of faces, the fluxes:
    - positive
    - negative
    - net
    - absolute

Usage
    \verbatim
    fluxSummary1
    {
        type        fluxSummary;
        libs        ("libfieldFunctionObjects.so");
        ...
        write       yes;
        log         yes;
        mode        cellZoneAndDirection;
        cellZoneAndDirection
        (
            (porosity (1 0 0))
        );
        scaleFactor 1.2;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: fluxSummary  | yes         |
        write        | write flux data to file |  no  | yes
        log          | write flux data to standard output | no | yes
        mode         | mode to generate faces to test | yes |
        scaleFactor  | optional factor to scale result | no | 1
    \endtable

    The mode is one of:
    - faceZone
    - faceZoneAndDirection
    - cellZoneAndDirection
    - surface
    - surfaceAndDirection

    Output data is written to files of the form \<timeDir\>/<faceZoneName>.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    fluxSummary.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fluxSummary_H
#define functionObjects_fluxSummary_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "vector.H"
#include "DynamicList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dimensionSet;

namespace functionObjects
{
/*---------------------------------------------------------------------------*\
                         Class fluxSummary Declaration
\*---------------------------------------------------------------------------*/

class fluxSummary
:
    public fvMeshFunctionObject,
    public writeFile
{
public:

    // Public enumerations

        //- Face mode type
        enum modeType
        {
            mdFaceZone,                 //!< face zone
            mdFaceZoneAndDirection,     //!< face zone with prescribed direction
            mdCellZoneAndDirection,     //!< cell zone with prescribed direction
            mdSurface,                  //!< surfMesh
            mdSurfaceAndDirection       //!< surfMesh with prescribed direction
        };

        //- Mode type names
        static const NamedEnum<modeType, 5> modeTypeNames_;


protected:

    // Protected data

        //- Mode for face determination
        modeType mode_;

        //- Scale factor
        scalar scaleFactor_;

        //- Name of flux field, default = phi
        word phiName_;


        // Per-faceZone/surface information

            //- Region (zone/surface) names
            List<word> zoneNames_;

            //- Face IDs
            List<List<label>> faceID_;

            //- Face patch IDs
            List<List<label>> facePatchID_;

            //- Face flip map signs
            List<boolList> faceFlip_;

            //- Output file per face zone
            PtrList<OFstream> filePtrs_;


        //- Tolerance applied when matching face normals
        scalar tolerance_;


    // Private Member Functions

        //- Check if surface mode instead of zone mode
        bool isSurfaceMode() const;

        //- Check flowType (mass or volume).
        //  Return name on success, fatal error on failure.
        word checkFlowType
        (
            const dimensionSet& fieldDims,
            const word& fieldName
        ) const;

        //- Initialise for given surface name
        void initialiseSurface
        (
            const word& surfName,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise for given surface name and direction
        void initialiseSurfaceAndDirection
        (
            const word& surfName,
            const vector& refDir,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise face set from face zone
        void initialiseFaceZone
        (
            const word& faceZoneName,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<List<label>>& faceID,
            DynamicList<List<label>>& facePatchID,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise face set from face zone and direction
        void initialiseFaceZoneAndDirection
        (
            const word& faceZoneName,
            const vector& refDir,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<List<label>>& faceID,
            DynamicList<List<label>>& facePatchID,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Initialise face set from cell zone and direction
        void initialiseCellZoneAndDirection
        (
            const word& cellZoneName,
            const vector& refDir,
            DynamicList<word>& names,
            DynamicList<vector>& dir,
            DynamicList<List<label>>& faceID,
            DynamicList<List<label>>& facePatchID,
            DynamicList<boolList>& faceFlip
        ) const;

        //- Calculate the total area for the surface or derived faceZone
        scalar totalArea(const label idx) const;

        //- Output file header information
        virtual void writeFileHeader
        (
            const word& zoneName,
            const scalar area,
            const vector& refDir,
            Ostream& os
        ) const;

        //- Specialized write for surfaces
        bool surfaceModeWrite();

        //- Disallow default bitwise copy construct
        fluxSummary(const fluxSummary&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fluxSummary&) = delete;


public:

    //- Runtime type information
    TypeName("fluxSummary");


    // Constructors

        //- Construct from Time and dictionary
        fluxSummary
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fluxSummary();


    // Member Functions

        //- Read the field fluxSummary data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fluxSummary
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
