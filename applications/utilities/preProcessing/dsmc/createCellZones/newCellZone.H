/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    newCellZone

Description

SourceFiles
    newCellZoneI.H
    newCellZone.C
    newCellZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef newCellZone_H
#define newCellZone_H

#include "IOdictionary.H"
#include "polyMesh.H"
#include "boundedBox.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class newCellZone Declaration
\*---------------------------------------------------------------------------*/

class newCellZone
{
    // Private data

        const polyMesh& mesh_;
        
        const dictionary& dict_;
        
        word regionName_;

        labelField cells_;

        word option_;

        bool writeCellSet_;
        
        
    // Private Member Functions


public:

    // Constructors

        //- Construct from properties
        newCellZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        ~newCellZone();

    // Member Functions

        void setZone();

        // Access

        inline const labelField& cells() const;

        inline const word& name() const;
        
        inline const bool& writeCellSet() const;        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "newCellZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
