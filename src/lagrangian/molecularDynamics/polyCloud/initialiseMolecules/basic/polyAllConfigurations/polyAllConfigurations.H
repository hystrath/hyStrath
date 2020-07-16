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

Description



\*---------------------------------------------------------------------------*/

#ifndef polyAllConfigurations_H
#define polyAllConfigurations_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "polyConfiguration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyAllConfigurations Declaration
\*---------------------------------------------------------------------------*/

class polyAllConfigurations
{

    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary mdInitialiseDict_;

		PtrList<entry>configurationList_;
        List<label> ids_;
        List< autoPtr<polyConfiguration> > configurations_;

public:

    // Constructors

        //- Null Constructor
        polyAllConfigurations
        (
            const polyMesh& mesh
        );

        //- Constructor for gnemdFOAM
        polyAllConfigurations
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud
        );


    // Member Functions

        //- initial configuration
        void setInitialConfig();

        // Access

        inline List< autoPtr<polyConfiguration> >& configurations();
        inline const List< autoPtr<polyConfiguration> >& configurations() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyAllConfigurationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
