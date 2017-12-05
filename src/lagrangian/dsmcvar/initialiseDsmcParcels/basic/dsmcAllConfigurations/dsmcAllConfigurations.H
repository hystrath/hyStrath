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

Description

    

\*---------------------------------------------------------------------------*/

#ifndef dsmcAllConfigurations_H
#define dsmcAllConfigurations_H


#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"

#include "dsmcConfiguration.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcAllConfigurations Declaration
\*---------------------------------------------------------------------------*/

class dsmcAllConfigurations
{

    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        const IOdictionary& dsmcInitialiseDict_;

        PtrList<entry>configurationList_;
        List<label> ids_;
        List< autoPtr<dsmcConfiguration> > configurations_;


        
public:

    // Constructors

        //- Null Constructor 
        dsmcAllConfigurations
        (
//             const polyMesh& mesh,
            const IOdictionary& dsmcInitialiseDict,
            dsmcCloud& cloud
        );


    // Member Functions

        //- initial configuration
        void setInitialConfig();

        // Access


        inline List< autoPtr<dsmcConfiguration> >& configurations();
        inline const List< autoPtr<dsmcConfiguration> >& configurations() const;


//             inline const List<word>& stateControllersNames() const; 
    


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcAllConfigurationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
