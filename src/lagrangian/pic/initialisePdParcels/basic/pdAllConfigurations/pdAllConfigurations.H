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

#ifndef pdAllConfigurations_H
#define pdAllConfigurations_H


#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"

#include "pdConfiguration.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdAllConfigurations Declaration
\*---------------------------------------------------------------------------*/

class pdAllConfigurations
{

    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        const IOdictionary& pdInitialiseDict_;

		PtrList<entry>configurationList_;
        List<label> ids_;
        List< autoPtr<pdConfiguration> > configurations_;


        
public:

    // Constructors

        //- Null Constructor 
        pdAllConfigurations
        (
//             const polyMesh& mesh,
            const IOdictionary& pdInitialiseDict,
            pdCloud& cloud
        );


    // Member Functions

        //- initial configuration
        void setInitialConfig();

        // Access


            inline List< autoPtr<pdConfiguration> >& configurations();
            inline const List< autoPtr<pdConfiguration> >& configurations() const;


//             inline const List<word>& stateControllersNames() const; 
    


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdAllConfigurationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
