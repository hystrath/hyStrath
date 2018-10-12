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
    pdControllers

Description

	Stores all the information for the pdControllers models defined within 
    the pdControllersDict, and selects & builds the models automatically.

    

\*---------------------------------------------------------------------------*/

#ifndef pdControllers_H
#define pdControllers_H


#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
#include "timeData.H"

#include "pdStateController.H"
#include "pdFluxController.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdControllers Declaration
\*---------------------------------------------------------------------------*/

class pdControllers
{

    // Private data

        //- The entire dictionary (containing multiple subDictionaries)


        Time& time_;

        IOdictionary pdControllersDict_;

        // number of pdControllers
        label nStateControllers_;
        label nFluxControllers_;

        //- state pdControllers
		PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List< autoPtr<pdStateController> > stateControllers_;

        //- flux pdControllers
		PtrList<entry> fluxControllersList_;
        List<word> fCNames_;
        List<label> fCIds_;
        List<fileName> fCFixedPathNames_;
		List< autoPtr<pdFluxController> > fluxControllers_;

        
public:

    // Constructors

        //- Null Constructor 
        pdControllers
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFOAM
        pdControllers
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud
        );


        ~pdControllers();

    // Member Functions

        //- initial configuration
        void initialConfig();

        //- update time information
        void updateTimeInfo();

        //- different control stages 
        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        //- calculate controller properties
        void calculateProps();

        //- output controller properties
        void outputResults();
        
        // Access

            const label& nStateControllers() const;

            inline List< autoPtr<pdStateController> >& stateControllers();
            inline const List< autoPtr<pdStateController> >& stateControllers() const;

            inline List< autoPtr<pdFluxController> >& fluxControllers();
            inline const List< autoPtr<pdFluxController> >& fluxControllers() const;

            inline const List<word>& stateControllersNames() const; 
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdControllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
