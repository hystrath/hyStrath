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
    polyControllers

Description

	Stores all the information for the polyControllers models defined within 
    the polyControllersDict, and selects & builds the models automatically.

    

\*---------------------------------------------------------------------------*/

#ifndef polyControllers_H
#define polyControllers_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "polyStateController.H"
#include "polyFluxController.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyControllers Declaration
\*---------------------------------------------------------------------------*/

class polyControllers
{

    // Private data

        Time& time_;

        const polyMesh& mesh_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary polyControllersDict_;

        // number of polyControllers
        label nStateControllers_;
        label nFluxControllers_;

        //- state polyControllers
        PtrList<entry> stateControllersList_;
        List<word> sCNames_;
        List<label> sCIds_;
        List<fileName> sCFixedPathNames_;
        List< autoPtr<polyStateController> > stateControllers_;

        //- flux polyControllers
        PtrList<entry> fluxControllersList_;
        List<word> fCNames_;
        List<label> fCIds_;
        List<fileName> fCFixedPathNames_;
        List< autoPtr<polyFluxController> > fluxControllers_;

        labelList controllersDuringForceComp_;
        
public:

    // Constructors

        //- Null Constructor 
        polyControllers
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFOAM
        polyControllers
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud
        );


        ~polyControllers();

    // Member Functions

        //- initial configuration
        void initialConfig();

        void controlVelocitiesI();
        void controlBeforeMove();
        void controlBeforeForces();
        void controlDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );
        void controlAfterForces();
        void controlVelocitiesII();
        
        
        //- calculate controller properties
        void calculateStateProps();

        //- output controller properties
        void outputStateResults();
        
        // Access

            const label& nStateControllers() const;

            inline List< autoPtr<polyStateController> >& stateControllers();
            inline const List< autoPtr<polyStateController> >& stateControllers() const;

            inline List< autoPtr<polyFluxController> >& fluxControllers();
            inline const List< autoPtr<polyFluxController> >& fluxControllers() const;

            inline const List<word>& stateControllersNames() const; 
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyControllersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
