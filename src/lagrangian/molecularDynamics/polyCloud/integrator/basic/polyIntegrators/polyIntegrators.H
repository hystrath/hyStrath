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
    polyIntegrators

Description

	Stores all the information for the polyIntegrators models defined within 
    the polyIntegratorsDict, and selects & builds the models automatically.

    

\*---------------------------------------------------------------------------*/

#ifndef polyIntegrators_H
#define polyIntegrators_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "polyIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyIntegrators Declaration
\*---------------------------------------------------------------------------*/

class polyIntegrators
{

    // Private data

        Time& time_;

        const polyMesh& mesh_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary controlDict_;

        autoPtr<polyIntegrator> int_;
        
public:

    // Constructors

        //- Null Constructor 
        polyIntegrators
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFOAM
        polyIntegrators
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud
        );


        ~polyIntegrators();

    // Member Functions
        
        autoPtr<polyIntegrator>& integrator();    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
