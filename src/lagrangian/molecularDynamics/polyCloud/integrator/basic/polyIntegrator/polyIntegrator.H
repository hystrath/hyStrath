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
    polyIntegrator

Description
    Basic/abstract class of a state controller


SourceFiles
    polyIntegratorI.H
    polyIntegrator.C
    polyIntegratorIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyIntegrator_H
#define polyIntegrator_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"
#include "writeTimeData.H"
#include "selectIds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyIntegrator Declaration
\*---------------------------------------------------------------------------*/

class polyIntegrator
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        Time& time_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyIntegrator");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyIntegrator,
            dictionary,
            (
                Time& t,
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (t, molCloud, dict)
        );

    // Constructors

        //- Construct from components
        polyIntegrator
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<polyIntegrator> New
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~polyIntegrator();


    // Member Functions


        //- create an initial configuration
        virtual void init() = 0;

        //- evolve MD simulation
        virtual void evolve() = 0;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
