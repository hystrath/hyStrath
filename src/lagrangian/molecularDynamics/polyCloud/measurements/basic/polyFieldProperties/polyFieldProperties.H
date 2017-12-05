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
    polyFieldProperties

Description

	Stores all the information for the polyFieldProperties models (state and fluxes),
    defined within the polyFieldPropertiesDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef polyFieldProperties_H
#define polyFieldProperties_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "polyField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyFieldProperties Declaration
\*---------------------------------------------------------------------------*/

class polyFieldProperties
{

    // Private data

        //- The entire field-properties dictionary 
        //  (it contains multiple subDictionaries)
        Time& time_;

        IOdictionary polyFieldPropertiesDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List< autoPtr<polyField> > fields_;

        labelList measurementsDuringForceComp_;
        labelList measurementsDuringForceCompSite_;

public:

    // Constructors

        //- Null Constructor (mdInitialise)
        polyFieldProperties
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFoam
        polyFieldProperties
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud
        );

        ~polyFieldProperties();

    // Member Functions


        void createFields();

        void calculateFields();

        void writeFields();

        void measurementsDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

//         void measurementsDuringForceComputation
//         (
//             polyMolecule* molReal,
//             polyReferredMolecule* molRef
//         );

        void measurementsDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

//         void measurementsDuringForceComputationSite
//         (
//             polyMolecule* molReal,
//             polyReferredMolecule* molRef,
//             label sReal,
//             label sRef
//         );

        inline List< autoPtr<polyField> >& fields();
        inline const List< autoPtr<polyField> >& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;

    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const polyFieldProperties&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyFieldPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
