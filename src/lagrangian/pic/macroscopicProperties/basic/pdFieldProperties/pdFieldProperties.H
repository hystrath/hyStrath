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
    pdFieldProperties

Description

	Stores all the information for the pdFieldProperties models
    defined within the fieldPropertiesDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef pdFieldProperties_H
#define pdFieldProperties_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"

#include "timeData.H"

#include "pdField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdFieldProperties Declaration
\*---------------------------------------------------------------------------*/

class pdFieldProperties
{

    // Private data

        //- The entire field-properties dictionary 
        //  (it contains multiple subDictionaries)

        Time& time_;

        IOdictionary pdFieldPropertiesDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List< autoPtr<pdField> > fields_;

public:

    // Constructors

        //- Null Constructor (mdInitialise)
        pdFieldProperties
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFoam
        pdFieldProperties
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud
        );

        ~pdFieldProperties();

    // Member Functions


        void updateTimeInfo();

        void createFields();

        void calculateFields();

        void writeFields();


        inline List< autoPtr<pdField> >& fields();
        inline const List< autoPtr<pdField> >& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;

    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const pdFieldProperties&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdFieldPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
