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
    dsmcFieldProperties

Description

	Stores all the information for the dsmcFieldProperties models
    defined within the fieldPropertiesDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef dsmcFieldProperties_H
#define dsmcFieldProperties_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"

#include "timeData.H"

#include "dsmcField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcFieldProperties Declaration
\*---------------------------------------------------------------------------*/

class dsmcFieldProperties
{

    // Private data

        //- The entire field-properties dictionary 
        //  (it contains multiple subDictionaries)

        Time& time_;

        IOdictionary dsmcFieldPropertiesDict_;

        PtrList<entry> fieldList_;
        List<word> fieldNames_;
        List<label> fieldIds_;
        List< autoPtr<dsmcField> > fields_;

public:

    // Constructors

        //- Null Constructor (mdInitialise)
        dsmcFieldProperties
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFoam
        dsmcFieldProperties
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );

        ~dsmcFieldProperties();

    // Member Functions


        void updateTimeInfo();

        void createFields();

        void calculateFields();

        void setReset // For hybridVolFields
        (
            bool& stopReset
        );

        void setStartTime // For hybridSubrelax
        (
            fileName& startTimeName
        );

        void setSubReset // For hybridSubrelax
        (
            bool& subreset
        );

        void setSubrelax // For hybridSubrelax
        (
            bool& subrelax
        );

        void setCellZone // for hybridSubrelax
        (
            word& zoneToModify
        );

        void setItNumber // for hybridSubrelax
        (
            label& itNumber
        );

        void getFields
        (
            label& nIter,
            volScalarField& rho,
            volScalarField& Ttr,
            volVectorField& U,
            volScalarField& Tt,
            volScalarField& Tr,
            PtrList<volScalarField>& Tv,
            PtrList<volScalarField>& Y
        );

        void set0Fields
        (
            volScalarField& rho,
            volScalarField& Ttr,
            volVectorField& U,
            PtrList<volScalarField>& Tv,
            PtrList<volScalarField>& Y
        );

        void writeFields();


        inline List< autoPtr<dsmcField> >& fields();
        inline const List< autoPtr<dsmcField> >& fields() const;
        inline const List<word>& fieldNames() const;
        inline const List<label>& fieldIds() const;

    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const dsmcFieldProperties&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcFieldPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
