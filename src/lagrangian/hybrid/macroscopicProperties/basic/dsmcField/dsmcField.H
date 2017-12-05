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
    dsmcField

Description

SourceFiles
    dsmcFieldI.H
    dsmcField.C
    dsmcFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcField_H
#define dsmcField_H

#include "IOdictionary.H"
// #include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"

#include "Random.H"
#include "dsmcParcel.H"
#include "timeDataMeas.H"
// #include "propertyField.H"
#include "writeTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcField Declaration
\*---------------------------------------------------------------------------*/

class dsmcField
{

protected:

    // Protected data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        //- subDictionary containing the time properties
        dictionary timeDict_;

        timeDataMeas time_;

//         bool resetFieldsAtOutput_;

//         label resetIndex_;
//         Random rndGen_;
        fileName casePath_;

        fileName timePath_;


    // Private Member Functions

public:

    //- Runtime type information
    TypeName("dsmcField");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcField,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (t, mesh, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        dsmcField
        (
			Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcField> New
        (
			Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

/*
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const scalarField& yData
        );

        void writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const List< Pair<scalar> >& data
        );

        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const List< Pair<scalar> >& data,
            const label& dummy
        );

        //- vector field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData
        );

        //- tensor field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const tensorField& yData
        );*/

    // Destructor

        virtual ~dsmcField();


    // Member Functions


        void updateTime();

        const fileName& casePath() const;
        fileName& casePath();

        const fileName& timePath() const;
        fileName& timePath();

        //- create initial field
 		virtual void createField() = 0;

        //- calculate field
 		virtual void calculateField() = 0;

        // for hybridVolFields
        virtual void setReset(bool&) {};

        // for hybridSubrelax
        virtual void setStartTime(fileName&) {};
        virtual void setSubrelax(bool&) {};
        virtual void setSubReset(bool&) {};

        virtual void setCellZone(word&) {};
        virtual void getFields
        (
            label&,
            volScalarField&,
            volScalarField&,
            volVectorField&,
            volScalarField&,
            volScalarField&,
            PtrList<volScalarField>&,
            PtrList<volScalarField>&
        )
        {};
        virtual void set0Fields
        (
            volScalarField&,
            volScalarField&,
            volVectorField&,
            PtrList<volScalarField>&,
            PtrList<volScalarField>&
        )
        {};

        virtual void setItNumber(label&) {};

        //- write field
        virtual void writeField() = 0;

        void updateBasicFieldProperties(const dictionary&);

        virtual void updateProperties(const dictionary&) = 0;

//         virtual const propertyField& fields() const = 0;
        //- return fields

//         virtual const propertyField& field() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
