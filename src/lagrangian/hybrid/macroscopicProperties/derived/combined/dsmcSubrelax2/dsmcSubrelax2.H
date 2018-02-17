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
    dsmcSubrelax2

Description
    Geometric field for temperature, density, velocity, pressure.



SourceFiles
    dsmcSubrelax2I.H
    dsmcSubrelax2.C
    dsmcSubrelax2IO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcSubrelax2_H
#define dsmcSubrelax2_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcSubrelax2 Declaration
\*---------------------------------------------------------------------------*/

class dsmcSubrelax2
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;
        
        scalar mfpReferenceTemperature_;
        word fieldName_;
        label itBeforeRelaxing_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        PtrList<volScalarField> vibrationalTs_;
        PtrList<volScalarField> rhoNType_;
        volVectorField UMean_;
        volScalarField rhoN0_;
        volScalarField rhoM0_;
        volScalarField p0_;
        volScalarField translationalT0_;
        volScalarField rotationalT0_;
        PtrList<volScalarField> vibrationalTs0_;
        PtrList<volScalarField> rhoNType0_;
        volVectorField UMean0_;
        label nTimeSteps_;
        labelList typeIds_;


        // load initial variables
        fileName startingTime_;

        const word pName_;
        const word UName_;
        const word TtName_;
        const word TrName_;
        const List<word> rhoNNames_;
        const List<word> TvNames_;
        const List<scalar> mass_;

        bool stopSamplingReset_;

    // Private Member Functions

        void writeFiles();


public:
    //- Runtime type information
    TypeName("dsmcSubrelax2");

    // Constructors

        //- Construct from components
    dsmcSubrelax2
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    ~dsmcSubrelax2();


    // Member Functions

    //- create initial field
    void createField();

    //- calculate field
    void calculateField();

    // For dsmcSubrelax2
    void setReset
    (
        bool& stopReset
    );
    void setInitialFields
    (
        const List<word>& molecules,
        fileName& timeToReadFrom
    );

    //- write field
    void writeField();

//         const propertyField& fields() const;
    void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

