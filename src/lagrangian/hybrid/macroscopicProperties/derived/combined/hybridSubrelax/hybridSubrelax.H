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
    hybridSubrelax

Description
    Geometric field for temperature, density, velocity, pressure.



SourceFiles
    hybridSubrelaxI.H
    hybridSubrelax.C
    hybridSubrelaxIO.C

\*---------------------------------------------------------------------------*/

#ifndef hybridSubrelax_H
#define hybridSubrelax_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hybridSubrelax Declaration
\*---------------------------------------------------------------------------*/

class hybridSubrelax
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;
        
        scalar mfpReferenceTemperature_;
        word fieldName_;
        word zoneToModify_;
        label itBeforeRelaxing_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
//        PtrList<volScalarField> translationalTs_;
//        PtrList<volScalarField> rotationalTs_;
        PtrList<volScalarField> vibrationalTs_;
        PtrList<volScalarField> rhoNType_;
        volVectorField UMean_;
        volScalarField translationalT0_;
        volScalarField rotationalT0_;
//        PtrList<volScalarField> translationalTs0_;
//        PtrList<volScalarField> rotationalTs0_;
        PtrList<volScalarField> vibrationalTs0_;
        PtrList<volScalarField> rhoNType0_;
        volVectorField UMean0_;
        scalarField rhoNMean_;
        scalarField rhoMMean_;
        scalarField rhoNMeanInt_;
        scalarField linearKEMean_;
        vectorField momentumMean_;
        scalarField rotationalEMean_;
        scalarField rotationalDofMean_;
        scalarField muu_;
        scalarField muv_;
        scalarField muw_;
        scalarField mvv_;
        scalarField mvw_;
        scalarField mww_;
        scalarField mcc_;
        scalarField mccu_;
        scalarField mccv_;
        scalarField mccw_;
        scalarField eu_;
        scalarField ev_;
        scalarField ew_;
        scalarField e_;
        List<scalarField> vibrationalETotal_;
        List<scalarField> nParcels_;
        label nTimeSteps_;
        labelList typeIds_;


        // load initial variables
        fileName startingTime_;

/*        const word pName_;
        const word UName_;
        const word TtName_;
        const word TrName_;
        const List<word> rhoNNames_;
        const List<word> TvNames_;
        const List<scalar> mass_;*/

//        bool meshChangeing_;
        bool subResetting_;
        bool subrelaxation_;

    // Private Member Functions

        void writeFiles();


public:
    //- Runtime type information
    TypeName("hybridSubrelax");

    // Constructors

        //- Construct from components
    hybridSubrelax
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    ~hybridSubrelax();


    // Member Functions

    //- create initial field
    void createField();

    //- calculate field
    void calculateField();

    // For hybridSubrelax
    void setStartTime
    (
        fileName& startTimeName
    );

    void setSubReset
    (
        bool& subReset
    );

    void setSubrelax
    (
        bool& subrelax
    );

    void setCellZone
    (
        word& zoneToModify
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

    void setInitialFields
    (
        const List<word>& molecules,
        fileName& timeToReadFrom
    );

    void setItNumber
    (
        label& itValue
    );

    //- write field
    void writeField();

//         const propertyField& fields() const;
    void updateProperties(const dictionary& dict);

    template<typename T>
    void assignFields
    (
        const fvMesh& mesh1,
        T& field1,
        T& field2
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

