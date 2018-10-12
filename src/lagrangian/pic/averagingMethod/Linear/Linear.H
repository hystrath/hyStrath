/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::AveragingMethods::Linear

Description
    Linear-mesh lagrangian averaging procedure.

    Uses a first order (linear) inverse weighting scheme. Arbrirary
    quadralaterials are transformed to a unit cude through a trilinear
    interpolation technique. Inverse trilinear calculation calculated using a
    Guass-Newton method.

SourceFiles
    Linear.C

\*---------------------------------------------------------------------------*/

#ifndef Linear_H
#define Linear_H

#include "AveragingMethod.H"

#include "fvCFD.H"
#include "pointPatchField.H"
#include "pointMesh.H"
#include "pointField.H"


#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace AveragingMethods
{

/*---------------------------------------------------------------------------*\
                            Class Linear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Linear
:
    public AveragingMethod<Type>
{
public:

    //- Public typedefs

        //- Gradient type
        typedef typename AveragingMethod<Type>::TypeGrad TypeGrad;


private:

    //- Private data

        //- Volume of the cell-centered regions
        const Field<scalar>& volumeCell_;

        //- Data on the cells
        GeometricField<Type, fvPatchField, volMesh> cellValue_;

        GeometricField<scalar, fvPatchField, volMesh> cellCounter_;

        GeometricField<scalar, pointPatchField, pointMesh> pointVolume_;

        labelListList cellUniqueList_;

        labelList cellTotalSize_;

        mutable vector S_;

        //outer list, cell id, inner list, points owned by cell
        List< vectorField > pL_;

        //outer list, cell id, inner list, point ID in same order as pL_
        List< labelList > pOrder_;

    //- Private static member functions

        //- Return the size of the FieldField parts
        static autoPtr<labelList> size(const fvMesh &mesh);

    //- Private member functions

        //- Order points of hexahedral in clockwise direction
        void order_hexa
        (
        );

        void XtoL
        (
            const point position,
            const label& cI
        ) const;

public:

    //- Runtime type information
    TypeName("Linear");

    //- Constructors

        //- Construct from components
        Linear
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        Linear(const Linear<Type>& am);

        //- Construct and return a clone
        virtual autoPtr<AveragingMethod<Type> > clone() const
        {
            return autoPtr<AveragingMethod<Type> >
            (
                new Linear<Type>(*this)
            );
        }

    //- Destructor
    virtual ~Linear();


    //- Member Functions

        void resetFields();

        //- Add point value to interpolation
        void add
        (
            const point position,
            const tetIndices& tetIs,
            const Type& value
        );

        //- Add point value to interpolation
        List<List<scalar> > stencil
        (
            const point position,
            const tetIndices& tetIs
        );

        //- Interpolate
        Type interpolate
        (
            const point position,
            const tetIndices& tetIs
        );

        //- Interpolate gradient
        TypeGrad interpolateGrad
        (
            const point position,
            const tetIndices& tetIs
        ) const;

        //- Calculate the average
        void average();
        void average(const AveragingMethod<scalar>& weight);

        void updateField
        (
            GeometricField<Type, fvPatchField, volMesh> cField
        );

        //- Return an internal field of the average
        //tmp<Field<Type> > internalField() const;

        tmp<Field<Type> > internalField() const;

        //- Return an internal field of the gradient
        tmp<Field<TypeGrad> > internalFieldGrad() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AveragingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Linear.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
