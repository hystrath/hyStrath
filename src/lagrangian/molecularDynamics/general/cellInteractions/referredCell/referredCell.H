/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::referredCell

Description


SourceFiles
    referredCellI.H
    referredCell.C

\*---------------------------------------------------------------------------*/

#ifndef referredCell_H
#define referredCell_H

// #include "face.H"
// #include "pointField.H"
#include "boundedBox.H"
#include "fvMesh.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class referredCell;
Istream& operator>>(Istream&, referredCell&);
Ostream& operator<<(Ostream&, const referredCell&);


/*---------------------------------------------------------------------------*\
                        Class referredCell Declaration
\*---------------------------------------------------------------------------*/

class referredCell
:
    public boundedBox
{
    // Private data

        // original cell label (on original mesh)
        label sourceCell_;

        // original processor number
        label origProcNo_;

        // offset distance to expand boundBox
        scalar offset_;

        // translation vector for referrals
        vector translation_;

        // neighbouring cells of a referred cell
        labelList neighbouringCells_;

        // rotations of referredCell
        bool rotate_;
        vector rotationPt_;

        // rotational tensor
        tensor R_;
    



        pointField boundingCellPoints
        (
            const polyMesh& mesh,
            const label& cell
        );

//         void resetBoundBox
//         (
//             const vector& vMinNew,
//             const vector& vMaxNew
//         );


public:

    // Constructors

        //- Construct null
        referredCell();

        //- Construct from components
        referredCell
        (
            const polyMesh& mesh,
            const scalar& offset,
            const label& cell
        );

        //- Construct as copy
        referredCell(const referredCell&);


    //- Destructor
    ~referredCell();


    // Member Functions
    
        void initialise
        (
            const polyMesh& mesh,
            const scalar& offset,
            const label& cell
        );

        // Access

            void setOffsetBoundBox();

//             void removeOffsetBoundBox();

            void translate(const vector& translate);

            void setRotate
            (
                const bool& rotate,
                const vector& rotationPt,
//                 const vector& rotationAxis,   
//                 const scalar& angle
                const tensor& R
            );

//             void transformCell();

            void transformPoint(vector& pt);

            void setNeighbouringCells(const labelList& cellList);

            //- Return access to the stored points
            const scalar& offset() const;
            scalar& offset();

            const label& origProcNo() const;

            const label& sourceCell() const;

            const labelList& neighbouringCells() const;

            const vector& translation() const;
            
            const bool& rotate() const;            
            //- Return non-const access to the stored points
//             inline pointField& points();

            //- Return access to the patch index
//             inline label patchIndex() const;

            //- Return non-const access to the patch index
//             inline label& patchIndex();


    // Member Operators

        bool operator==(const referredCell&) const;
        bool operator!=(const referredCell&) const;

    // IOstream Operators

        friend Istream& operator>>(Istream&, referredCell&);
        friend Ostream& operator<<(Ostream&, const referredCell&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
