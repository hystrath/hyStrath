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

\*---------------------------------------------------------------------------*/

#include "referredCell.H"
#include "transform.H"


namespace Foam
{





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

referredCell::referredCell()
:
    boundedBox(),
    sourceCell_(0),
    origProcNo_(Pstream::myProcNo()),
    offset_(0.0),
    translation_(vector::zero),
    neighbouringCells_(),
    rotate_(false),
    rotationPt_(vector::zero),
    R_(tensor::zero)
{}


referredCell::referredCell
(
    const polyMesh& mesh,
    const scalar& offset,
    const label& cell
)
:
    boundedBox(boundingCellPoints(mesh, cell), false),
    sourceCell_(cell),
    origProcNo_(Pstream::myProcNo()),
    offset_(offset),
    translation_(vector::zero),
    neighbouringCells_(),
    rotate_(false),
    rotationPt_(vector::zero),
    R_(tensor::zero)
{
//     setOffsetBoundBox();
}

referredCell::referredCell(const referredCell& rWF)
:
    boundedBox(rWF),
    sourceCell_(rWF.sourceCell_),
    origProcNo_(rWF.origProcNo_),
    offset_(rWF.offset_),
    translation_(rWF.translation_),
    neighbouringCells_(rWF.neighbouringCells_),
    rotate_(rWF.rotate_),
    rotationPt_(rWF.rotationPt_),
    R_(rWF.R_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

referredCell::~referredCell()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

pointField referredCell::boundingCellPoints
(
    const polyMesh& mesh,
    const label& cellI
)
{
    const labelList& points = mesh.cellPoints()[cellI];
    const labelList& faces = mesh.cells()[cellI];

    label sizeOfList = points.size() + faces.size();

    pointField vectorPoints(sizeOfList, vector::zero);

    label counter = 0;

    forAll(points, p)
    {
        vectorPoints[counter] = mesh.points()[points[p]];
        counter++;
    }

    forAll(faces, f)
    {
        vectorPoints[counter] = mesh.faceCentres()[faces[f]];
        counter++;
    }


    return vectorPoints;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void referredCell::initialise
(
    const polyMesh& mesh,
    const scalar& offset,
    const label& cell
)
{
    boundedBox bb (boundingCellPoints(mesh, cell), false);
    
    this->min() = bb.min();
    this->max() = bb.max();

    sourceCell_ = cell;
    offset_ = offset;
}



void referredCell::transformPoint(vector& pt)
{
    if(rotate_)
    {
        pt = (R_ & (pt - rotationPt_) ) + rotationPt_;
    }

    pt += translation_;
}

void referredCell::setNeighbouringCells(const labelList& cellList)
{
    neighbouringCells_.clear();
    neighbouringCells_.setSize(cellList.size());

    forAll(cellList, c)
    {
        neighbouringCells_[c] = cellList[c];
    }
}

void referredCell::setOffsetBoundBox()
{
//     this->expand(offset_);
    this->expandII(offset_);
}

void referredCell::translate(const vector& translate)
{
    this->min() += translate;
    this->max() += translate;
    translation_ += translate;
}


void referredCell::setRotate
(
    const bool& rotate,
    const vector& rotationPt,
//     const vector& rotationAxis,   
//     const scalar& angle
    const tensor& R
)
{
    if(rotate)
    {
        rotate_ = rotate;
        rotationPt_ = rotationPt;
        R_ = R;

        vector vMinNew = (R_ & (this->min() - rotationPt_) ) + rotationPt_;
        vector vMaxNew = (R_ & (this->max() - rotationPt_) ) + rotationPt_;

        // reset min max points
        this->resetBoundedBox(vMinNew, vMaxNew);
    }
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool referredCell::operator==(const referredCell& rhs) const
{
    return
    (
        static_cast<const boundedBox&>(rhs) == static_cast<boundedBox>(*this)
//      && rhs.pts_ == pts_
//      && rhs.patchI_ == patchI_
    );
}


bool referredCell::operator!=(const referredCell& rhs) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Istream& operator>>(Istream& is, referredCell& rWF)
{
    is  >> static_cast<boundedBox&>(rWF) 
        >> rWF.sourceCell_ 
        >> rWF.origProcNo_
        >> rWF.offset_ 
        >> rWF.translation_
        >> rWF.rotate_
        >> rWF.rotationPt_
        >> rWF.R_;

    // Check state of Istream
    is.check
    (
        "Istream& "
        "operator>>(Istream&, referredCell&)"
    );

    return is;
}


Ostream& operator<<(Ostream& os, const referredCell& rWF)
{
    os  << static_cast<const boundedBox&>(rWF) << token::SPACE
        << rWF.sourceCell_ << token::SPACE
        << rWF.origProcNo_<< token::SPACE
        << rWF.offset_ << token::SPACE
        << rWF.translation_<< token::SPACE
        << rWF.rotate_<< token::SPACE
        << rWF.rotationPt_<< token::SPACE
        << rWF.R_;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const referredCell&)"
    );

    return os;
}

const scalar& referredCell::offset() const
{
    return offset_;
}

scalar& referredCell::offset()
{
    return offset_;
}

const label& referredCell::origProcNo() const
{
    return origProcNo_;
}

const label& referredCell::sourceCell() const
{
    return sourceCell_;
}

const labelList& referredCell::neighbouringCells() const
{
    return neighbouringCells_;
}

const vector& referredCell::translation() const
{
    return translation_;
}

const bool& referredCell::rotate() const
{
    return rotate_;
}

}
// ************************************************************************* //
