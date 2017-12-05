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

#include "boundedBox.H"
#include "PstreamReduceOps.H"
#include "tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::boundedBox::great(VGREAT);

const Foam::boundedBox Foam::boundedBox::greatBox
(
    point(-VGREAT, -VGREAT, -VGREAT),
    point(VGREAT, VGREAT, VGREAT)
);


const Foam::boundedBox Foam::boundedBox::invertedBox
(
    point(VGREAT, VGREAT, VGREAT),
    point(-VGREAT, -VGREAT, -VGREAT)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundedBox::calculate(const UList<point>& points, const bool doReduce)
{
    if (points.empty())
    {
        min_ = point::zero;
        max_ = point::zero;

        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
    }
    else
    {
        min_ = points[0];
        max_ = points[0];

        for (label i = 1; i < points.size(); i++)
        {
            min_ = ::Foam::min(min_, points[i]);
            max_ = ::Foam::max(max_, points[i]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundedBox::boundedBox(const UList<point>& points, const bool doReduce)
:
    min_(point::zero),
    max_(point::zero)
{
    calculate(points, doReduce);
}


Foam::boundedBox::boundedBox(const tmp<pointField>& points, const bool doReduce)
:
    min_(point::zero),
    max_(point::zero)
{
    calculate(points(), doReduce);
    points.clear();
}


Foam::boundedBox::boundedBox
(
    const UList<point>& points,
    const labelUList& indices,
    const bool doReduce
)
:
    min_(point::zero),
    max_(point::zero)
{
    if (points.empty() || indices.empty())
    {
        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
    }
    else
    {
        min_ = points[indices[0]];
        max_ = points[indices[0]];

        for (label i=1; i < indices.size(); ++i)
        {
            min_ = ::Foam::min(min_, points[indices[i]]);
            max_ = ::Foam::max(max_, points[indices[i]]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundedBox::points() const
{
    tmp<pointField> tPts = tmp<pointField>(new pointField(8));
    //pointField& pt = tPts(); // TODO VINCENT
    pointField pt = tPts();

    pt[0] = min_;                                   // min-x, min-y, min-z
    pt[1] = point(max_.x(), min_.y(), min_.z());    // max-x, min-y, min-z
    pt[2] = point(max_.x(), max_.y(), min_.z());    // max-x, max-y, min-z
    pt[3] = point(min_.x(), max_.y(), min_.z());    // min-x, max-y, min-z
    pt[4] = point(min_.x(), min_.y(), max_.z());    // min-x, min-y, max-z
    pt[5] = point(max_.x(), min_.y(), max_.z());    // max-x, min-y, max-z
    pt[6] = max_;                                   // max-x, max-y, max-z
    pt[7] = point(min_.x(), max_.y(), max_.z());    // min-x, max-y, max-z

    return tPts;
}


void Foam::boundedBox::inflate(const scalar s)
{
    vector ext = vector::one*s*mag();

    min_ -= ext;
    max_ += ext;
}

void Foam::boundedBox::expand(const scalar s)
{
    vector v = max_-min_;
    vector unitV = v/Foam::mag(v);

    min_ -= unitV*s*Foam::sqrt(3.0);
    max_ += unitV*s*Foam::sqrt(3.0);
}

void Foam::boundedBox::expandII(const scalar s)
{
//     Info << "(before) min_ " << min_ << endl;
//     Info << "(before) max_ " << max_ << endl;
    min_ -= vector(1,0,0)*s;
    min_ -= vector(0,1,0)*s;
    min_ -= vector(0,0,1)*s;    
    
    max_ += vector(1,0,0)*s;
    max_ += vector(0,1,0)*s;
    max_ += vector(0,0,1)*s;

//     Info << "(after) min_ " << min_ << endl;
//     Info << "(after) max_ " << max_ << endl;

}

void Foam::boundedBox::contractII(const scalar s)
{
//     Info << "(before) min_ " << min_ << endl;
//     Info << "(before) max_ " << max_ << endl;
    min_ += vector(1,0,0)*s;
    min_ += vector(0,1,0)*s;
    min_ += vector(0,0,1)*s;    
    
    max_ -= vector(1,0,0)*s;
    max_ -= vector(0,1,0)*s;
    max_ -= vector(0,0,1)*s;

//     Info << "(after) min_ " << min_ << endl;
//     Info << "(after) max_ " << max_ << endl;

}

bool Foam::boundedBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (!contains(points[i]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundedBox::contains
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (!contains(points[indices[i]]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundedBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (contains(points[i]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::boundedBox::containsAny
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (contains(points[indices[i]]))
        {
            return true;
        }
    }

    return false;
}


void Foam::boundedBox::resetBoundedBox
(
    const vector& startPoint,
    const vector& endPoint
)
{
    if(startPoint.x() < endPoint.x())
    {
        min_.x() = startPoint.x();
        max_.x() = endPoint.x();
    }
    else
    {
        min_.x() = endPoint.x();
        max_.x() = startPoint.x();
    }
    if(startPoint.y() < endPoint.y())
    {
        min_.y() = startPoint.y();
        max_.y() = endPoint.y();
    }
    else
    {
        min_.y() = endPoint.y();
        max_.y() = startPoint.y();
    }
    if(startPoint.z() < endPoint.z())
    {
        min_.z() = startPoint.z();
        max_.z() = endPoint.z();
    }
    else
    {
        min_.z() = endPoint.z();
        max_.z() = startPoint.z();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundedBox& bb)
{
    if (os.format() == IOstream::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundedBox)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const boundedBox&)");
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundedBox& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundedBox)
        );
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, boundedBox&)");
    return is;
}


// ************************************************************************* //
