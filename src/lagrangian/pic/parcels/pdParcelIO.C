/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is partArt of OpenFOAM.

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

#include "pdParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "pdCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//template <class ParcelType>
Foam::pdParcel::pdParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    U_(vector::zero),
    A_(vector::zero),
    EPot_(0.0),
    ERot_(0.0),
    EVib_(0.0),
    typeId_(-1),
    newParcel_(0),
    classification_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            is >> A_;
            EPot_ = readScalar(is);
            ERot_ = readScalar(is);
            EVib_ = readScalar(is);
            typeId_ = readLabel(is);
            newParcel_ = readLabel(is);
            classification_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
                + sizeof(A_)
                + sizeof(EPot_)
                + sizeof(ERot_)
                + sizeof(EVib_)
                + sizeof(typeId_)
                + sizeof(newParcel_)
                + sizeof(classification_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "Foam::pdParcel::pdParcel"
        "(const Cloud<pdParcel>& cloud, Foam::Istream&), bool"
    );
}


//template <class ParcelType>
void Foam::pdParcel::readFields(Cloud<pdParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<vector> A(c.fieldIOobject("A", IOobject::MUST_READ));
    c.checkFieldIOobject(c, A);

    IOField<scalar> EPot(c.fieldIOobject("EPot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, EPot);

    IOField<scalar> ERot(c.fieldIOobject("ERot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ERot);

    IOField<scalar> EVib(c.fieldIOobject("EVib", IOobject::MUST_READ));
    c.checkFieldIOobject(c, EVib);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newParcel);

    IOField<label> classification(c.fieldIOobject("classification", IOobject::MUST_READ));
    c.checkFieldIOobject(c, classification);

    label i = 0;
    forAllIter(pdCloud, c, iter)
    {
        pdParcel& p = iter();

        p.U_ = U[i];
        p.A_ = A[i];
        p.EPot_ = EPot[i];
        p.ERot_ = ERot[i];
        p.EVib_ = EVib[i];
        p.typeId_ = typeId[i];
        p.newParcel_ = newParcel[i];
        p.classification_ = classification[i];
        i++;
    }
}


//template <class ParcelType>
void Foam::pdParcel::writeFields(const Cloud<pdParcel>& c)
{
    particle::writeFields(c);

    label np =  c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> A(c.fieldIOobject("A", IOobject::NO_READ), np);
    IOField<scalar> EPot(c.fieldIOobject("EPot", IOobject::NO_READ), np);
    IOField<scalar> ERot(c.fieldIOobject("ERot", IOobject::NO_READ), np);
    IOField<scalar> EVib(c.fieldIOobject("EVib", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::NO_READ), np);
    IOField<label> classification(c.fieldIOobject("classification", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(pdCloud, c, iter)
    {
        const pdParcel& p = iter();

        U[i] = p.U();
        A[i] = p.A();
        EPot[i] = p.EPot();
        ERot[i] = p.ERot();
        EVib[i] = p.EVib();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();
        i++;
    }

    U.write();
    A.write();
    EPot.write();
    ERot.write();
    EVib.write();
    typeId.write();
    newParcel.write();
    classification.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//template <class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const pdParcel& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.A()
            << token::SPACE << p.EPot()
            << token::SPACE << p.ERot()
            << token::SPACE << p.EVib()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << p.classification();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U()),
            sizeof(p.U())
            + sizeof(p.A())
            + sizeof(p.EPot())
            + sizeof(p.ERot())
            + sizeof(p.EVib())
            + sizeof(p.typeId())
            + sizeof(p.newParcel())
            + sizeof(p.classification())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::pdParcel&)"
    );

    return os;
}


// ************************************************************************* //
