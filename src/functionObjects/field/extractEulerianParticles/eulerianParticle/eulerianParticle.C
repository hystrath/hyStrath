/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "eulerianParticle.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::functionObjects::eulerianParticle::eulerianParticle()
:
    faceIHit(-1),
    VC(vector::zero),
    VU(vector::zero),
    V(0),
    time(0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const functionObjects::eulerianParticle& p
)
{
    os  << p.faceIHit << token::SPACE
        << p.VC << token::SPACE
        << p.VU << token::SPACE
        << p.V << token::SPACE
        << p.time;

    return os;
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    functionObjects::eulerianParticle& p
)
{
    is  >> p.faceIHit
        >> p.VC
        >> p.VU
        >> p.V
        >> p.time;

    return is;
}


void Foam::functionObjects::eulerianParticle::write(Ostream& os) const
{
    scalar pDiameter = cbrt(6*V/constant::mathematical::pi);
    vector U = VU/(V + ROOTVSMALL);
    vector C = VC/(V + ROOTVSMALL);

    os  << time << token::SPACE
        << faceIHit << token::SPACE
        << C << token::SPACE
        << pDiameter << token::SPACE
        << U << token::SPACE
        << endl;
}


Foam::dictionary Foam::functionObjects::eulerianParticle::writeDict() const
{
    scalar pDiameter = cbrt(6*V/constant::mathematical::pi);
    vector U = VU/(V + ROOTVSMALL);
    vector C = VC/(V + ROOTVSMALL);

    dictionary dict;
    dict.add("time", time);
    dict.add("meshFace", faceIHit);
    dict.add("position", C);
    dict.add("diameter", pDiameter);
    dict.add("U", U);

    return dict;
}


// ************************************************************************* //
