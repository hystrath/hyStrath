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

Class
    Foam::eulerianParticle

Description
    Lightweight class to store particle data derived from VOF calculations,
    with special handling for input, output and parallel reduction.

SourceFiles
    eulerianParticle.H
    eulerianParticle.C
    eulerianParticleTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_eulerianParticle_H
#define functionObjects_eulerianParticle_H

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class Istream;
class Ostream;

namespace functionObjects
{
    class eulerianParticle;
}

// Forward declaration of friend functions and operators
Istream& operator>>(Istream&, functionObjects::eulerianParticle&);
Ostream& operator<<(Ostream&, const functionObjects::eulerianParticle&);

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class eulerianParticle Declaration
\*---------------------------------------------------------------------------*/

class eulerianParticle
{

public:

    // Public data

        //- Index of face in faceZone that this particle hits.  Also used to
        //  identify the index of the coarse face of the surface agglomeration
        //  Note: value of -1 used to indicate that the particle has not
        //        been initialised
        label faceIHit;

        //- Volume multiplied by face centres [m4]
        vector VC;

        //- Volume multiplied by velocity [m4/s]
        vector VU;

        //- Volume [m3]
        scalar V;

        //- Injection time - set at collection [s]
        scalar time;


    //- Constructor
    eulerianParticle();


    // Public Member Functions

        //- Write to stream
        void write(Ostream& os) const;

        //- Write to dictionary
        Foam::dictionary writeDict() const;


    // Operators

        friend bool operator==
        (
            const eulerianParticle& a,
            const eulerianParticle& b
        )
        {
            return
                a.faceIHit == b.faceIHit
             && a.VC == b.VC
             && a.VU == b.VU
             && a.V == b.V
             && a.time == b.time;
        }

        friend bool operator!=
        (
            const eulerianParticle& a,
            const eulerianParticle& b
        )
        {
            return !(a == b);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "eulerianParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
