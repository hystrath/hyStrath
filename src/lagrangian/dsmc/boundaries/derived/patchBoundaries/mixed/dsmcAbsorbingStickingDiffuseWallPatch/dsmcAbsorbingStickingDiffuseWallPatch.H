/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    dsmcAbsorbingStickingDiffuseWallPatch

Description
    Wall interaction where a particle can be adsorbed according to a
    probability. Othwerwise, its properties are drawn from a Maxwellian
    distribution corresponding to a specified temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict

SourceFiles
    dsmcAbsorbingStickingDiffuseWallPatch.C
\*---------------------------------------------------------------------------*/

#ifndef dsmcAbsorbingStickingDiffuseWallPatch_H
#define dsmcAbsorbingStickingDiffuseWallPatch_H

#include "dsmcAbsorbingWallPatch.H"
#include "dsmcStickingWallPatch.H"
#include "dsmcDiffuseWallPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class dsmcAbsorbingStickingDiffuseWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcAbsorbingStickingDiffuseWallPatch
:
    public dsmcDiffuseWallPatch,
    public dsmcAbsorbingWallPatch,
    public dsmcStickingWallPatch
{

private:

    // Private data



protected:

    // Protected member functions



public:

    //- Runtime type information
    TypeName("dsmcAbsorbingStickingDiffuseWallPatch");

    // Constructors

        //- Construct from components
        dsmcAbsorbingStickingDiffuseWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcAbsorbingStickingDiffuseWallPatch();


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        );

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        void writeField();

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&);


    // Access


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
