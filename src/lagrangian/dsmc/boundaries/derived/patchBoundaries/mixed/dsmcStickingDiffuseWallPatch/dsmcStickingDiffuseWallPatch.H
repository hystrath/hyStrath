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
    dsmcStickingDiffuseWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcStickingDiffuseWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcStickingDiffuseWallPatch_H
#define dsmcStickingDiffuseWallPatch_H

#include "dsmcDiffuseWallPatch.H"
#include "dsmcStickingWallPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dsmcStickingDiffuseWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcStickingDiffuseWallPatch
:
	  public dsmcDiffuseWallPatch,
	  public dsmcStickingWallPatch
{

private:

    // Private data

        dictionary propsDict_;

protected:

    // Protected data



    // Protected member functions

        virtual void setProperties();



public:

    //- Runtime type information
    TypeName("dsmcStickingDiffuseWallPatch");

    // Constructors

        //- Construct from components
        dsmcStickingDiffuseWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcStickingDiffuseWallPatch();


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
