/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
    dsmcDynamicLoadBalancing

Description

SourceFiles
    dsmcDynamicLoadBalancing.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcDynamicLoadBalancing_H
#define dsmcDynamicLoadBalancing_H

#include "vector.H"
#include "fvMesh.H"
#include "dsmcParcel.H"
#include "decompositionMethod.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dsmcDynamicLoadBalancing Declaration
\*---------------------------------------------------------------------------*/

class dsmcDynamicLoadBalancing
{
    // Private data

        Time& time_;

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        IOdictionary dsmcLoadBalanceDict_;

        bool performBalance_;

        bool enableBalancing_;

        scalar balanceUntilTime_;

        scalar originalEndTime_;

        scalar maxImbalance_;

        // allow specifying a limit for the number of time directories that will
        // be backuped up in the resultFolders directory.
        //  - -1 (default): backup all time directories
        //  - N >= 0: backup the N most recent time directories and delete
        //            oldest directories in case the limit is exceeded.
        label limitTimeDirBackups_;

    // Private Member Functions


public:

    friend class dsmcParcel;

    // Constructors

         //- Constructor
        dsmcDynamicLoadBalancing
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );

    // Destructor

        ~dsmcDynamicLoadBalancing();

    // Member Functions

        void update();

        void perform (const label noRefinement = 0);

        void copyPolyMeshToLatestTimeFolder() const;

        void updateProperties();

        // Access

          inline bool performBalance() const
          {
              return performBalance_;
          }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
