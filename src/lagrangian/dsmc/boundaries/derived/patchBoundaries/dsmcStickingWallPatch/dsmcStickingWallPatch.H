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
    dsmcStickingWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcStickingWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcStickingWallPatch_H
#define dsmcStickingWallPatch_H

#include "dsmcPatchBoundary.H"

#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dsmcStickingWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcStickingWallPatch
:
	  public virtual dsmcPatchBoundary
{

private:

    // Private data

        dictionary propsDict_;

protected:

    // Protected data

        //- typeIds to be adsorbed
        labelList typeIds_;

        //- Adsorption probabilities for all species
        scalarList adsorptionProbs_;

        //- Residence time for a particle stuck to the wall
        scalarList residenceTimes_;

        //- Number of stuck particles per face
        scalarList nStuckParticles_;

        //- Saturation limit per face [#]
        scalarList saturationLimit_;


    // Protected member functions

        virtual void setProperties();

        //- Read the number of stuck particles on the patch
        void readPatchField();

        //- Check if adsorption can be still operated on this face
        bool isNotSaturated(const label facei);

        //- Adsorb a particle
        void adsorbParticle
        (
            dsmcParcel& p,
            scalar localTemperature = 0
        );

        //- Check if a particle p should be released
        void testForDesorption(dsmcParcel& p);

        //- Update boundaryFluxMeasurements after the particle desorbed
        void measurePropertiesAfterDesorption
        (
            dsmcParcel& p,
            const scalar& heatOfReaction = 0
        );



public:

    //- Runtime type information
    TypeName("dsmcStickingWallPatch");

    // Constructors

        //- Construct from components
        dsmcStickingWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcStickingWallPatch();


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        ) = 0;

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

        //- Return const access to all adsorption probabilities
        const scalarList& adsorptionProbs() const
        {
            return adsorptionProbs_;
        }

        //- Return const access to all residence times
        const scalarList& residenceTimes() const
        {
            return residenceTimes_;
        }

        //- Return const access to the number of stuck parcels
        const scalarList& nStuckParticles() const
        {
            return nStuckParticles_;
        }

        //- Return const access to the number of stuck parcels for a given face
        const scalar& nStuckParticles(const label& facei) const
        {
            return nStuckParticles_[facei];
        }

        //- Return const access to the surface coverage
        const scalarList& surfaceCoverage() const
        {
            return nStuckParticles()/saturationLimit();
        }

        //- Return const access to the surface coverage for a given face
        scalar surfaceCoverage(const label& facei) const
        {
            return nStuckParticles(facei)/saturationLimit(facei);
        }

        //- Return const access to the patch saturation limit
        const scalarList& saturationLimit() const
        {
            return saturationLimit_;
        }

        //- Return const access to the patch saturation limit for a given face
        const scalar& saturationLimit(const label& facei) const
        {
            return saturationLimit_[facei];
        }

        //- Return const access to a molecule adsorption probability
        //  given by typeId
        scalar adsorptionProb(const label& typeId) const
        {
            const label& iD = findIndex(typeIds_, typeId);

            if (iD != -1)
            {
                return adsorptionProbs_[iD];
            }

            FatalErrorIn
            (
                "dsmcStickingDiffuseWallPatch::const scalar&\
                 adsorptionProb(const label& typeId) const"
            )
                << "typeId is not listed."
                << exit(FatalError);

            return 0;
        }

        //- Return const access to a molecule residence time
        //  given by typeId
        scalar residenceTime(const label& typeId) const
        {
            const label& iD = findIndex(typeIds_, typeId);

            if (iD != -1)
            {
                return residenceTimes_[iD];
            }

            FatalErrorIn
            (
                "dsmcStickingWallPatch::const scalar&\
                 residenceTime(const label& typeId) const"
            )
                << "typeId is not listed."
                << exit(FatalError);

            return 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
