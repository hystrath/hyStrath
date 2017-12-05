/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "polyMassFluxMultiSurfaces.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMassFluxMultiSurfaces, 0);
addToRunTimeSelectionTable(polyField, polyMassFluxMultiSurfaces, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMassFluxMultiSurfaces::polyMassFluxMultiSurfaces
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    fields_(t, mesh, "dummy"),
    regionIds_(),
    faceZoneNames_(),
    molIds_(),
    fluxDirection_(propsDict_.lookup("fluxDirection")),
    molsZone_(),
    massZone_(),
    timeIndex_(0),
    molFluxZone_(),
    massFluxZone_(),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
{
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput")); 

   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    fluxDirection_ /= mag(fluxDirection_);

    // read in list of zone names - test for one or more defined zones and 
    //                              for multiple definitions

    const List<word> zoneNames (propsDict_.lookup("faceZoneNames"));

    if(!zoneNames.size())
    {
        FatalErrorIn("polyMassFluxMultiSurfaces::polyMassFluxMultiSurfaces()")
            << "Define faceZoneNames " << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    DynamicList<word> regionNames(0);

    forAll(zoneNames, i)
    {
        const word& zoneName(zoneNames[i]);

        if(findIndex(regionNames, zoneName) == -1)
        {
            regionNames.append(zoneName);
        }
        else
        {
            FatalErrorIn("polyMassFluxMultiSurfaces::polyMassFluxMultiSurfaces()")
                << "Zone name: " << zoneName << " cannot be defined twice "
                << nl << "in: "
                << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }
    }

    //regionNames.shrink();

    faceZoneNames_.transfer(regionNames);

    regionIds_.setSize(faceZoneNames_.size(), -1);

    // check facezones are on faceZoneMesh
    const faceZoneMesh& faceZones = mesh_.faceZones();

    forAll(regionIds_, r)
    {
        regionIds_[r] = faceZones.findZoneID(faceZoneNames_[r]);

        if(regionIds_[r] == -1)
        {
            FatalErrorIn("polyMassFluxMultiSurfaces::polyMassFluxMultiSurfaces()")
                << "Cannot find faceZone: " << faceZoneNames_[r] << nl << "in: "
                << time_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }
    }

    molsZone_.setSize(regionIds_.size(), 0.0);
    massZone_.setSize(regionIds_.size(), 0.0);

    molFluxZone_.setSize(regionIds_.size());
    massFluxZone_.setSize(regionIds_.size());

    forAll(molFluxZone_, r)
    {
    	molFluxZone_[r].setSize(1, 0.0);
		massFluxZone_[r].setSize(1, 0.0);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMassFluxMultiSurfaces::~polyMassFluxMultiSurfaces()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMassFluxMultiSurfaces::createField()
{}

//- call this function every time-step before the state and flux objects are cleaned
void polyMassFluxMultiSurfaces::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    const List<scalarField>& molIdFlux = molCloud_.tracker().molIdFlux();
    const List<scalarField>& massIdFlux = molCloud_.tracker().massIdFlux();

    const faceZoneMesh& faceZones = mesh_.faceZones();

    forAll(regionIds_, r)
    {
        scalar molFlux = 0.0;
        scalar massFlux = 0.0;

        const labelList& faces = faceZones[regionIds_[r]];

        forAll(faces, f)
        {
            const label& faceI = faces[f];
            vector nF = mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

            forAll(molIdFlux, id)
            {
                if(findIndex(molIds_, id) != -1)
                {
                    molFlux += (molIdFlux[id][faceI]*nF) & fluxDirection_;
                    massFlux += (massIdFlux[id][faceI]*nF) & fluxDirection_;
                }
            }
        }

        molsZone_[r] += molFlux;
        massZone_[r] += massFlux;
    }

    // -average measurement and calculate properties
    if(time_.outputTime()) 
    {
        scalarField molsZone = molsZone_;
        scalarField massZone = massZone_;

        if(Pstream::parRun())
        {
            forAll(regionIds_, r)
            {
                reduce(molsZone[r], sumOp<scalar>());
                reduce(massZone[r], sumOp<scalar>());
            }
        }

        scalar averagingTime = nAvTimeSteps_*time_.deltaT().value();

        forAll(regionIds_, r)
        {
            molFluxZone_[r][timeIndex_] = molsZone[r]/averagingTime;
            massFluxZone_[r][timeIndex_] = massZone[r]/averagingTime;
        }

        if(resetAtOutput_)
        {
            nAvTimeSteps_ = 0.0;
            molsZone_ = 0.0;
            massZone_ = 0.0;
        }
    }
}

void polyMassFluxMultiSurfaces::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, runTime.timeOutputValue());

            forAll(regionIds_, r)
            {
            	writeTimeData
                (
                    casePath_,
                    "faceFlux_"+faceZoneNames_[r]+fieldName_+"_N.xy",
                    timeField,
                    molFluxZone_[r],
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "faceFlux_"+faceZoneNames_[r]+fieldName_+"_M.xy",
                    timeField,
                    massFluxZone_[r],
                    true
                );

            	writeTimeData
				(
					casePath_,
					"faceFlux_"+faceZoneNames_[r]+fieldName_+"_N.xy",
					molFluxZone_[r],
					true
				);

				writeTimeData
				(
					casePath_,
					"faceFlux_"+faceZoneNames_[r]+fieldName_+"_M.xy",
					massFluxZone_[r],
					true
				);
    
                const reducedUnits& rU = molCloud_.redUnits();
        
                if(rU.outputSIUnits())
                {

                	writeTimeData
                    (
                        casePath_,
                        "faceFlux_"+faceZoneNames_[r]+fieldName_+"_N_SI.xy",
                        timeField*rU.refTime(),
                        molFluxZone_[r]*rU.refMolFlux(),
                        true
                    );
        
                    writeTimeData
                    (
                        casePath_,
                        "faceFlux_"+faceZoneNames_[r]+fieldName_+"_M_SI.xy",
                        timeField*rU.refTime(),
                        massFluxZone_[r]*rU.refMassFlux(),
                        true
                    );

                	writeTimeData
					(
						casePath_,
						"faceFlux_"+faceZoneNames_[r]+fieldName_+"_N_SI.xy",
						molFluxZone_[r]*rU.refMolFlux(),
						true
					);

					writeTimeData
					(
						casePath_,
						"faceFlux_"+faceZoneNames_[r]+fieldName_+"_M_SI.xy",
						massFluxZone_[r]*rU.refMassFlux(),
						true
					);
                }
            }
        }
    }
}

void polyMassFluxMultiSurfaces::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyMassFluxMultiSurfaces::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyMassFluxMultiSurfaces::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
