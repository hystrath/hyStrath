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

#include "polyMappingZone.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMappingZone, 0);

addToRunTimeSelectionTable(polyMappingModel, polyMappingZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMappingZone::polyMappingZone
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMappingModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    translation_(propsDict_.lookup("translationalVector")),
    molIds_()
{
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    findMolsToMap();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMappingZone::~polyMappingZone()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyMappingZone::findMolsToMap()
{
    DynamicList<polyMolecule*> molsToDelete;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                vector newPos = mol().position() + translation_;

                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    newPos,
                    cell,
                    tetFace,
                    tetPt
                );
                
                if(cell == -1)
                {
                    polyMolecule* molI = &mol();
                    molsToDelete.append(molI);
                }
                
                mol().position() = newPos;
                mol().cell() = cell;
                mol().tetFace() = tetFace;
                mol().tetPt() = tetPt;

                // shift the sites too
                forAll( mol().sitePositions(), i )
                {
                    mol().sitePositions()[i] += translation_;
                }
            }
        }
    }

    //molsToDelete.shrink();

    Info<< " deleting " << molsToDelete.size() << " molecules "
        << endl;

    forAll (molsToDelete, mTD)
    {
        molCloud_.deleteParticle(*(molsToDelete[mTD]));
    }
}

} // End namespace Foam

// ************************************************************************* //
