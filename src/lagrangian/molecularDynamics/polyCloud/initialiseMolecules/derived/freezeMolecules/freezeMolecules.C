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

Description

\*---------------------------------------------------------------------------*/

#include "freezeMolecules.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(freezeMolecules, 0);

addToRunTimeSelectionTable(polyConfiguration, freezeMolecules, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
freezeMolecules::freezeMolecules
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyConfiguration(molCloud, dict)
{

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        mdInitialiseDict_
    );

    molIds_ = ids.molIds();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

freezeMolecules::~freezeMolecules()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void freezeMolecules::setInitialConfiguration()
{
//     label initialSize = molCloud_.size();

    boundedBox bb;

    setBoundBox(mdInitialiseDict_, bb, "boundBox");

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    label nFrozen = 0;

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            if(bb.contains(mol().position()))
            {
                mol().special() = -2;
                nFrozen++;
            }
        }
    }

    Info << "No of frozen atoms = " << nFrozen << endl;

}

void freezeMolecules::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name
)
{
    const dictionary& dict(propsDict.subDict(name));

    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}




} // End namespace Foam

// ************************************************************************* //
