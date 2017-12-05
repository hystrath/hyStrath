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

#include "polyDelListOfAtoms.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelListOfAtoms, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelListOfAtoms, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelListOfAtoms::polyDelListOfAtoms
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_()
{
    
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    

    molPoints_ = List<vector>(propsDict_.lookup("molPoints"));  
    
    
    findMolsToDel();
  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDelListOfAtoms::~polyDelListOfAtoms()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyDelListOfAtoms::findMolsToDel()
{

    label initialSize = molCloud_.size();
    label deletedMols = 0;
    
    Info << nl << "Deleting the following molecules... " << nl << endl;

    DynamicList<vector> positions;
    DynamicList<scalar> rMinCollect;    
    
    forAll(molPoints_, i)
    {    
        DynamicList<polyMolecule*> molsToDel;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
       
        scalar rMin = GREAT;
        
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                vector rT = molPoints_[i];
                scalar magRIJ = mag(rT-mol().position());
                
                if(magRIJ < rMin)
                {
                    molsToDel.clear();
                    rMin = magRIJ;
                    polyMolecule* molI = &mol();
                    molsToDel.append(molI);
                }
            }
        }
        
       

        forAll(molsToDel, m)
        {
            positions.append(molsToDel[m]->position());
            rMinCollect.append(rMin);     
            
            Info << molsToDel[m]->position()
                << endl;            


                
            deletedMols++;
            deleteMolFromMoleculeCloud(*molsToDel[m]);
        }
    }
    
    Info << nl << "more details ... " << nl << endl;
    
    forAll(molPoints_, i)
    {

        Info << "Deleting molecule at position = "
                << positions[i]
                << ", requested position = " << molPoints_[i]
                << ", residual = " << rMinCollect[i]
                << endl;        
        
    }
    

    label molsKept = initialSize - deletedMols;

    Info<< tab << " initial polyMolecules: " <<  initialSize 
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << deletedMols
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareInteractions();
}

} // End namespace Foam

// ************************************************************************* //
