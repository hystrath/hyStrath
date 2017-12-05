/*---------------------------------------------------------------------------*\
 *  =========                 |
 *  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 *   \\    /   O peration     |
 *    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
 *     \\/     M anipulation  |
 * -------------------------------------------------------------------------------
 * License
 *    This file is part of OpenFOAM.
 * 
 *    OpenFOAM is free software; you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published by the
 *    Free Software Foundation; either version 2 of the License, or (at your
 *    option) any later version.
 * 
 *    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *    for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with OpenFOAM; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * Description
 * 
 * \*---------------------------------------------------------------------------*/

#include "polyZoneRdf.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    
defineTypeNameAndDebug(polyZoneRdf, 0);

addToRunTimeSelectionTable(polyField, polyZoneRdf, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void polyZoneRdf::setRadii()
{
    for(label i = 0; i < nBins_; i++)
    {
        magRadii_[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
        
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyZoneRdf::polyZoneRdf
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
    :
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("rdfName")),
    rMax_(readScalar(propsDict_.lookup("rMax"))),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    binWidth_(rMax_/nBins_),
    molIdA_(-1),
    molIdB_(-1),
    rdf_(nBins_, scalar(0.0)),
    sumRdf_(nBins_, scalar(0.0)),
    magRadii_(nBins_, 0.0),
    binVolume_(nBins_, -1.0),
    totalVolume_(-1.0),
    mols_(0.0),
    nMols_(0.0),
    rdfAveraged_(nBins_, scalar(0.0)),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)    
{
    
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));      

    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);
    
    if(regionId_ == -1)
    {
        FatalErrorIn("polyZoneRdf::polyZoneRdf()")
        << "Cannot find region: " << regionName_ << nl << "in: "
        << time_.time().system()/"fieldPropertiesDict"
        << exit(FatalError);
    }
    
    //-set the total volume
    const labelList& cells = cellZones[regionId_];
    
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        totalVolume_ += mesh_.cellVolumes()[cellI];
    }
    
    if(Pstream::parRun())
    {
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << totalVolume_;
                }
            }
        }
        
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalar totalVolumeProc;
                
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> totalVolumeProc;
                }
                
                totalVolume_ += totalVolumeProc;
            }
        }
    }
    
    forAll(binVolume_, n)
    {
        binVolume_[n] = (4.0/3.0)*constant::mathematical::pi*(pow((n+1.0),3.0)-pow(n,3.0))*pow(binWidth_,3.0);
    }
    
    setRadii();

    const List<word>& idList(molCloud_.cP().molIds());
    const word molIdA = propsDict_.lookup("molIdA");
    molIdA_ = findIndex(idList, molIdA);
    const word molIdB = propsDict_.lookup("molIdB");
    molIdB_ = findIndex(idList, molIdB);
    
    if(molIdA_ == -1)
    {
        FatalErrorIn("polyZoneRdf::polyZoneRdf()")
        << "Cannot find molIdA: " << molIdA << nl << "in: "
        << time_.time().system()/"fieldPropertiesDict"
        << exit(FatalError);
    }
    
    if(molIdB_ == -1)
    {
        FatalErrorIn("polyZoneRdf::polyZoneRdf()")
        << "Cannot find molIdB: " << molIdB << nl << "in: "
        << time_.time().system()/"fieldPropertiesDict"
        << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyZoneRdf::~polyZoneRdf()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyZoneRdf::createField()
{
    
}

void polyZoneRdf::evaluatePair
(
    polyMolecule* molI,
    polyMolecule* molJ
)

{
    label idI = molI->id();
    label idJ = molJ->id();
    
    vector rsIsJ = molI->position() - molJ->position();
    
    if(idI == molIdA_ && idJ == molIdB_)
    {
        scalar rsIsJMag = mag(rsIsJ);
        
        if(rsIsJMag < rMax_)
        {
            label ig = label(rsIsJMag/binWidth_);
            rdf_[ig] += 2.0;
        }
        
    }
}


void polyZoneRdf::evaluatePairS
(
    polyMolecule* molReal,
    polyMolecule* molRef
)

{
    label idReal = molReal->id();
    label idRef = molRef->id();
    
    vector rsRealsRef =molReal->position() - molRef->position();
    
    if(idReal == molIdA_ && idRef == molIdB_)
    {
        scalar rsRealsRefMag = mag(rsRealsRef);
        
        if(rsRealsRefMag < rMax_)
        {
            label ig = label(rsRealsRefMag/binWidth_);
            rdf_[ig] ++;
        }
    }
}

void polyZoneRdf::calculateField()
{
    
    nAvTimeSteps_ += 1.0;
    
    // - sampling density within the spherical zone 
    const List< DynamicList<polyMolecule*> >& cellOccupancy = molCloud_.cellOccupancy();
    
    createField();
    mols_ = 0;
    
    const labelList& cells = mesh_.cellZones()[regionId_];
    
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];
        
        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];
            
            //As far as the two species have the same number of particles this works fine,
            // otherwise I have to do as in molecular where we count both individually 
            // and multiply at the end the 2 number molsA*molsB 
            if(molI->id() == molIdA_)
            {
                mols_ += 1.0; 
            }
        }
    }
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    {
        // Real-Real rdf
        
        polyMolecule* molI = &mol();
        polyMolecule* molJ = &mol();
        
        const labelListList& dil=molCloud_.il().dil();
        
        molCloud_.il().setRIPL();
        
        const List<DynamicList<polyMolecule*> >& ril=molCloud_.il().ripl();
            
        forAll(dil, d)
        {   
            bool flag = false;
            
            forAll(cells, c)
            {
                if(cells[c] == d)
                {
                    flag = true;
                    break;
                }
            }
            
            if(flag)
            {
                forAll(cellOccupancy[d],cellIMols)
                {
                    molI = cellOccupancy[d][cellIMols];
                    
                    forAll(dil[d], interactingCells)
                    {
                        List< polyMolecule* > cellJ =
                        cellOccupancy[dil[d][interactingCells]];
                        
                        forAll(cellJ, cellJMols)
                        {
                            
                            molJ = cellJ[cellJMols];
                            
                            evaluatePair(molI, molJ);
                        }
                    }
                    
                    forAll(cellOccupancy[d],cellIOtherMols)
                    {
                        molJ = cellOccupancy[d][cellIOtherMols];
                        
                        if (molJ > molI)
                        {
                            evaluatePair(molI, molJ);
                        }
                    }
                    
                    forAll(ril[d], r)
                    {
                        molJ= ril[d][r];
                        evaluatePairS(molI, molJ);
                    }
                }
            }
        }
    }
    
    
    
    if (Pstream::parRun())
    {
        reduce(mols_, sumOp<scalar>());
        reduce(rdf_, sumOp<scalarField>());
    }
    
    nMols_ += mols_;
    sumRdf_ += rdf_;
    rdf_ = scalar(0.0);
    

    
    if(time_.outputTime())
    {
        scalar nAvTimeSteps = nAvTimeSteps_;
        scalar density = nMols_/(totalVolume_*nAvTimeSteps);
        
        forAll(sumRdf_, ic)
        {
            rdfAveraged_[ic] = sumRdf_[ic]/(density*nMols_*binVolume_[ic]);
        }    
        
        if(resetAtOutput_)
        {
            nAvTimeSteps_ = 0.0;
            nMols_ = 0.0;
            sumRdf_ = scalar(0.0);
//                Info << "Reseting field !!!" << endl; 
        }
    }
}

void polyZoneRdf::writeField()
{
    if(time_.outputTime())
    {
        if(Pstream::master())
        {
            fileName timePath(time_.path()/time_.timeName()/"uniform");
            
            writeTimeData(timePath, fieldName_+"_rdfBin", magRadii_, rdfAveraged_);
            
            const reducedUnits& rU = molCloud_.redUnits();
            
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    timePath,
                    fieldName_+"_rdfBin_SI",
                    magRadii_*rU.refLength(),
                    rdfAveraged_
                );
            }
        }
    }
}


void polyZoneRdf::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void polyZoneRdf::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

const propertyField& polyZoneRdf::fields() const
{
    return  fields_;
}

/*
void polyZoneRdf::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);
    
}
   */ 
} // End namespace Foam

// ************************************************************************* //
