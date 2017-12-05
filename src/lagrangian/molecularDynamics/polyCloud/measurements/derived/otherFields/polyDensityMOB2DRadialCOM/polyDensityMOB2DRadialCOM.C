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

#include "polyDensityMOB2DRadialCOM.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDensityMOB2DRadialCOM, 0);

addToRunTimeSelectionTable(polyField, polyDensityMOB2DRadialCOM, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void polyDensityMOB2DRadialCOM::setRadii()
// {
//     for(label i = 0; i < nBinsX_; i++)
//     {
//         for(label j = 0; j < nBinsY_; j++)
//         {
//             magRadii_[i][j] = 0.5*binWidthX_ + scalar(i)*binWidthX_;
//         }
// //         radii_[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDensityMOB2DRadialCOM::polyDensityMOB2DRadialCOM
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
    fieldName_(propsDict_.lookup("fieldName")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    molIds_(),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    h_(mag(endPoint_ - startPoint_)),
    unitVectorX_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    radius_(readScalar(propsDict_.lookup("radius"))),
    nBinsX_(readLabel(propsDict_.lookup("nBinsX"))),
    nBinsY_(readLabel(propsDict_.lookup("nBinsY"))),
    binWidthX_(h_/(nBinsX_)),
    binWidthY_(radius_/nBinsY_),
    mols_(nBinsX_),
    oldCenterOfMass_(propsDict_.lookup("startPoint")),

    densityField_(nBinsX_),

    magRadii_(),
    binWidths_(),
    volume_(),
    avVolume_(0.0),

    minBinWidth_(0.0),
    n_(),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
{
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));  
    
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    { 
        FatalErrorIn("polyDensityRadialMOBZone::polyDensityRadialMOBZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    avVolume_ = radius_*radius_*constant::mathematical::pi*binWidthX_/nBinsY_;

    DynamicList<scalar> radii(0);
    DynamicList<scalar> binWidths(0);

    binWidths.append(sqrt(avVolume_/(constant::mathematical::pi*binWidthX_)));

    radii.append(0.5*binWidths[0]);

    for(label i = 1; i < nBinsY_; i++)
    {
        scalar r2 = radii[i-1] + 0.5*binWidths[i-1];
        scalar r1 = sqrt((avVolume_/(constant::mathematical::pi*binWidthX_)) + (r2*r2));

        if(r2 <= radius_)
        {
            radii.append(0.5*(r1+r2));
            binWidths.append(r1-r2);
        }
        else
        {
            break;
        }
        
    }
    
    binWidths.shrink();
    radii.shrink(); 
    
    nBinsY_ = binWidths.size();

    magRadii_.setSize(nBinsY_, 0.0);
    volume_.setSize(nBinsY_, 0.0);
    binWidths_.setSize(nBinsY_, 0.0);

    forAll(mols_, x)
    {
        mols_[x].setSize(nBinsY_, 0.0);
        densityField_[x].setSize(nBinsY_, 0.0);
    }

    forAll(magRadii_, n)
    {
        magRadii_[n] = radii[n];
        binWidths_[n] = binWidths[n];
        
        if(n == 0)
        {
            volume_[n] = constant::mathematical::pi*binWidths_[n]*binWidths_[n]*binWidthX_;
        }
        else
        {
            volume_[n] = 2.0*constant::mathematical::pi*magRadii_[n]*binWidths_[n]*binWidthX_;
        }
    }

    minBinWidth_ = binWidths_[nBinsY_-1]/3.0;

    n_.setSize(label(radius_/minBinWidth_), 0);

    forAll(n_, n)
    {
        scalar r = 0.5*minBinWidth_ + scalar(n)*minBinWidth_;

        for(label i = 0; i < nBinsY_; i++)
        {
            scalar rLimit1 = magRadii_[i] - 0.5*binWidths_[i];
            scalar rLimit2 = magRadii_[i] + 0.5*binWidths_[i];

            if((r >= rLimit1) && (r < rLimit2))
            {
                n_[n] = i;
            }
        }
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDensityMOB2DRadialCOM::~polyDensityMOB2DRadialCOM()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label polyDensityMOB2DRadialCOM::findBin(const scalar& r)
{
    label n = -1;

    label n2 = label(r/minBinWidth_) - 1;

    if(n2 < 0)
    {
        n2 = 0;
    }

    if(n2 < n_.size())
    {
        n = n_[n2];

        scalar rLimit1 = magRadii_[n] - 0.5*binWidths_[n];
        scalar rLimit2 = magRadii_[n] + 0.5*binWidths_[n];

        if((r >= rLimit1) && (r < rLimit2))
        {}
        else
        {
            n++;
        }
    }
    
    return n;
}




void polyDensityMOB2DRadialCOM::createField()
{}


void polyDensityMOB2DRadialCOM::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    const boundBox& globalBb = mesh_.bounds();
    vector domainLength = globalBb.max() - globalBb.min();
    scalar domainLengthX = fabs(domainLength.x());
    scalar domainLengthY = fabs(domainLength.y());
    scalar domainLengthZ = fabs(domainLength.z());
        

    {
        const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];
        scalar mass = 0.0;
        vector centreOfMass = vector::zero;

        forAll(cells, c)
        {
            const label& cellI = cells[c];
            const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

            if(molsInCell.size() > 0)
            {
                forAll(molsInCell, mIC)
                {                    
                    polyMolecule* molI = molsInCell[mIC];

                    vector rI = molI->position(); 
                    vector roC = rI - oldCenterOfMass_;
                    
                    if(fabs(roC.x())>domainLengthX/2.0)
                    {
                        if(rI.x() > oldCenterOfMass_.x())
                        {    
                            rI.x() -= domainLengthX;
                        }
                        else
                        {    
                            rI.x() += domainLengthX;
                        }
                    }
                    
                    if(fabs(roC.y())>domainLengthY/2.0)
                    {                        
                        if(rI.y() > oldCenterOfMass_.y())
                        {    
                            rI.y() -= domainLengthY;
                        }
                        else
                        {    
                            rI.y() += domainLengthY;
                        }   
                    }
                    
                    if(fabs(roC.z())>domainLengthZ/2.0)
                    {                        
                        if(rI.z() > oldCenterOfMass_.z())
                        {    
                            rI.z() -= domainLengthZ;
                        }
                        else
                        {    
                            rI.z() += domainLengthZ;
                        }   
                    }
                    
                    if(findIndex(molIds_, molI->id()) != -1)
                    {
                        const scalar& massI = molCloud_.cP().mass(molI->id());
                        mass += massI;
                        centreOfMass += massI*rI;
                    }
                }
            }
        }

        if(Pstream::parRun())
        {
            //-sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << mass << centreOfMass;
                    }
                }
            }

            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar massProc;
                    vector centreOfMassProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> massProc >> centreOfMassProc;
                    }
                    mass += massProc;
                    centreOfMass += centreOfMassProc;
                }
            }
        }

        centreOfMass /= mass;

        centreOfMass =  ((centreOfMass ^ unitVectorX_)^ unitVectorX_)*(-1.0)
                        + (startPoint_ & unitVectorX_) * unitVectorX_;
           
        if(centreOfMass.x() > globalBb.max().x())
        {
            centreOfMass.x() -= domainLengthX;
        }
        else if(centreOfMass.x() < globalBb.min().x())
        {
            centreOfMass.x() += domainLengthX;
        }
        if(centreOfMass.y() > globalBb.max().y())
        {
            centreOfMass.y() -= domainLengthY;
        }
        else if(centreOfMass.y() < globalBb.min().y())
        {
            centreOfMass.y() += domainLengthY;
        }
        if(centreOfMass.z() > globalBb.max().z())
        {
            centreOfMass.z() -= domainLengthZ;
        }
        else if(centreOfMass.z() < globalBb.min().z())
        {
            centreOfMass.z() += domainLengthZ;
        }
        
         oldCenterOfMass_ = centreOfMass;
            
        forAll(cells, c)
        {
            const label& cellI = cells[c];
            const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

            if(molsInCell.size() > 0)
            {
                forAll(molsInCell, mIC)
                {
                    polyMolecule* molI = molsInCell[mIC];

                    vector rI = molI->position(); 

                    vector rSI = rI - centreOfMass;
                    
                    if(fabs(rSI.x())>domainLengthX/2.0)
                    {   
                        if(rI.x() > centreOfMass.x())
                        {    
                            rI.x() -= domainLengthX;
                        }
                        else
                        {    
                            rI.x() += domainLengthX;
                        }
                    }
                    
                    if(fabs(rSI.y())>domainLengthY/2.0)
                    {
                        if(rI.y() > centreOfMass.y())
                        {    
                            rI.y() -= domainLengthY;
                        }
                        else
                        {    
                            rI.y() += domainLengthY;
                        }
                    }
                    
                    if(fabs(rSI.z())>domainLengthZ/2.0)
                    { 
                        if(rI.z() > centreOfMass.z())
                        {    
                            rI.z() -= domainLengthZ;
                        }
                        else
                        {    
                            rI.z() += domainLengthZ;
                        }
                    }
                    
                    scalar rD = rSI & unitVectorX_;

                    if((rD <= h_) && (rD >= 0.0))
                    {   
                        scalar rN = mag((rD*unitVectorX_ + centreOfMass) - rI);
                        
                        if(rN <= radius_)
                        {
                            label nY = findBin(rN);
                            // nY defines the No. in the radial direction.
                            if
                            (
                                nY != -1
                            )
                            {
                                label nX = label(rD/binWidthX_);
                                // nX defines the No. in the axial direction.

                                if(nX == nBinsX_)
                                {
                                    nX--;
                                }

                                if(findIndex(molIds_, molI->id()) != -1)
                                {
                                    mols_[nX][nY] += 1.0;                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(time_.outputTime())
    {
        Field<scalarField> mols = mols_;

        //- parallel communication
        if(Pstream::parRun())
        {
            //-sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << mols;
                    }
                }
            }

            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    Field<scalarField> molsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> molsProc;
                    }

                    forAll(molsProc, x)
                    {
                        mols[x] += molsProc[x];
                    }
                }
            }
        }
        
//         const scalar& nAvTimeSteps = time_.nAvTimeSteps().value();



        forAll(densityField_, x)
        {
            forAll(densityField_[x], y)
            {
                densityField_[x][y] = mols[x][y]/(nAvTimeSteps_*volume_[y]);
            }
        }
        
        //- reset fields
        if(resetAtOutput_)
        {
            nAvTimeSteps_ = 0.0;
            
            forAll(mols_, x)
            {
                mols_[x] = 0.0;
            }
        }
    }
}

void polyDensityMOB2DRadialCOM::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            OFstream file(timePath_/fieldName_+"_"+regionName_+"_density2DRadialCOM.xyz");

            if(file.good())
            {
                forAll(densityField_, x)
                {
                    scalar rX = 0.5*binWidthX_ + scalar(x)*binWidthX_;

                    forAll(densityField_[x], y)
                    {
                        scalar rY = magRadii_[y];

                        file
                            << rY << "\t"
                            << rX << "\t"
                            << densityField_[x][y]
                            << endl;
                    }

                    file << endl;
                }
            }
            else
            {
                FatalErrorIn("void polyDensityMOB2DRadialCOM::writeField()")
                    << "Cannot open file " << file.name()
                    << abort(FatalError);
            }
        }
    }
}

void polyDensityMOB2DRadialCOM::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

// void polyDensityMOB2DRadialCOM::measureDuringForceComputation
// (
//     polyMolecule* molReal,
//     polyReferredMolecule* molRef
// ){}

void polyDensityMOB2DRadialCOM::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

// void polyDensityMOB2DRadialCOM::measureDuringForceComputationSite
// (
//     polyMolecule* molReal,
//     polyReferredMolecule* molRef,
//     label sReal,
//     label sRef
// ){}

const propertyField& polyDensityMOB2DRadialCOM::fields() const
{
    return  fields_;
}

// void polyDensityMOB2DRadialCOM::updateProperties(const dictionary& newDict)
// {
//     //- the main properties should be updated first
//     updateBasicFieldProperties(newDict);
// 
// }

} // End namespace Foam

// ************************************************************************* //
