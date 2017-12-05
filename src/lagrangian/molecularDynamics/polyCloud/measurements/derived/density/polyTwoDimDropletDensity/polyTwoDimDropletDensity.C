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

#include "polyTwoDimDropletDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyTwoDimDropletDensity, 0);

addToRunTimeSelectionTable(polyField, polyTwoDimDropletDensity, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyTwoDimDropletDensity::polyTwoDimDropletDensity
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
//     regionName_(propsDict_.lookup("zoneName")),
//     regionId_(-1),
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

    oldCenterOfMass_(propsDict_.lookup("startPoint")), //****

    magRadii_(),
    binWidths_(),
    volume_(),
    avVolume_(0.0),

    minBinWidth_(0.0),
    n_()

{
    {
        dropletMolIds_.clear();

        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "dropletMolIds"
        );

        dropletMolIds_ = ids.molIds();
    }
    
    {
        molIds_.clear();

        selectIds ids
        (
            molCloud_.cP(),
            propsDict_
        );

        molIds_ = ids.molIds();
    }
    
    // bins 
    
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

    mols_.setSize(nBinsX_);
    mass_.setSize(nBinsX_);    
    rhoN_.setSize(nBinsX_);
    rhoM_.setSize(nBinsX_);
    
    magRadii_.setSize(nBinsY_, 0.0);
    volume_.setSize(nBinsY_, 0.0);
    binWidths_.setSize(nBinsY_, 0.0);

    forAll(mols_, x)
    {
        mols_[x].setSize(nBinsY_, 0.0);
        mass_[x].setSize(nBinsY_, 0.0);
        rhoN_[x].setSize(nBinsY_, 0.0);
        rhoN_[x].setSize(nBinsY_, 0.0);
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
    
    const boundBox& globalBb = mesh_.bounds();
    
    box_.resetBoundedBox(globalBb.max(), globalBb.min());
    
    x_ = box_.span().x();
    y_ = box_.span().y();
    z_ = box_.span().z();
    
    nAvTimeSteps_ = 0.0;
    
    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));     
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyTwoDimDropletDensity::~polyTwoDimDropletDensity()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label polyTwoDimDropletDensity::findBin(const scalar& r)
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




void polyTwoDimDropletDensity::createField()
{}


void polyTwoDimDropletDensity::calculateField()
{
    nAvTimeSteps_ += 1.0; 

    scalar mass = 0.0;
    vector centreOfMass = vector::zero;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(dropletMolIds_, mol().id()) != -1)            
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());

                vector rI = mol().position(); 
                vector roC = rI - oldCenterOfMass_;
                
                if(fabs(roC.x())>x_/2.0)
                {
                    if(rI.x() > oldCenterOfMass_.x())
                    {    
                        rI.x() -= x_;
                    }
                    else
                    {    
                        rI.x() += x_;
                    }
                }
                
                if(fabs(roC.y())>y_/2.0)
                {                        
                    if(rI.y() > oldCenterOfMass_.y())
                    {    
                        rI.y() -= y_;
                    }
                    else
                    {    
                        rI.y() += y_;
                    }   
                }
                
                if(fabs(roC.z())>z_/2.0)
                {                        
                    if(rI.z() > oldCenterOfMass_.z())
                    {    
                        rI.z() -= z_;
                    }
                    else
                    {    
                        rI.z() += z_;
                    }   
                }

                mass += massI;
                centreOfMass += massI*rI;
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(mass, sumOp<scalar>());
        reduce(centreOfMass, sumOp<vector>());
    }

    centreOfMass /= mass;

    centreOfMass =  ((centreOfMass ^ unitVectorX_)^ unitVectorX_)*(-1.0)
                    + (startPoint_ & unitVectorX_) * unitVectorX_;
        
    if(centreOfMass.x() > box_.max().x())
    {
        centreOfMass.x() -= x_;
    }
    else if(centreOfMass.x() < box_.min().x())
    {
        centreOfMass.x() += x_;
    }
    if(centreOfMass.y() > box_.max().y())
    {
        centreOfMass.y() -= y_;
    }
    else if(centreOfMass.y() < box_.min().y())
    {
        centreOfMass.y() += y_;
    }
    if(centreOfMass.z() > box_.max().z())
    {
        centreOfMass.z() -= z_;
    }
    else if(centreOfMass.z() < box_.min().z())
    {
        centreOfMass.z() += z_;
    }
    
    oldCenterOfMass_ = centreOfMass;
        
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(dropletMolIds_, mol().id()) != -1)            
            {
                vector rI = mol().position(); 

                vector rSI = rI - centreOfMass;
                
                if(fabs(rSI.x())>x_/2.0)
                {   
                    if(rI.x() > centreOfMass.x())
                    {    
                        rI.x() -= x_;
                    }
                    else
                    {    
                        rI.x() += x_;
                    }
                }
                
                if(fabs(rSI.y())>y_/2.0)
                {
                    if(rI.y() > centreOfMass.y())
                    {    
                        rI.y() -= y_;
                    }
                    else
                    {    
                        rI.y() += y_;
                    }
                }
                
                if(fabs(rSI.z())>z_/2.0)
                { 
                    if(rI.z() > centreOfMass.z())
                    {    
                        rI.z() -= z_;
                    }
                    else
                    {    
                        rI.z() += z_;
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

                            if(findIndex(molIds_, mol().id()) != -1)
                            {
                                const scalar& massI = molCloud_.cP().mass(mol().id());
                                
                                mols_[nX][nY] += 1.0;
                                mass_[nX][nY] += massI;
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
        Field<scalarField> mass = mass_;
        
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
                        toNeighbour << mols << mass;
                    }
                }
            }

            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    Field<scalarField> molsProc;
                    Field<scalarField> massProc;
    
                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> molsProc >> massProc;
                    }

                    forAll(molsProc, x)
                    {
                        mols[x] += molsProc[x];
                        mass[x] += massProc[x];
                    }
                }
            }
        }
        
        

        forAll(mols, x)
        {
            forAll(mols[x], y)
            {
                rhoN_[x][y] = mols[x][y]/(nAvTimeSteps_*volume_[y]);
                rhoM_[x][y] = mass[x][y]/(nAvTimeSteps_*volume_[y]);
            }
        }
        
        //- reset fields
        if(resetAtOutput_)
        {
            forAll(mols_, x)
            {
                mass_[x] = 0.0;
                mols_[x] = 0.0;
            }
        }
    }
}

void polyTwoDimDropletDensity::writeField()
{
    if(time_.outputTime())
    {
        if(Pstream::master())
        {
            {
                OFstream file(timePath_/"twoDim_"+fieldName_+"_rhoN.xyz");

                if(file.good())
                {
                    forAll(rhoN_, x)
                    {
                        scalar rX = 0.5*binWidthX_ + scalar(x)*binWidthX_;

                        forAll(rhoN_[x], y)
                        {
                            scalar rY = magRadii_[y];

                            file
                                << rY << "\t"
                                << rX << "\t"
                                << rhoN_[x][y]
                                << endl;
                        }

                        file << endl;
                    }
                }
                else
                {
                    FatalErrorIn("void polyTwoDimDropletDensity::writeField()")
                        << "Cannot open file " << file.name()
                        << abort(FatalError);
                }
            }
            
            {
                OFstream file(timePath_/"twoDim_"+fieldName_+"_rhoN.xyz");

                if(file.good())
                {
                    forAll(rhoM_, x)
                    {
                        scalar rX = 0.5*binWidthX_ + scalar(x)*binWidthX_;

                        forAll(rhoM_[x], y)
                        {
                            scalar rY = magRadii_[y];

                            file
                                << rY << "\t"
                                << rX << "\t"
                                << rhoM_[x][y]
                                << endl;
                        }

                        file << endl;
                    }
                }
                else
                {
                    FatalErrorIn("void polyTwoDimDropletDensity::writeField()")
                        << "Cannot open file " << file.name()
                        << abort(FatalError);
                }
            }    
        }
    }
}

void polyTwoDimDropletDensity::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}


void polyTwoDimDropletDensity::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& polyTwoDimDropletDensity::fields() const
{
    return  fields_;
}

// void polyTwoDimDropletDensity::updateProperties(const dictionary& newDict)
// {
//     //- the main properties should be updated first
//     updateBasicFieldProperties(newDict);
// 
// }

} // End namespace Foam

// ************************************************************************* //
