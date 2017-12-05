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
    for more deciliumMotions.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "ciliumMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ciliumMotion, 0);

addToRunTimeSelectionTable(polyStateController, ciliumMotion, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ciliumMotion::ciliumMotion
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),    
    molIds_()
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();


    nC_ = readLabel(propsDict_.lookup("numberOfCiliumPoints"));
    nT_ = readLabel(propsDict_.lookup("numberOfTimeSteps"));      

    const reducedUnits& rU = molCloud_.redUnits();
    
    x_.setSize(nT_);
    y_.setSize(nT_);
//     t_.setSize(nT_);
    
    forAll(x_, i)
    {
        x_[i].setSize(nC_);
    }

    forAll(y_, i)
    {
        y_[i].setSize(nC_);
    }
    
    {
        ifstream file("X.xy");
        
        if(file.is_open())
        {
            for (label i = 0; i < nT_; i++)
            {
                for (label j = 0; j < nC_; j++)
                {
                    file >> x_[i][j];
                    x_[i][j] /= rU.refLength();
                }
            }
        }
    }
    
//     Pout << "X_" << X_[1][0] << endl;
    
//     Info << "x = " << x_ << endl;
  
    {
        ifstream file("Y.xy");
        
        if(file.is_open())
        {
            for (label i = 0; i < nT_; i++)
            {
                for (label j = 0; j < nC_; j++)
                {
                    file >> y_[i][j];
                    y_[i][j] /= rU.refLength();
                }
            }
        }
    }
    

   // read in tracking numbers
    
    trackingNumbers_ = List<label>(propsDict_.lookup("trackingNumbers"));
    startPoint_ = propsDict_.lookup("startPoint");
    
    deltaTMD_ = time_.deltaT().value(); 
    
    tI_ = 0;
    
    // allow start from latest time
    if(time_.time().timeOutputValue() > 0.0)
    {
        scalar nTimeSteps = time_.time().timeOutputValue()/deltaTMD_;
        label factor = label(nTimeSteps/scalar(nT_));
        tI_ = label(nTimeSteps-factor*nT_);
        Info << "starting at time index = " << tI_ << endl;
    }
    
    
    nWrite_ = 0;
    nAvTimeSteps_ = 0.0;
    velocities_.setSize(trackingNumbers_.size(), vector::zero);
    forces_.setSize(trackingNumbers_.size(), vector::zero);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ciliumMotion::~ciliumMotion()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ciliumMotion::initialConfiguration()
{
}
    




void ciliumMotion::controlBeforeVelocityI()
{
    
}

void ciliumMotion::controlBeforeMove()
{
    tI_++;

    if(tI_ == nT_)
    {
        tI_ = 0;
    }
    
    nAvTimeSteps_ += 1.0;
    
//     Info << "ciliumMotion: control, t = " << tI_ << endl;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label id = findIndex(trackingNumbers_, mol().trackingNumber());
            
            if(id != -1)
            {
                vector rNew = startPoint_ + vector(x_[tI_][id], y_[tI_][id], 0.0);
                
                vector deltaR = rNew - mol().position();
                
                mol().a() = vector::zero;
                mol().v() = deltaR/deltaTMD_;
                
                velocities_[id] += mol().v();
            }
        }
    }    
}


void ciliumMotion::controlBeforeForces()
{}

void ciliumMotion::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void ciliumMotion::controlAfterForces()
{
    
}

void ciliumMotion::controlAfterVelocityII()
{}

void ciliumMotion::calculateProperties()
{
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label id = findIndex(trackingNumbers_, mol().trackingNumber());
            
            if(id != -1)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());                
                forces_[id] += massI*mol().a();
            }
        }
    }
}

void ciliumMotion::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    if(time_.outputTime())
    {
        nWrite_++;
        
        if(Pstream::parRun())
        {
            forAll(forces_, i)
            {
                reduce(forces_[i], sumOp<vector>());
                reduce(velocities_[i], sumOp<vector>());
            }
        }
        
        forAll(forces_, i)
        {
            forces_[i] /= nAvTimeSteps_;
            velocities_[i] /= nAvTimeSteps_;
        }

        
        if(Pstream::master())
        {
            List<vectorField> forces(1);
            List<vectorField> velocities(1);
            
            forces[0].setSize(forces_.size());
            velocities[0].setSize(velocities_.size());
            
            forAll(forces_, i)
            {
                forces[0][i]=forces_[i];
                velocities[0][i]=velocities_[i];
            }


            

            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_Fx.xy",
                forces,
                "x",
                true
            );
            
            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_Fy.xy",
                forces,
                "y",
                true
            );  
            
            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_Fz.xy",
                forces,
                "z",
                true
            );                

            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_Ux.xy",
                velocities,
                "x",
                true
            );
           
            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_Uy.xy",
                velocities,
                "y",
                true
            );
            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_Uz.xy",
                velocities,
                "z",
                true
            );            
          
//             std::string s;
//             std::stringstream out;
//             out << nWrite_;
//             s = out.str();
            
//                 writeTimeData
//                 (
//                     timePath_,
//                     "bins_twoDim_"+fieldName_+"_"+s+"_rhoM.xy",
//                     binsY,
//                     rhoM_[i]
//                 );        
   
/*            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_"+ s +"_forces.xyz",
                forces_
            );            
            writeTimeData
            (
                fixedPathName,
                "ciliumMotion_"+fieldName_+"_"+ s +"_velocities.xyz",
                velocities_
            ); */           
        }
        
        // reset
        nAvTimeSteps_ = 0.0;
        forces_ = vector::zero;
        velocities_ = vector::zero;
    }
}

void ciliumMotion::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");


}

} // End namespace Foam

// ************************************************************************* //
