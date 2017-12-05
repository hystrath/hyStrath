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

#include "pairPotentialModel.H"
#include "energyScalingFunction.H"
#include "IOstreams.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pairPotentialModel, 0);

defineRunTimeSelectionTable(pairPotentialModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pairPotentialModel::pairPotentialModel
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud, 
    const reducedUnits& redUnits,
    const word& name,
    const dictionary& dict
)
:
    mesh_(mesh),
    molCloud_(molCloud),
    rU_(redUnits),
    pairPotentialProperties_(dict),    
    name_(name),
    idList_(),
    rCut_(readScalar(dict.lookup("rCut"))),
    rMin_(readScalar(dict.lookup("rMin"))),
    dr_(readScalar(dict.lookup("dr"))),
    useTables_(true),
    forceLookup_(0),
    energyLookup_(0),
    esfPtr_(NULL),
    writeTables_(false),
    exclusions_(false)   
{
    writeTables_ = false;  

    if (dict.found("writeTables"))
    {
        writeTables_ = Switch(dict.lookup("writeTables"));
    }
    
    if(rU_.runReducedUnits())
    {
        rCut_ /= rU_.refLength();
        rMin_ /= rU_.refLength();
        dr_ /= rU_.refLength();
        rCutSqr_ = rCut_*rCut_;
    }
    
    // splitting the name using a delimeter "A-B" => "A" and "B"
    idList_.setSize(2);
    
//     Info << nl << "name = " << name_ << endl;
    
    std::string s = name_;
    std::string delimiter = "-";
    
    size_t pos = 0;
    std::string token;

    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        idList_[0]=token;
        s.erase(0, pos + delimiter.length());
        idList_[1]=s;
    }
    
//     Info << " idList = " << idList_ << endl;
    
    
    // exclusions 
    if(dict.found("exclusionModel"))
    {
        exclusionModel_ = autoPtr<exclusionModel>
        (
            exclusionModel::New(mesh, molCloud_, dict)
        );
        
        if(exclusionModel_->type() != "noExclusions")
        {
            exclusions_ = true;
        }
    }
    
/*    Info << "pairPotentialModel, " << name_ <<" exclusionModel =  " 
    << exclusions_ << endl;  */  
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pairPotentialModel> pairPotentialModel::New
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& redUnits, 
    const word& name, 
    const dictionary& dict
)
{
    word pairPotentialModelName
    (
        dict.lookup("pairPotential")
    );

    Info<< "Selecting model: "
         << pairPotentialModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pairPotentialModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pairPotentialModel::New(const dictionary&) : " << endl
            << "    unknown pairPotential type "
            << pairPotentialModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pairPotentialModel>
    (
        cstrIter()(mesh, molCloud, redUnits, name, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pairPotentialModel::~pairPotentialModel()
{}

void Foam::pairPotentialModel::scaleEnergy
(
    scalar& e,
    const scalar r
) const
{
    if (!esfPtr_)
    {
        esfPtr_ = energyScalingFunction::New
        (
            name_, pairPotentialProperties_, *this, rU_
        ).ptr();
    }

    esfPtr_->scaleEnergy(e, r);
}

void pairPotentialModel::setLookupTables()
{
    label N = label((rCut_ - rMin_)/dr_) + 1;

    Info << "Number of entries: " << N << endl;

    forceLookup_.setSize(N);
    energyLookup_.setSize(N);

    forAll(forceLookup_, k)
    {
        energyLookup_[k] = scaledEnergy(k*dr_ + rMin_);
        forceLookup_[k] = -energyDerivative((k*dr_ + rMin_), true);
    }
    
    fMin_ = forceLookup_[0];
    energyMin_ = energyLookup_[0];
}

scalar pairPotentialModel::forceLookUpFromTable(const scalar r) const
{
    if(r < rMin_)
    {
        return fMin_;
    }
    else
    {
        scalar k_rIJ = (r - rMin_)/dr_;

        label k = label(k_rIJ);

        scalar f = 0.0;

        if(k < forceLookup_.size()-1)
        {
            f =
                (k_rIJ - k)*forceLookup_[k+1]
            + (k + 1 - k_rIJ)*forceLookup_[k];
        }
    
        return f;
    }
}

List< Pair< scalar > >
pairPotentialModel::forceTable() const
{
    List<Pair<scalar> > forceTab(forceLookup_.size());

    forAll(forceLookup_,k)
    {
        forceTab[k].first() = rMin_ + k*dr_;

        forceTab[k].second() = forceLookup_[k];
    }

    return forceTab;
}

scalar pairPotentialModel::energyLookUpFromTable(const scalar r) const
{
    if(r < rMin_)
    {
        return energyMin_;
    }
    else
    {
        scalar k_rIJ = (r - rMin_)/dr_;

        label k = label(k_rIJ);

        scalar e = 0.0;

        if(k < forceLookup_.size()-1)
        {
            e = (k_rIJ - k)*energyLookup_[k+1]
                + (k + 1 - k_rIJ)*energyLookup_[k];
        }

        return e;
    }
}

List< Pair< scalar > >
    pairPotentialModel::energyTable() const
{
    List<Pair<scalar> > energyTab(energyLookup_.size());

    forAll(energyLookup_,k)
    {
        energyTab[k].first() = rMin_ + k*dr_;

        energyTab[k].second() = energyLookup_[k];
    }

    return energyTab;
}


scalar pairPotentialModel::scaledEnergy
(
    const scalar r
) const
{
    scalar e = unscaledEnergy(r);

    scaleEnergy(e, r);

    return e;
}


scalar pairPotentialModel::energyDerivative
(
    const scalar r,
    const bool scaledEnergyDerivative
) const
{
    // Local quadratic fit to energy: E = a0 + a1*r + a2*r^2
    // Differentiate to give f = -dE/dr = -a1 - 2*a2*r

    scalar ra = r - dr_;
    scalar rf = r;
    scalar rb = r + dr_;

    scalar Ea, Ef, Eb;

    if (scaledEnergyDerivative)
    {
        Ea = scaledEnergy(ra);
        Ef = scaledEnergy(rf);
        Eb = scaledEnergy(rb);
    }
    else
    {
        Ea = unscaledEnergy(ra);
        Ef = unscaledEnergy(rf);
        Eb = unscaledEnergy(rb);
    }

    scalar denominator = (ra - rf)*(ra - rb)*(rf - rb);

    scalar a1 =
    (
        rb*rb*(Ea - Ef) + ra*ra*(Ef - Eb) + rf*rf*(Eb - Ea)
    ) / denominator;

    scalar a2 =
    (
        rb*(Ef - Ea) + rf*(Ea - Eb) + ra*(Eb - Ef)
    ) / denominator;

    return a1 + 2.0*a2*r;
}


// bool pairPotentialModel::read
// (
//     const dictionary& pairPotentialProperties
// )
// {
//     pairPotentialProperties_ = pairPotentialProperties;
// 
//     return true;
// }

void pairPotentialModel::output(const fileName& pathName)
{
    writeEnergyAndForceTables(pathName);
    write(pathName);
}

void pairPotentialModel::writeEnergyAndForceTables(const fileName& pathName)
{
    Info<< "Writing energy and force to file for potential "
            << name_ << endl;
            
    label nBins = label((rCut_ - rMin_)/dr_) + 1;            
    
    scalarField U(nBins, 0.0);
    scalarField f(nBins, 0.0);
    
    for (label i=0; i<nBins; ++i)
    {
        scalar r = rMin_+dr_*i;
        
        U[i] = energy(r);
        f[i] = force(r);
    }
    {
        OFstream file(pathName/name_+"-RU.xy");

        if(file.good())
        {
            forAll(U, i)
            {
                file 
                    << dr_*i << "\t"
                    << U[i] << "\t"
                    << f[i]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void pairPotentialModel::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    
    {
        OFstream file(pathName/name_+"-SI.xy");

        if(file.good())
        {
            forAll(U, i)
            {
                file 
                    << dr_*i*rU_.refLength() << "\t"
                    << U[i]*rU_.refEnergy() << "\t"
                    << f[i]*rU_.refForce()
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void pairPotentialModel::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }  
    } 
}

const dictionary& pairPotentialModel::pairPotentialProperties() const
{
    return pairPotentialProperties_;
}



} // End namespace Foam

// ************************************************************************* //
