/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    brennerPotential

Description

\*----------------------------------------------------------------------------*/

#include "brennerPotential.H"
#include "OFstream.H"
#include "IFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void brennerPotential::setRadius()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//        radius_[i] = (0.5 + scalar(i)) * binWidth();
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
brennerPotential::brennerPotential
(
    const polyMesh& mesh,
    const reducedUnits& redUnits
)
:
    mesh_(mesh),
    rU_(redUnits)
{
    setPotential();
}

void brennerPotential::setPotential()
{
    IOdictionary potentialDict
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word potentialName = "brenner";

    const dictionary& dict = potentialDict.subDict(potentialName);

//     dictionary brennerCoeffs = dict.subDict(potentialName + "Coeffs");


    sIJ_ = readScalar(dict.lookup("sIJ"));
    cIJ_ = readScalar(dict.lookup("cIJ"));
    betaIJ_ = readScalar(dict.lookup("betaIJ"));
    rIJ0_ = readScalar(dict.lookup("rIJ0"));
    rIJ1_ = readScalar(dict.lookup("rIJ1"));
    rIJ2_ = readScalar(dict.lookup("rIJ2"));
    deltaI_ = readScalar(dict.lookup("deltaI"));
    alphaIJK_ = readScalar(dict.lookup("alphaIJK"));
    a0_ = readScalar(dict.lookup("a0"));
    c0_ = readScalar(dict.lookup("c0"));
    d0_ = readScalar(dict.lookup("d0"));

    dr_= readScalar(dict.lookup("dr"));
    rMin_= readScalar(dict.lookup("rMin"));
    rCut_= readScalar(dict.lookup("rCut"));


    // modify to reduced units
    

    cIJ_ /= rU_.refEnergy();
    betaIJ_ *= rU_.refLength();
    rIJ0_ /= rU_.refLength();
    rIJ1_ /= rU_.refLength();
    rIJ2_ /= rU_.refLength();

    dr_ /= rU_.refLength();
    rMin_ /= rU_.refLength();
    rCut_ /= rU_.refLength();


    constantRepulsive_ = -sqrt(2.0*sIJ_)*betaIJ_;
    constantAttractive_ = -sqrt(2.0/sIJ_)*betaIJ_;
    constantBrenner_ = cIJ_/(sIJ_-1.0);

    outputProperties();

    Info << nl << "brenner properties in reduced units: " << endl;

    Info<< "rIJ0 = " << rIJ0_ 
        << nl << "rIJ1 = " << rIJ1_
        << nl << "rIJ2 = " << rIJ2_
        << nl << "betaIJ = " << betaIJ_
        << nl << "cIJ = " << cIJ_
        << nl << "constantBrenner c/(s-1) = " << constantBrenner_
        << endl;

}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

brennerPotential::~brennerPotential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void brennerPotential::outputProperties()
{
    // output fij
    label nBins = label((rCut_-rMin_)/dr_);

    scalarField rField(nBins, 0.0);
    scalarField fIJField(nBins, 0.0);
    
    for (label i=0; i<nBins; ++i)
    {
        rField[i] = rMin_+(i*dr_);
        fIJField[i]=fIJ(rField[i]);
    }

    OFstream fileFij("fij");

    if(fileFij.good())
    {
        forAll(rField, n)
        {
            fileFij 
                << rField[n] << "\t" 
                << fIJField[n]
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void brennerPotential::brennerPotential()")
            << "Cannot open file " << fileFij.name()
            << abort(FatalError);
    }


    // output repulsive part of potential

    scalarField URepField(nBins, 0.0);
    scalarField fRepField(nBins, 0.0);

    for (label i=0; i<nBins; ++i)
    {
        const scalar& r = rField[i];
        URepField[i]=Urep(r);
        fRepField[i]=fRep(r);
    }

    OFstream fileUrep("Urepulsive");

    if(fileUrep.good())
    {
        forAll(rField, n)
        {
            fileUrep 
                << rField[n] << "\t" 
                << URepField[n] << "\t" 
                << fRepField[n]
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void brennerPotential::brennerPotential()")
            << "Cannot open file " << fileUrep.name()
            << abort(FatalError);
    }


    // output attractive part of potential 

    scalarField UAttField(nBins, 0.0);
    scalarField fAttField(nBins, 0.0);

    for (label i=0; i<nBins; ++i)
    {
        const scalar& r = rField[i];
        UAttField[i]=-Uatt(r);
        fAttField[i]=-fAtt(r);
    }

    OFstream fileUatt("Uattractive");

    if(fileUatt.good())
    {
        forAll(rField, n)
        {
            fileUatt 
                << rField[n] << "\t" 
                << UAttField[n] << "\t" 
                << fAttField[n]
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void brennerPotential::brennerPotential()")
            << "Cannot open file " << fileUatt.name()
            << abort(FatalError);
    }

    // estimate of overall potential

    scalarField UField(nBins, 0.0);
    scalarField fField(nBins, 0.0);

    // arbitrarily choosing bond order term to be equal to 1.0

    for (label i=0; i<nBins; ++i)
    {
        const scalar& r = rField[i];
        UField[i]=Urep(r)-Uatt(r);
        fField[i]=fRep(r)-fAtt(r);
    }

    OFstream fileU("Ubrenner");

    if(fileU.good())
    {
        forAll(rField, n)
        {
            fileU 
                << rField[n] << "\t" 
                << UField[n] << "\t" 
                << fField[n]
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void brennerPotential::brennerPotential()")
            << "Cannot open file " << fileU.name()
            << abort(FatalError);
    }

}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void brennerPotential::operator=(const brennerPotential& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("brennerPotential::operator=(const brennerPotential&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// 
//     Map<label>::operator=(rhs);
// 
//     binWidth_ = rhs.binWidth();
// }


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const brennerPotential& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const brennerPotential&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
