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
    fenePotential

Description

\*----------------------------------------------------------------------------*/

#include "fenePotential.H"
#include "OFstream.H"
#include "IFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void fenePotential::setRadius()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//        radius_[i] = (0.5 + scalar(i)) * binWidth();
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
fenePotential::fenePotential
(
    const polyMesh& mesh,
    const reducedUnits& redUnits
)
:
    mesh_(mesh),
    rU_(redUnits)
{
}

void fenePotential::setPotential()
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

    const word potentialName = "fene";

    const dictionary& dict = potentialDict.subDict(potentialName);

//     dictionary brennerCoeffs = dict.subDict(potentialName + "Coeffs");


    rCut_ = readScalar(dict.lookup("rCut"));

    dr_ = readScalar(dict.lookup("dr"));
    
    r0_ = readScalar(dict.lookup("Ro"));

    r0in_ = readScalar(dict.lookup("Roin"));
    
    k_ = readScalar(dict.lookup("k"));    

    // modify to reduced units
    


    r0_ /= rU_.refLength();
    r0in_ /= rU_.refLength();

    dr_ /= rU_.refLength();
    rCut_ /= rU_.refLength();


    rCutSqr_ = rCut_*rCut_;
    
    setLookUpTables();
    
    outputProperties();
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fenePotential::~fenePotential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fenePotential::setLookUpTables()
{
    nBins_ = label(rCut_ / dr_)+1; 
    
    U_.setSize(nBins_, 0.0);
    f_.setSize(nBins_, 0.0);
    
    for (label i=0; i<nBins_; ++i)
    {
        scalar r = dr_*i;
        
        U_[i] = Ufene(r);
        f_[i] = fFene(r);
    }
    
}
    
void fenePotential::outputProperties()
{
    if(Pstream::master())
    {
        // output fij
        Info << "output properties of fene potential " << endl;

        OFstream fileFij("fene");

        if(fileFij.good())
        {
            forAll(U_, i)
            {
                fileFij 
                    << dr_*i << "\t"
                    << U_[i] << "\t"
                    << f_[i]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void fenePotential::fenePotential()")
                << "Cannot open file " << fileFij.name()
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void fenePotential::operator=(const fenePotential& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("fenePotential::operator=(const fenePotential&)")
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

// Ostream& operator<<(Ostream& os, const fenePotential& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const fenePotential&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
