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

#include "harmonicPotential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(harmonicPotential, 0);
addToRunTimeSelectionTable(pairPotentialModel, harmonicPotential, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
harmonicPotential::harmonicPotential
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud, 
    const reducedUnits& redUnits,
    const word& name, 
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict),
    propsDict_(dict.subDict(typeName + "Coeffs")),
    k_(readScalar(propsDict_.lookup("k"))),
    r0_(readScalar(propsDict_.lookup("r0")))
{
    if(redUnits.runReducedUnits())
    {
        k_ /= (redUnits.refMass()/ (redUnits.refTime()*redUnits.refTime()));
        r0_ /= redUnits.refLength();
       
        Info << "k = " << k_ << endl;
    }
    
//     useTables_ = false;
    setLookupTables();    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

harmonicPotential::~harmonicPotential()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar harmonicPotential::unscaledEnergy(const scalar r) const
{
    return 0.5*k_*(r-r0_)*(r-r0_);
}

scalar harmonicPotential::energy(const scalar r) const
{
    return energyLookUpFromTable(r);  
}

scalar harmonicPotential::force(const scalar r) const
{
    return forceLookUpFromTable(r);
}

const dictionary& harmonicPotential::dict() const
{
    return propsDict_;
}

void harmonicPotential::write(const fileName& pathName)
{
//     Info<< "Writing energy and force to file for potential "
//             << name_ << endl;
//             
//     label nBins = 10000;
//     scalar dr = r0_/nBins;
//     scalarField U(nBins, 0.0);
//     scalarField f(nBins, 0.0);
//     
//     for (label i=0; i<nBins; ++i)
//     {
//         scalar r = dr*i;
//         
//         U[i] = energy(r);
//         f[i] = force(r);
//     }
//     {
//         OFstream file(pathName/name_+"-harmonicPotential-RU.xy");
// 
//         if(file.good())
//         {
//             forAll(U, i)
//             {
//                 file 
//                     << dr*i << "\t"
//                     << U[i] << "\t"
//                     << f[i]
//                     << endl;
//             }
//         }
//         else
//         {
//             FatalErrorIn("void harmonicPotential::write()")
//                 << "Cannot open file " << file.name()
//                 << abort(FatalError);
//         }
//     }
//     
//     {
//         OFstream file(pathName/name_+"-harmonicPotential-SI.xy");
// 
//         if(file.good())
//         {
//             forAll(U, i)
//             {
//                 file 
//                     << dr*i*rU_.refLength() << "\t"
//                     << U[i]*rU_.refEnergy() << "\t"
//                     << f[i]*rU_.refForce()
//                     << endl;
//             }
//         }
//         else
//         {
//             FatalErrorIn("void harmonicPotential::write()")
//                 << "Cannot open file " << file.name()
//                 << abort(FatalError);
//         }  
//     }
}


} // End namespace Foam

// ************************************************************************* //
