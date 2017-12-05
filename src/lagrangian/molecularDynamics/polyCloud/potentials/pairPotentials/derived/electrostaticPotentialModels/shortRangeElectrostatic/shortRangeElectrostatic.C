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

#include "shortRangeElectrostatic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(shortRangeElectrostatic, 0);
addToRunTimeSelectionTable(pairPotentialModel, shortRangeElectrostatic, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
shortRangeElectrostatic::shortRangeElectrostatic
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
    constant_(1.0/(4.0 * constant::mathematical::pi * 8.854187817e-12)),
    G_(readScalar(propsDict_.lookup("G")))    
{
 
    if(redUnits.runReducedUnits())
    {
        constant_ = (1.0/(4.0 * constant::mathematical::pi * redUnits.epsilonPermittivity()));
    }
    else
    {
        constant_ = 1.0/(4.0*constant::mathematical::pi*8.854187817e-12);
    }
    
    // we can override this in the future 
    
    useTables_ = false;

//     setLookupTables();   
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

shortRangeElectrostatic::~shortRangeElectrostatic()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar shortRangeElectrostatic::unscaledEnergy(const scalar r) const
{
    return constant_/r;
}



scalar shortRangeElectrostatic::force(const scalar r) const
{
    scalar force = 0.0;
    
    scalar term1 = erfc(G_*r/sqrt(2.0));
    scalar term2 = 2.0*G_*r*exp(-r*r)/sqrt(2* constant::mathematical::pi);
    force = constant_*(term1 + term2)/(r*r);
    
    return force;
}
    
scalar shortRangeElectrostatic::energy(const scalar r) const
{
    scalar energy = constant_*erfc(G_*r/sqrt(2.0))/r;
    
    return energy;
}

const dictionary& shortRangeElectrostatic::dict() const
{
    return pairPotentialProperties_;
}

void  shortRangeElectrostatic::write(const fileName& pathName)
{
    Info<< "Writing energy and force to file for potential "
            << name_ << endl;
            
    label nBins = label((rCut_ - rMin_)/dr_) + 1;   
    
//     scalar dr = (rCut_-rMin_)/nBins;
    scalarField U(nBins, 0.0);
    scalarField f(nBins, 0.0);
    
    for (label i=0; i<nBins; ++i)
    {
        scalar r = rMin_+dr_*i;
        
        U[i] = energy(r);
        f[i] = force(r);
    }
    {
        OFstream file(pathName/name_+"-electrostatics-RU.xy");

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
            FatalErrorIn("void shortRangeElectrostatic::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    
    {
        OFstream file(pathName/name_+"-electrostatics-SI.xy");

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
            FatalErrorIn("void shortRangeElectrostatic::write()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }  
    }    
}

} // End namespace Foam

// ************************************************************************* //
