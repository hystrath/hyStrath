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

#include "lennardJonesEqn.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(lennardJonesEqn, 0);
addToRunTimeSelectionTable(pairPotentialModel, lennardJonesEqn, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
lennardJonesEqn::lennardJonesEqn
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
    sigma_(readScalar(propsDict_.lookup("sigma"))),
    epsilon_(readScalar(propsDict_.lookup("epsilon")))    
{
    if(redUnits.runReducedUnits())
    {
        sigma_ /= redUnits.refLength();
        epsilon_ /= redUnits.refEnergy();
    }

    useTables_ = false;
    

    F_at_Rmin_ = rawForce(rMin_);
    E_at_Rmin_ = rawEnergy(rMin_);    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lennardJonesEqn::~lennardJonesEqn()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar lennardJonesEqn::unscaledEnergy(const scalar r) const
{
    return 0.0;
}

scalar lennardJonesEqn::rawEnergy(const scalar r) const
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;

    return 4.0*epsilon_*(ir6*(ir6 - 1.0));
}

scalar lennardJonesEqn::rawForce(const scalar r) const
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;    
        
    return 24.0*epsilon_*ir6*(2.0*ir6 - 1.0)*sqrt(ir2);
}

   
scalar lennardJonesEqn::energy(const scalar r) const
{
    scalar energy = E_at_Rmin_;
    
    if(r > rMin_)
    {
        energy = rawEnergy(r);
    }    
    
    return energy;
}

scalar lennardJonesEqn::force(const scalar r) const
{
    scalar force = F_at_Rmin_;
    
    if(r > rMin_)
    {
        force = rawForce(r);
    }    
    
    return force;
}



// bool lennardJonesEqn::read
// (
//     const dictionary& pairPotentialProperties,
//     const reducedUnits& rU
// )
// {
//     pairPotentialModel::read(pairPotentialProperties, rU);
// 
//     lennardJonesEqnCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");
// 
//     lennardJonesEqnCoeffs_.lookup("sigma") >> sigma_;
//     lennardJonesEqnCoeffs_.lookup("epsilon") >> epsilon_;
// 
//     if(rU.runReducedUnits())
//     {
//         sigma_ /= rU.refLength();
//         epsilon_ /= rU.refEnergy();
//     }
// 
//     return true;
// }

void lennardJonesEqn::write(const fileName& pathName)
{
    
}

const dictionary& lennardJonesEqn::dict() const
{
    return propsDict_;
}


} // End namespace Foam

// ************************************************************************* //
