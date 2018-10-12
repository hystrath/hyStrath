/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "constBackground.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(constBackground, 0);
    addToRunTimeSelectionTable(ElectronModel, constBackground, dictionary);

};

Foam::constBackground::constBackground
(
    const dictionary& dict,
    pdCloud& cloud
)
:
    ElectronModel(cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    rhoB_
    (
    dimensionedScalar
        (
            "rhoB",
            dimensionSet(0, -3, 0, 0, 0, 0, 0),
            readScalar(coeffDict_.lookup("density"))
        )
    ),
    TB_
    (
    dimensionedScalar
        (
            "TB",
            dimensionSet(0, 0, 0, 1, 0, 0, 0),
            readScalar(coeffDict_.lookup("temperature"))
        )
    ),
    phi0_(
    dimensionedScalar
        (
            "phi0",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            readScalar(coeffDict_.lookup("potential"))
        )
    ),
    Zb_(readScalar(coeffDict_.lookup("charge")))
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::constBackground::~constBackground()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::constBackground::active() const
{
    return false;
}

void Foam::constBackground::calculateRhoe()
{
    //- constants
    const scalar e      = electromagnetic::e.value();           //- elementary charge
    const scalar eps0   = electromagnetic::epsilon0.value();    //- permitivity of space
    const scalar kB     = physicoChemical::k.value();           //- boltzmann constant

    //- non-constant reference to field objects
    pdEmFields& emF         = cloud_.emFields();
    pdStandardFields& stdF  = cloud_.stdFields();

    //- non-constant reference fields for modifying
    scalarField& rhoEF     = emF.rhoEF_.primitiveFieldRef();
    //scalarField& lambdaD   = emF.lambdaD_.primitiveFieldRef();

    //- constant reference to standard fields *DO NOT MODIFY HERE*
    //const scalarField transT  = stdF.transT();
    const scalarField rhoN    = stdF.rhoN();

    //- debye constant
    scalar K = eps0*kB/(e*e);

    /** calculate debye length **/
    forAll(rhoN, cI)
    {
        //- set background fluid
        rhoEF[cI]  = rhoB_.value();

        /*if(rhoN[cI] > VSMALL)
        {
            lambdaD[cI]    = Foam::sqrt( K / ( (rhoEF[cI]/TB_.value()) + Z*rhoN[cI]/transT[cI] ) );
        }
        else
        {
            if(rhoEF[cI] > VSMALL)
            {
                lambdaD[cI]    = Foam::sqrt( K / ( (rhoEF[cI]/TB_.value()) ) );
            }
            else
            {
                lambdaD[cI]    = 0.0; //this may give div(0) errors if used
                //lambdaD[cI]    = GREAT; //this would not average out
            }
        }*/
    }
}


void Foam::constBackground::calculatePotential()
{
    //- constants
    const dimensionedScalar e = electromagnetic::e;
    const dimensionedScalar eps0 = electromagnetic::epsilon0; // permitivity of space

    //- non-constant reference to electromagnetic field object
    pdEmFields& emF = cloud_.emFields();

    //- calculate background density
    calculateRhoe();

    /** assemble linear system **/
    fvScalarMatrix phiEqn
    (
        fvm::laplacian(emF.phiE())
        + (emF.rhoQ_ + emF.rhoEF_*Zb_*e)/eps0
    );

    //- if there are no fixed value of gradient BCs, force a
    phiEqn.setReference(0, phi0_.value()); //- <-- should add in option to change referance potential

    /** solve linear system **/
    phiEqn.solve();
    emF.phiE_.correctBoundaryConditions();

}

const Foam::dictionary& Foam::constBackground::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
