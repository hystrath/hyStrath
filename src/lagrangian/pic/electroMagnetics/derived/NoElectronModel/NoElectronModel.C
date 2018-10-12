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

#include "NoElectronModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(NoElectronModel, 0);
    addToRunTimeSelectionTable(ElectronModel, NoElectronModel, dictionary);

};



Foam::NoElectronModel::NoElectronModel
(
    const dictionary& dict,
    pdCloud& cloud
)
:
    ElectronModel(cloud)
{
    Info    << "No electron model selected. Running Full-PD method" << nl
            << "Remember to intialise electrons correctly as particles." << nl << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::NoElectronModel::~NoElectronModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::NoElectronModel::active() const
{
    return false;
}

void Foam::NoElectronModel::calculateRhoe()
{

    //- constants
    //const scalar e      = electromagnetic::e.value();           //- elementary charge
    //const scalar eps0   = electromagnetic::epsilon0.value();    //- permitivity of space
    //const scalar kB     = physicoChemical::k.value();           //- boltzmann constant

    //- non-constant referance to field objects
    //pdEmFields& emF         = cloud_.emFields();
    //pdStandardFields& stdF  = cloud_.stdFields();

    //- non-constant referance fields for modifying
    //scalarField& lambdaD   = emF.lambdaD_.internalField();

    //- constant referance to standard fields *DO NOT MODIFY HERE*
    //const scalarField transT  = stdF.transT();
    //const scalarField rhoN    = stdF.rhoN();

    //- debye constant
    //scalar K = eps0*kB/(e*e);

    //- calculate debye length
    /** calculate debye length **/
    /*forAll(rhoN, cI)
    {
        if(rhoN[cI] > VSMALL)
        {
            lambdaD[cI]    = Foam::sqrt( K / ( Z*rhoN[cI]/transT[cI] ) );
        }
        else
        {
            lambdaD[cI]    = 0.0; //this may give div(0) errors if used
            //lambdaD[cI]    = GREAT; //this would not average out
        }
    }*/
}


void Foam::NoElectronModel::calculatePotential()
{
    //- constants
    const dimensionedScalar eps0 = electromagnetic::epsilon0; // permitivity of space

    //- non-constant referance to electromagnetic field object
    pdEmFields& emF = cloud_.emFields();

    //- calculate debye length
    //calculateRhoe();

    /** assemble linear equation system **/
    fvScalarMatrix phiEqn
    (
        fvm::laplacian(emF.phiE()) + emF.rhoQ_/eps0
    );

    //- if there are no fixed value of gradient BCs, force a
    phiEqn.setReference(0, 0); //- <-- should add in option to change referance potential

    /** solve linear system **/
    phiEqn.solve();
    emF.phiE_.correctBoundaryConditions();

}

const Foam::dictionary& Foam::NoElectronModel::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
