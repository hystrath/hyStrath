/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "Boltz.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(Boltz, 0);
    addToRunTimeSelectionTable(ElectronModel, Boltz, dictionary);

};

Foam::Boltz::Boltz
(
    const dictionary& dict,
    pdCloud& cloud
)
:
    ElectronModel(cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    n0_
    (
    dimensionedScalar
        (
              "n0",
              dimensionSet(0, -3, 0, 0, 0, 0, 0),
              readScalar(coeffDict_.lookup("density"))
        )
    ),
    TeV_
    (
            "TeV",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            0.0
    ),
    Te_
    (
    dimensionedScalar
        (
            "Te",
            dimensionSet(0, 0, 0, 1, 0, 0, 0),
            readScalar(coeffDict_.lookup("temperature"))
        )
    ),
    phi0_
    (
    dimensionedScalar
        (
            "phi0",
            dimensionSet(1, 2, -3, 0, 0, -1, 0),
            readScalar(coeffDict_.lookup("potential"))
        )
    ),
    convTol_(readScalar(coeffDict_.lookup("tol"))),
    nCorr_(readScalar(coeffDict_.lookup("nCorr"))),
    phiN_
    (
        IOobject
        (
            "phiN",
            cloud_.mesh().time().timeName(),
            cloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud_.mesh(),
        dimensionedScalar("zero",dimensionSet(1, 2, -3, 0, 0, -1, 0),0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            cloud_.mesh().time().timeName(),
            cloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud_.mesh(),
        dimensionedScalar("zero",dimensionSet(0, -2, 0, 0, 0, 0, 0),0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    dphi_
    (
        IOobject
        (
            "dphi",
            cloud_.mesh().time().timeName(),
            cloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        cloud_.mesh(),
        dimensionedScalar("zero",dimensionSet(1, 2, -3, 0, 0, -1, 0),0.0),
        zeroGradientFvPatchField<scalar>::typeName
    )

{
    const dimensionedScalar kB = physicoChemical::k;          //- boltzmann constant
    const dimensionedScalar e     = electromagnetic::e;       //- elementary charge

    TeV_ = Te_ * kB/e;

    Info << "   Referance Temperature:  " << Te_.value() << "K (" << TeV_.value() << " eV)" << endl;
    Info << "             Potential:    " << phi0_.value() << "V" << endl;
    Info << "             Density:      " << n0_.value() << " m^-3" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::Boltz::~Boltz()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::Boltz::active() const
{
    return false;
}


void Foam::Boltz::calculateRhoe()
{
    //- constants
    const scalar e      = electromagnetic::e.value();           //- elementary charge
    const scalar eps0   = electromagnetic::epsilon0.value();    //- permitivity of space
    const scalar kB     = physicoChemical::k.value();           //- boltzmann constant

    //- non-constant reference to field objects
    pdEmFields& emF         = cloud_.emFields();
    //pdStandardFields& stdF  = cloud_.stdFields();

    //- non-constant reference to electromagnetic fields to modify
    scalarField& phiE       = emF.phiE_.primitiveFieldRef();
    scalarField& rhoEF      = emF.rhoEF_.primitiveFieldRef();
    //scalarField& lambdaD   = emF.lambdaD_.primitiveFieldRef();

    //- constant reference to standard fields
    //const scalarField transT  = stdF.transT();
    //const scalarField rhoN    = stdF.rhoN();

    //- debye constant
    //scalar K = eps0*kB/(e*e);

    forAll(rhoEF,cI)
    {
        /** calculate electron fluid density **/
        rhoEF[cI]   = n0_.value()*Foam::exp( e/(kB*Te_.value()) * (phiE[cI] - phi0_.value()));

        /** calculate debye length **/
        /*if(rhoN[cI] > VSMALL && transT[cI] > VSMALL)
        {
            lambdaD[cI]    = Foam::sqrt( K / ( (rhoEF[cI]/Te_.value()) + Z*rhoN[cI]/transT[cI] ) );
        }
        else // error catch for 0 temperature where there are 0 particles
        {
            lambdaD[cI]    = Foam::sqrt( K / ( (rhoEF[cI]/Te_.value()) ) );
        }*/
    }

    //- correct boundary conditions across processors
    emF.rhoEF_.correctBoundaryConditions();
    //emF.lambdaD_.correctBoundaryConditions();
}


void Foam::Boltz::calculatePotential()
{
    //- constants
    const dimensionedScalar e     = electromagnetic::e;             //- elementary charge
    const dimensionedScalar eps0  = electromagnetic::epsilon0;      //- permitivity of space
    const dimensionedScalar kB     = physicoChemical::k;            //- boltzmann constant

    //- non-constant refernce to electromagnetic fields
    pdEmFields& emF = cloud_.emFields();

    //- electron temperature (assumed constant)
    dimensionedScalar alpha = e/(kB*Te_);

    //- initialise computational variables
    label iCorr = 0;                    //- outer iteration counter
    Foam::solverPerformance solverPerf; //- reference to solver performance monitor
    scalar initialResidual = 1.0;       //- initial residual
    scalar delta = 1.0;

    /**Non-linear Poisson-Boltzmann Equations**/
    /* Newton-Raphson Method - page 172 Computer simulations using particles
    N(phi) = \nabla^2 \phi +  e/eps0 (n_i - n_0_e*exp(e/(kB*T) \phi)

    \frac{\partial N}{\partial \phi} = \nabla^2 - n_0_e*e^2/(eps0*kB*T)*exp(e/(kB*T) \phi)
                                     = \nabla^2 - 1/\lambda_D*exp(e/(kB*T) \phi)

    solve for phi^{n+1}

    \left( \frac{\partial N}{\partial \phi} \right) ^ n \phi^{n+1} = -N(\phi^n) + \left( \frac{\partial N}{\partial \phi} \right) ^ n \phi^{n}
    */

    Info   << "Solving for " << emF.phiE().name() << endl;

    //- set initial guess as previous field
    phiN_  = emF.phiE_;
    phiN_.correctBoundaryConditions();
    do
    {

        //- dN/dphi constant
        kappa_ = n0_*alpha*e/eps0 * Foam::exp(alpha * (phiN_ - phi0_) );
        kappa_.correctBoundaryConditions();

        /** assemble linear system **/
        fvScalarMatrix phiEqn
        (
            fvm::laplacian(emF.phiE_)
            - fvm::SuSp(kappa_,emF.phiE_)
            ==
            - emF.rhoQ_/eps0
            + e*n0_/eps0 * (1 - alpha*phiN_ )* Foam::exp(alpha * (phiN_ - phi0_ ))
        );

        //- if there are no fixed value of gradient BCs, force a
        phiEqn.setReference(0, 0); //- TODO: add in option to change reference potential

        /** solve system **/
        solverPerf = solve(phiEqn); //phiE^{n+1}
        emF.phiE_.correctBoundaryConditions();

        /** calculate error **/
        dphi_ = emF.phiE_ - phiN_;
        dphi_.correctBoundaryConditions();

        /** save phi^N **/
        phiN_ = emF.phiE_;

        //delta = max(mag(dphi_));
        delta = gMax(mag(dphi_.primitiveField()));
        /**----------------------------**/

        // initial outer residual
        if (iCorr == 0)
        {
            //initialResidual = delta.value();
            initialResidual = delta;
        }
    }
    while //- iterate until converged
    (
         delta > convTol_
            &&
        ++iCorr < nCorr_
    );

    calculateRhoe();

    // output final data
    //TODO: Add verbosity option
    Info   << ", Initial residual = "  << initialResidual
           << ", Final residual = "    << delta
           << ", No outer iterations " << iCorr << endl;
}

const Foam::dictionary& Foam::Boltz::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
