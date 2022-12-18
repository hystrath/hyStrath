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

#include "lowReMag.H"
#include "electricalConductivityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(lowReMag, 0);
        addToMhdRunTimeSelectionTables(lowReMag);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void lowReMag::updateElectricCurrentDensity(const volVectorField& U)
{
    //- Ohm's law
    j_ = sigma()&(E_ + (U^B_));
}


tmp<volTensorField> lowReMag::hallTensor() const
{
    tmp<volTensorField> thallTensor
    (
        new volTensorField
        (
            IOobject
            (
                "hall",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor
            (
                "nullTensor",
                dimless,
                tensor::zero
            )
        )
    );
    
    volTensorField& hall = thallTensor.ref();
    
    volScalarField sigma = electricalConductivity_->sigma();
    
    forAll(hall, cellI)
    {
        scalar beta = 0.0;
        
        const scalar magB = mag(B_[cellI]);
        const scalar magSqrB = sqr(magB);
        
        if (constBeta_ >= 0.0)
        {
            beta = constBeta_;
        }
        else
        {
            const scalar nDe = pe_[cellI]/(localkB_*T_[cellI]);
            if (nDe > SMALL)
            {
                beta = sigma[cellI]*magB/(localElecCharge_*nDe);
            }
        }
        
        const scalar sqrBeta = sqr(beta);
        const scalar Bx = B_[cellI].component(0);
        const scalar By = B_[cellI].component(1);
        const scalar Bz = B_[cellI].component(2);
        
        const scalar D = magSqrB*(1.0 + sqrBeta);
        
        if (D > VSMALL)
        {
            tensor& hParam = hall[cellI];
            
            hParam.replace(0, magSqrB + sqrBeta*sqr(Bx)); 
            hParam.replace(1, beta*(beta*By*Bx - magB*Bz));
            hParam.replace(2, beta*(beta*Bz*Bx + magB*By));

            hParam.replace(3, beta*(beta*By*Bx + magB*Bz));
            hParam.replace(4, magSqrB + sqrBeta*sqr(By));
            hParam.replace(5, beta*(beta*Bz*By + magB*Bx));
            
            hParam.replace(6, beta*(beta*Bz*Bx - magB*By));
            hParam.replace(7, beta*(beta*Bz*By + magB*Bx));
            hParam.replace(8, magSqrB + sqrBeta*sqr(Bz));
            
            hParam /= D;
        }
    }

    return thallTensor;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lowReMag::lowReMag(const rho2ReactionThermo& thermo)
:
    mhdModel(thermo),
    T_(thermo.T()),
    pe_(thermo.pe()),
    localkB_(constant::physicoChemical::k.value()),
    localElecCharge_(constant::electromagnetic::e.value()),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(1, 0, -2, 0, 0, -1, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "E",
            dimensionSet(1, 1, -3, 0, 0, -1, 0),
            vector::zero
        )
    ),
    elecPot_
    (
        IOobject
        (
            "elecPot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    j_
    (
        IOobject
        (
            "j",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(0, -2, 0, 0, 0, 1, 0)
    )
{
    initialise();
    
    if (active_)
    {
        B_.writeOpt() = IOobject::AUTO_WRITE;
        E_.writeOpt() = IOobject::AUTO_WRITE;
        elecPot_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


lowReMag::lowReMag
(
    const dictionary& dict,
    const rho2ReactionThermo& thermo
)
:
    mhdModel(thermo),
    T_(thermo.T()),
    pe_(thermo.pe()),
    localkB_(constant::physicoChemical::k.value()),
    localElecCharge_(constant::electromagnetic::e.value()),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(1, 0, -2, 0, 0, -1, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "E",
            dimensionSet(1, 1, -3, 0, 0, -1, 0),
            vector::zero
        )
    ),
    elecPot_
    (
        IOobject
        (
            "elecPot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    j_
    (
        IOobject
        (
            "j",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionSet(0, -2, 0, 0, 0, 1, 0)
    )
{
    initialise();
    
    if (active_)
    {
        B_.writeOpt() = IOobject::AUTO_WRITE;
        E_.writeOpt() = IOobject::AUTO_WRITE;
        elecPot_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lowReMag::~lowReMag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool lowReMag::read()
{
    return mhdModel::read();
}


void lowReMag::update(const volVectorField& U)
{
    //- Update electrical conductivity
    electricalConductivity_->update();
    
    //- Update electric current density
    updateElectricCurrentDensity(U);
}


tmp<volScalarField> lowReMag::jouleHeating(const volVectorField& U) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "jouleHeating",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            j_&E_
        )
    );
}


tmp<volVectorField> lowReMag::lorentzForce() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "lorentzForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            j_^B_
        )
    );
}


volTensorField lowReMag::sigma() const
{
    if (!electricalConductivity_.valid())
    {
        FatalErrorIn
        (
            "const Foam::mhd::electricalConductivityModel&"
            "Foam::mhd::mhdModel::electricalConductivity() const"
        )   << "Requested electrical conductivity model, but model is "
            << "not activated" << abort(FatalError);
    }
    
    if (hallEffect_)
    {
        //- tensor electric conductivity
        return electricalConductivity_->sigma()*hallTensor();
    }
    else
    {
        //- scalar electric conductivity
        return electricalConductivity_->sigma()*tensor::I;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// ************************************************************************* //
