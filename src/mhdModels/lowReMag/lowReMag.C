/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigma", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0)
    )
{
    Info << "Reading electric conductivity field" << endl;
    
    electricalConductivity_->update();
    Info<< "Max(sigma) = "
        << gMax(electricalConductivity_->sigma()) << endl;
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigma", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0)
    )
{
    electricalConductivity_->update();
    Info<< "Max(sigma) = "
        << gMax(electricalConductivity_->sigma()) << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lowReMag::~lowReMag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool lowReMag::read()
{
    return mhdModel::read();
}


void lowReMag::update()
{
    //- Update electrical conductivity
    electricalConductivity_->update();
}


tmp<volVectorField> lowReMag::j(const volVectorField& U) const
{
    volScalarField betaHall
    (
        IOobject
        (
            "betaHall",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "betaHall",
            dimless,
            0.0
        )
    );
    
    volTensorField hallParameter
    (
        IOobject
        (
            "hallParameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "hallParameter",
            dimless,
            Tensor<scalar>(0, 0, 0, 0, 0, 0, 0, 0, 0)
        )
    );
    
    forAll(betaHall, cellI)
    {
        const scalar nDe = pe_[cellI]/(localkB_*T_[cellI]);
        
        if (nDe > SMALL)
        {
            betaHall[cellI] = sigma_[cellI]*mag(B_[cellI])
              / (localElecCharge_*nDe);
        }
    }
    Info << gMin(betaHall) <<endl;

    Info<< "Calculating D\n" << endl;
    volScalarField D
    (
        IOobject
        (
          "D",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        magSqr(B_)*(1.0+sqr(betaHall))
    );
    Info << gMin(D) << endl;

    
    hallParameter.replace(0, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(0)))); 
    
    hallParameter.replace
    (
        1,
        betaHall
      * (
            betaHall*B_.component(0)*B_.component(1) - mag(B_)*B_.component(2)
        )
    );
    hallParameter.replace
    (
        2,
        betaHall
      * (
            betaHall*B_.component(0)*B_.component(2) + mag(B_)*B_.component(1)
        )
    );
    hallParameter.replace
    (
        3,
        betaHall
      * (
            betaHall*B_.component(1)*B_.component(0) + mag(B_)*B_.component(2)
        )
    );
    hallParameter.replace
    (
        4,
        magSqr(B_)+sqr(betaHall)*sqr(B_.component(1))
    );
    hallParameter.replace
    (
        5,
        betaHall
      * (
            betaHall*B_.component(1)*B_.component(2) - mag(B_)*B_.component(0)
        )
    );
    hallParameter.replace
    (
        6,
        betaHall
      * (
            betaHall*B_.component(2)*B_.component(0) - mag(B_)*B_.component(1)
        )
    );
    hallParameter.replace
    (
        7,
        betaHall
      * (
            betaHall*B_.component(2)*B_.component(1) + mag(B_)*B_.component(0)
        )
    );
    hallParameter.replace(8, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(2))));
   
    forAll(hallParameter, cellI)
    {
        if (D[cellI] != 0.0)
        {
            hallParameter[cellI] /= D[cellI];
        }
        else
        {
            hallParameter[cellI] = Tensor<scalar>(0, 0, 0, 0, 0, 0, 0, 0, 0);
        }
    }
    Info << gMax(hallParameter) << endl;

    if (hallEffect_)
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "j",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                electricalConductivity_->sigma()*hallParameter&(U^B_)
            )
        );
    }
    else
    {
        return tmp<volVectorField>
        (
            new volVectorField
            (
                IOobject
                (
                    "j",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                electricalConductivity_->sigma()*(U^B_)
            )
        );
    }
}


tmp<volTensorField> lowReMag::hallCorrection() const
{
    volScalarField betaHall
    (
        IOobject
        (
            "betaHall",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "betaHall",
            dimless,
            0.0
        )
    );
    
    forAll(betaHall, cellI)
    {
        const scalar nDe = pe_[cellI]/(localkB_*T_[cellI]);
        
        if (nDe > SMALL)
        {
            betaHall[cellI] = sigma_[cellI]*mag(B_[cellI])
                /(localElecCharge_*nDe);
        }
    }

    Info<< "Calculating D\n" << endl;
    volScalarField D
    (
        IOobject
        (
            "D",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        magSqr(B_)*(1.0+sqr(betaHall))
    );

    Info<< "Calculating Hall tensor\n" << endl;
    volTensorField hallParameter
    (
        IOobject
        (
            "hallParameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "hallParameter",
            dimless,
            Tensor<scalar>(0, 0, 0, 0, 0, 0, 0, 0, 0)
        )
    );
   
    forAll(hallParameter, cellI)
    {
        if (D[cellI] != 0.0)
        {
            hallParameter[cellI] /= D[cellI];
        }
        else
        {
            hallParameter[cellI] = Tensor<scalar>(0, 0, 0, 0, 0, 0, 0, 0, 0);
        }
    }

    return hallParameter;
}


tmp<volVectorField> lowReMag::F(const volVectorField& U) const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            j(U)^B_
        )
    );
}


tmp<volScalarField> lowReMag::Q(const volVectorField& U) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Q",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            j(U)&(U^B_)
        )
    );
}


tmp<volScalarField> lowReMag::Stuart(const volVectorField& U) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Stuart",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            magSqr(B_)*sigma_/(mag(U)*thermo_.rho())
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// ************************************************************************* //
