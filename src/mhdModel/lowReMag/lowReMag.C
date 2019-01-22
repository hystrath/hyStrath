/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "lowReMag.H"
#include "physicoChemicalConstants.H"
#include "conductivityModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
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
  Info << "Lowremag constructor started" << endl;
  calculateSigma();
}


lowReMag::lowReMag
(
    const dictionary& dict,
    const rho2ReactionThermo& thermo
)
:
    mhdModel(thermo),
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
  calculateSigma();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lowReMag::~lowReMag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool lowReMag::read()
{
    return mhdModel::read();
}
//- Calculating conductivity using the chosen model
void lowReMag::calculateSigma()
{
    Info << "Getting conductivity values from the model" << endl;
    sigma_ = conductivity_->sigma();
    Info << "Conductivity values retrieved. Max(sigma) = " << max(sigma_) << endl;
}
//- Calculating current density
tmp<volVectorField> lowReMag::j(const volVectorField& U) const
{
    
    if(mesh_.time().outputTime())
    {
        conductivity_->sigma().write();
    }

       volScalarField nDe = thermo_.composition().nD("e-");
    Info << max(nDe) << endl;
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
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );
    forAll(betaHall, cellI)
    {
        if(nDe[cellI] == 0.0)
        {
            betaHall[cellI] = 0.0;
        }
        else
        {
            betaHall[cellI] = sigma_[cellI]*mag(B_[cellI])/(1.69022e-19*nDe[cellI]);
        }
    }
    Info << min(betaHall) <<endl;

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
    Info << min(D) << endl;

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
     dimensionedTensor("hallParameter", dimensionSet(0,0,0,0,0,0,0), Tensor<double>(0, 0, 0, 0, 0, 0, 0, 0, 0))
   );
   Info<< "Hall tensor allocated\n" << endl;
           hallParameter.replace(0, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(0)))); // = sigma/D*(magSqr(B) + sqr(betaHall));
           hallParameter.replace(1, betaHall*(betaHall*B_.component(0)*B_.component(1) - mag(B_)*B_.component(2)));
           hallParameter.replace(2, betaHall*(betaHall*B_.component(0)*B_.component(2) + mag(B_)*B_.component(1)));
           hallParameter.replace(3, betaHall*(betaHall*B_.component(1)*B_.component(0) + mag(B_)*B_.component(2)));
           hallParameter.replace(4, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(1))));
           hallParameter.replace(5, betaHall*(betaHall*B_.component(1)*B_.component(2) - mag(B_)*B_.component(0)));
           hallParameter.replace(6, betaHall*(betaHall*B_.component(2)*B_.component(0) - mag(B_)*B_.component(1)));
           hallParameter.replace(7, betaHall*(betaHall*B_.component(2)*B_.component(1) + mag(B_)*B_.component(0)));
           hallParameter.replace(8, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(2))));
   forAll(hallParameter, cellI)
   {
       if (D[cellI] != 0)
       {
           hallParameter[cellI] /= D[cellI];
       }
       else
       {
           hallParameter[cellI] = Tensor<double>(0, 0, 0, 0, 0, 0, 0, 0, 0);
       }
   }
   Info << max(hallParameter) << endl;

    if (hallEffect_)
    {
    Info<< "Returning current density with Hall correction\n" << endl;
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
                    IOobject::AUTO_WRITE
                ),
                conductivity_->sigma()*hallParameter&(U^B_)
                //sigma_*(U^B_)
            )
        );
    }
    else
    {
    Info<< "Returning current density without Hall correction\n" << endl;
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
                    IOobject::AUTO_WRITE
                ),
                conductivity_->sigma()*(U^B_)
            )
        );
    }
}

//- Calculate Hall correction
tmp<volTensorField> lowReMag::hallCorrection() const
{
    volScalarField nDe = thermo_.composition().nD("e-");
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
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );
    forAll(betaHall, cellI)
    {
        if(nDe[cellI] == 0.0)
        {
            betaHall[cellI] = 0.0;
        }
        else
        {
            betaHall[cellI] = sigma_[cellI]*mag(B_[cellI])/(1.69022e-19*nDe[cellI]);
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
     dimensionedTensor("hallParameter", dimensionSet(0,0,0,0,0,0,0), Tensor<double>(0, 0, 0, 0, 0, 0, 0, 0, 0))
   );
   Info<< "Hall tensor allocated\n" << endl;
    /*       hallParameter.replace(0, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(0)))); // = sigma/D*(magSqr(B) + sqr(betaHall));
           hallParameter.replace(1, betaHall*(betaHall*B_.component(0)*B_.component(1) - mag(B_)*B_.component(2)));
           hallParameter.replace(2, betaHall*(betaHall*B_.component(0)*B_.component(2) + mag(B_)*B_.component(1)));
           hallParameter.replace(3, betaHall*(betaHall*B_.component(1)*B_.component(0) + mag(B_)*B_.component(2)));
           hallParameter.replace(4, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(1))));
           hallParameter.replace(5, betaHall*(betaHall*B_.component(1)*B_.component(2) - mag(B_)*B_.component(0)));
           hallParameter.replace(6, betaHall*(betaHall*B_.component(2)*B_.component(0) - mag(B_)*B_.component(1)));
           hallParameter.replace(7, betaHall*(betaHall*B_.component(2)*B_.component(1) + mag(B_)*B_.component(0)));
           hallParameter.replace(8, (magSqr(B_)+sqr(betaHall)*sqr(B_.component(2))));*/
   forAll(hallParameter, cellI)
   {
       if (D[cellI] != 0)
       {
           hallParameter[cellI] /= D[cellI];
       }
       else
       {
           hallParameter[cellI] = Tensor<double>(0, 0, 0, 0, 0, 0, 0, 0, 0);
       }
   }
   Info<< "Returning Hall tensor\n" << endl;
   return hallParameter;
   
}

//- Calculate MHD term for momentum equation
tmp<volVectorField> lowReMag::F(const volVectorField& U) const
{

   Info<< "Calculating momentum equation MHD term\n" << endl;
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
                IOobject::AUTO_WRITE
            ),
            //sigma_*(U*magSqr(B_) + B_*(U&B_))
            j(U)^B_
        )
    );
}

//- Calculate MHD term for energy equation
tmp<volScalarField> lowReMag::Q(const volVectorField& U) const
{
   Info<< "Calculating energy equation MHD term\n" << endl;
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
                IOobject::AUTO_WRITE
            ),
            //sigma_*(U^B_)&(U^B_)
            j(U)&(U^B_)
        )
    );
}

//- Calculate Stuart number equation
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
                IOobject::AUTO_WRITE
            ),
            magSqr(B_)*sigma_/(mag(U)*thermo_.rho())
        )
    );
}

} // End namespace mhd
} // End namespace FOAM

// ************************************************************************* //
