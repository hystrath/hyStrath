/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/
// #include "mixedFixedValueSlipFvPatchFields.H"
#include "viscosityModel.H"

#include "nonLinMaxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinMaxwellSlipUFvPatchVectorField::nonLinMaxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    Uwall_(p.size(), vector(0.0, 0.0, 0.0)),
//    refSlipLen_(1.0),
    rho2_(1.0),
    rho1_(1.0),
    rho0_(1.0),
//    gammaC5_(1.0),
//    gammaC4_(1.0), 
//    gammaC3_(1.0),
    gammaC2_(1.0),
    gammaC1_(1.0) 
//    gammaC0_(1.0) 
//    p00_(1.0),
//    p10_(1.0),
//    p01_(1.0),
//    p20_(1.0),
//    p11_(1.0), 
//    p30_(1.0),
//    p21_(1.0)
    // frictionCoeff_(1.0)    
{}


nonLinMaxwellSlipUFvPatchVectorField::nonLinMaxwellSlipUFvPatchVectorField
(
    const nonLinMaxwellSlipUFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFixedValueSlipFvPatchVectorField(tdpvf, p, iF, mapper),
    Uwall_(tdpvf.Uwall_),
//    refSlipLen_(tdpvf.refSlipLen_),
    rho2_(tdpvf.rho2_),
    rho1_(tdpvf.rho1_),
    rho0_(tdpvf.rho0_),
//    gammaC5_(tdpvf.gammaC5_),
//    gammaC4_(tdpvf.gammaC4_), 
//    gammaC3_(tdpvf.gammaC3_),
    gammaC2_(tdpvf.gammaC2_), 
    gammaC1_(tdpvf.gammaC1_)
//    gammaC0_(tdpvf.gammaC0_) 
//    p00_(tdpvf.p00_),
//    p10_(tdpvf.p10_),
//    p01_(tdpvf.p01_), 
//    p20_(tdpvf.p20_),
//    p11_(tdpvf.p11_), 
//    p30_(tdpvf.p30_),
//    p21_(tdpvf.p21_) 
    // frictionCoeff_(tdpvf.frictionCoeff_)
{}


nonLinMaxwellSlipUFvPatchVectorField::nonLinMaxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    // accommodationCoeff_(readScalar(dict.lookup("accommodationCoeff"))),
    Uwall_("Uwall", dict, p.size()),
//    refSlipLen_(readScalar(dict.lookup("refSlipLen"))),
    rho2_(readScalar(dict.lookup("rho2"))),
    rho1_(readScalar(dict.lookup("rho1"))),
    rho0_(readScalar(dict.lookup("rho0"))),
//    gammaC5_(readScalar(dict.lookup("gammaC5"))),
//    gammaC4_(readScalar(dict.lookup("gammaC4"))),
//    gammaC3_(readScalar(dict.lookup("gammaC3"))),
    gammaC2_(readScalar(dict.lookup("gammaC2"))),
    gammaC1_(readScalar(dict.lookup("gammaC1")))
//    gammaC0_(readScalar(dict.lookup("gammaC0")))
//    p00_("p00", dict, p.size()),
//    p10_("p10", dict, p.size()),
//    p01_("p01", dict, p.size()),
//    p20_("p20", dict, p.size()),
//    p11_("p11", dict, p.size()),
//    p30_("p30", dict, p.size()),
//    p21_("p21", dict, p.size())
//    p00_(readScalar(dict.lookup("p00"))),
//    p10_(readScalar(dict.lookup("p10"))),
//    p01_(readScalar(dict.lookup("p01"))),
//    p20_(readScalar(dict.lookup("p20"))),
//    p11_(readScalar(dict.lookup("p11"))),
//    p30_(readScalar(dict.lookup("p30"))),
//    p21_(readScalar(dict.lookup("p21")))

    // frictionCoeff_(readScalar(dict.lookup("frictionCoeff")))
{

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
        refValue() = vectorField("refValue", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        mixedFixedValueSlipFvPatchVectorField::evaluate();
    }
}


nonLinMaxwellSlipUFvPatchVectorField::nonLinMaxwellSlipUFvPatchVectorField
(
    const nonLinMaxwellSlipUFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(tdpvf, iF),
    Uwall_(tdpvf.Uwall_),
//    refSlipLen_(tdpvf.refSlipLen_),
    rho2_(tdpvf.rho2_),
    rho1_(tdpvf.rho1_),
    rho0_(tdpvf.rho0_),
//    gammaC5_(tdpvf.gammaC5_),
//    gammaC4_(tdpvf.gammaC4_),
//    gammaC3_(tdpvf.gammaC3_),
    gammaC2_(tdpvf.gammaC2_),
    gammaC1_(tdpvf.gammaC1_)
//    gammaC0_(tdpvf.gammaC0_)
//    p00_(tdpvf.p00_),
//    p10_(tdpvf.p10_),
//    p01_(tdpvf.p01_),
//    p20_(tdpvf.p20_),
//    p11_(tdpvf.p11_),
//    p30_(tdpvf.p30_),
//    p21_(tdpvf.p21_)
    // frictionCoeff_(tdpvf.frictionCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void nonLinMaxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volVectorField& U =
        this->db().objectRegistry::lookupObject<volVectorField>("U");
    const volScalarField& rho =
        this->db().objectRegistry::lookupObject<volScalarField>("rho");
    label patchi = this->patch().index();
    vectorField nHat = this->patch().nf();
    vectorField shear = (nHat & fvc::grad(U)().boundaryField()[patchi]) & (I - sqr(nHat));

//    scalarField shear = strainRate().boundaryField[patchi];
//    volScalarField shear2 = sqrt(2.0)*mag(symm(fvc::grad(U)));
//    scalarField shear = shear2.boundaryField[patchi];

//Info << mag(shear) << endl;

    scalar maxFaceI = fvc::grad(U)().boundaryField()[patchi].size();



//     LJ Cases

// linear density dependance

    for(scalar faceI = 0; faceI < maxFaceI ; faceI++)
    {
        if (mag(shear[faceI]) < (gammaC2_*rho[faceI]+gammaC1_))
        {     
            valueFraction()[faceI] = 
		(1.0/(1.0 + patch().deltaCoeffs()[faceI]
		*(rho2_* pow(rho[faceI],2) + rho1_*rho[faceI]+rho0_)*
	(1.0/sqrt(1.0 - (mag(shear[faceI])/gammaC1_)))));
        }
        else
        {
          valueFraction()[faceI] = 0.0; 
        }

    }

////	Water Case

// no strain dependance

//    for(scalar faceI = 0; faceI < maxFaceI ; faceI++)
//    {
//        valueFraction()[faceI] = (1.0/(1.0 + patch().deltaCoeffs() / (rho1_*rho[faceI] + rho0_)));
//    }	




    refValue() = Uwall_;

    mixedFixedValueSlipFvPatchVectorField::updateCoeffs();
}


// Write
void nonLinMaxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    Uwall_.writeEntry("Uwall", os);
     os.writeKeyword("rho2")
        << rho2_ << token::END_STATEMENT << nl;
     os.writeKeyword("rho1")
        << rho1_ << token::END_STATEMENT << nl;
     os.writeKeyword("rho0")
        << rho0_ << token::END_STATEMENT << nl;
     os.writeKeyword("gammaC2")
        << gammaC2_ << token::END_STATEMENT << nl;
   os.writeKeyword("gammaC1")
        << gammaC1_ << token::END_STATEMENT << nl;

    refValue().writeEntry("refValue", os);
    valueFraction().writeEntry("valueFraction", os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, nonLinMaxwellSlipUFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
