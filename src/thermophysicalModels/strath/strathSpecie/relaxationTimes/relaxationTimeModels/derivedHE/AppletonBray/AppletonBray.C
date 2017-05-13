/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "AppletonBray.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::AppletonBray<ThermoType>::AppletonBray
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    relaxationTimeModelHE(thermo, turbulence),
    
    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),
    
    electronListPosition_(species()["e-"]),
    
    RR(constant::physicoChemical::R.value()),
    NA(constant::physicoChemical::NA.value()),
    kB(RR/NA),
    ec(Foam::constant::electromagnetic::e.value()),
    pi(constant::mathematical::pi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::AppletonBray<ThermoType>::correct()
{
    // Scalabrin PhD thesis, Eq.(2.65)
    const volScalarField& Tt = thermo_.Tt();
    const volScalarField& pDe = thermo_.composition().pD(electronListPosition_);
    const volScalarField& nDe = thermo_.composition().nD(electronListPosition_);
    const volScalarField& Te = thermo_.composition().Tv(electronListPosition_);
    
    const scalarField& TtCells = Tt.internalField();
    const scalarField& pDeCells = pDe.internalField();
    const scalarField& nDeCells = nDe.internalField();
    const scalarField& TeCells = Te.internalField();
    scalarField& QHECells = this->QHE_.internalField();
    
    QHECells = 0.0;
    
    forAll(Tt.boundaryField(), patchi)
    {
        fvPatchScalarField& pQHE = this->QHE_.boundaryField()[patchi];
        pQHE = 0.0;
    }

    forAll(species(), specier)
    {
        if(specier != electronListPosition_) 
        { 
            const volScalarField& pDr = thermo_.composition().pD(specier);
            const scalarField& pDrCells = pDr.internalField();
            
            if(speciesThermo_[specier].particleType() < 3)
            {
                const scalar sigma_er = 1.0e-20;
                
                forAll(pDrCells, celli)
                {        
                    QHECells[celli] += sigma_er*pDr[celli]/sqr(W(electronListPosition_));
                }
                
                forAll(pDr.boundaryField(), patchi)
                {
                    const fvPatchScalarField& ppDr = pDr.boundaryField()[patchi];
                    fvPatchScalarField& pQHE = this->QHE_.boundaryField()[patchi];
                    
                    forAll(ppDr, facei)
                    {
                        pQHE[facei] += sigma_er*ppDr[facei]/sqr(W(electronListPosition_));
                    }
                }
            }
            else
            {
                forAll(pDrCells, celli)
                {        
                    scalar sigma_eIon = 8.0*pi*pow4(ec)/(27.0*sqr(kB*TeCells[celli]));
                    
                    if(nDeCells[celli] != 0.0)
                    {
                        sigma_eIon *= log(1.0+(9.0*pow3(kB*TeCells[celli]))/(4.0*pi*nDeCells[celli]*pow6(ec)));
                    }
                    
                    QHECells[celli] += sigma_eIon*pDr[celli]/sqr(W(electronListPosition_));
                }
                
                forAll(pDr.boundaryField(), patchi)
                {
                    const fvPatchScalarField& pnDe = nDe.boundaryField()[patchi];
                    const fvPatchScalarField& pTe = Te.boundaryField()[patchi];
                    const fvPatchScalarField& ppDr = pDr.boundaryField()[patchi];
                    fvPatchScalarField& pQHE = this->QHE_.boundaryField()[patchi];
                    
                    forAll(ppDr, facei)
                    {
                        scalar sigma_eIon = 8.0*pi*pow4(ec)/(27.0*sqr(kB*pTe[facei]));
                    
                        if(pnDe[facei] != 0.0)
                        {
                            sigma_eIon *= log(1.0+(9.0*pow3(kB*pTe[facei]))/(4.0*pi*pnDe[facei]*pow6(ec)));
                        }
                        
                        pQHE[facei] += sigma_eIon*ppDr[facei]/sqr(W(electronListPosition_));
                    }
                }
            }
        } 
    }//end species loop

    forAll(TtCells, celli)
    {        
        QHECells[celli] *= 3.0*RR*pDeCells[celli]*(TtCells[celli]-TeCells[celli])*NA
            *sqrt(8.0*RR*TeCells[celli]/(pi*W(electronListPosition_)));
      //Info<< "QHECells[celli]:"<< tab << QHECells[celli]<<endl; 
    }
    
    forAll(Tt.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
        const fvPatchScalarField& pTe = Te.boundaryField()[patchi];
        const fvPatchScalarField& ppDe = pDe.boundaryField()[patchi];
        fvPatchScalarField& pQHE = this->QHE_.boundaryField()[patchi];
        
        forAll(pTt, facei)
        {        
            pQHE[facei] *= 3.0*RR*ppDe[facei]*(pTt[facei]-pTe[facei])*NA
                *sqrt(8.0*RR*pTe[facei]/(pi*W(electronListPosition_)));
        }
    }
}  
    

template<class ThermoType>
bool Foam::AppletonBray<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}
   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
