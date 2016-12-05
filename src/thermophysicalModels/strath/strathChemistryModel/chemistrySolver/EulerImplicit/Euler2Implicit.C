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

#include "Euler2Implicit.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"

#include "rho2ReactionThermo.H" // NEW VINCENT 21/03/2016

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Chemistry2Model>
Foam::Euler2Implicit<Chemistry2Model>::Euler2Implicit
(
    const fvMesh& mesh
)
:
    chemistry2Solver<Chemistry2Model>(mesh),
    coeffsDict_(this->subDict("Euler2ImplicitCoeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter")),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Chemistry2Model>
Foam::Euler2Implicit<Chemistry2Model>::~Euler2Implicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Chemistry2Model>
void Foam::Euler2Implicit<Chemistry2Model>::updateRRInReactionI
(
    const label index,
    const scalar pr,
    const scalar pf,
    const scalar corr,
    const label lRef,
    const label rRef,
    const scalar p,
    const scalar T,
    const scalar Tv, // MODIFIED VINCENT
    simpleMatrix<scalar>& RR
) const
{
    // Info<< "info statement in void Foam::Euler2Implicit<Chemistry2Model>::updateRRInReactionI" << endl; // PRINTED = YES
    // CALLED FROM FUNCTION void Foam::Euler2Implicit<Chemistry2Model>::solve DOWN BELOW
    
    //Info<< "prev["<<index<<"]: " << pr << endl;
    //Info<< "pfwd["<<index<<"]: " << pf << endl;
    
    const Reaction2<typename Chemistry2Model::thermoType>& R =
        this->reactions_[index];

    forAll(R.lhs(), s)
    {
        label si = R.lhs()[s].index;
        scalar sl = R.lhs()[s].stoichCoeff;
        RR[si][rRef] -= sl*pr*corr;
        RR[si][lRef] += sl*pf*corr;
    }

    forAll(R.rhs(), s)
    {
        label si = R.rhs()[s].index;
        scalar sr = R.rhs()[s].stoichCoeff;
        RR[si][lRef] -= sr*pf*corr;
        RR[si][rRef] += sr*pr*corr;
    }
    
    // NEW VINCENT 29/03/2016 *************************************************
    /*if(R.controlT() == impactIonisation)
    {
        if(R.species().contains("N+"))
        {
            //this->setReactionRateNplusiir_(0.0);
        }
        else if(R.species().contains("O+"))
        {
            //this->setReactionRateOplusiir_(0.0);
        }
    }*/
    // END NEW VINCENT 29/03/2016 *********************************************
    
}


template<class Chemistry2Model>
void Foam::Euler2Implicit<Chemistry2Model>::solve
(
    scalarField& c,
    scalar& T,
    scalar& Tv, // NEW VINCENT
    List<scalar>& spTv, // NEW VINCENT
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    //Info<< "info statement in void Foam::Euler2Implicit<Chemistry2Model>::solve" << endl; // PRINTED = NO
    
    const label nSpecie = this->nSpecie();
    simpleMatrix<scalar> RR(nSpecie, 0, 0);

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }
    // Calculate the absolute enthalpy // DELETED VINCENT
    // Calculate the trans-rotational enthalpy // NEW VINCENT
    // Calculate the vibrational enthalpy // NEW VINCENT
    scalar cTot = sum(c);
    
    typename Chemistry2Model::thermoType mixture
    (
        (c[0]/cTot)*this->specieThermo_[0]
    );
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i];
    }
    // scalar ha = mixture.Ha(p, T); // DELETED VINCENT
    // NEW VINCENT ************************************************************
    scalar ht = mixture.HEt(p, T);
    List<scalar> sphv(nSpecie);
    forAll(sphv, speciei)
    {
        sphv[speciei] = this->specieThermo_[speciei].HEvel(p, spTv[speciei]);
    }
    // END NEW VINCENT ********************************************************

    scalar deltaTEst = min(deltaT, subDeltaT);

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pr, cr;
        label lRef, rRef;

        scalar omegai = this->omegaI(i, c, T, Tv, spTv, p, pf, cf, lRef, pr, cr, rRef); //MODIFIED VINCENT

        scalar corr = 1.0;
        if (eqRateLimiter_)
        {
            if (omegai < 0.0)
            {
                corr = 1.0/(1.0 + pr*deltaTEst);
            }
            else
            {
                corr = 1.0/(1.0 + pf*deltaTEst);
            }
        }

        updateRRInReactionI(i, pr, pf, corr, lRef, rRef, p, T, Tv, RR); // MODIFIED VINCENT (Tv)
    }
    
    // Calculate the stable/accurate time-step
    scalar tMin = GREAT;

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = 0;
        for (label j=0; j<nSpecie; j++)
        {
            d -= RR[i][j]*c[j];
        }

        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(cTot - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }

    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);

    // Add the diagonal and source contributions from the time-derivative
    for (label i=0; i<nSpecie; i++)
    {
        RR[i][i] += 1.0/deltaT;
        RR.source()[i] = c[i]/deltaT;
    }

    // Solve for the new composition
    c = RR.LUsolve();

    // Limit the composition
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    // Update the temperatures
    cTot = sum(c);
    
    mixture = (c[0]/cTot)*this->specieThermo_[0];
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i];
    }
    // T = mixture.THa(ha, p, T); // DELETED VINCENT
    // NEW VINCENT ************************************************************
    T = mixture.TtHEt(ht, p, T);
    
    scalar vibCutOffTemp = Foam::rho2ReactionThermo::vibrationalCutOffTemp;

    forAll(sphv, speciei)
    {
        if (this->specieThermo_[speciei].particleType() == 2) 
        {
            //if(composition().noVibrationalTemp(speciei) == 1 or downgradeSingleVibMode) TODO
            //{
                if(T > vibCutOffTemp)
                {
                    spTv[speciei] = this->specieThermo_[speciei].TvelHEvel(sphv[speciei], p, spTv[speciei]);
                }
                else
                {
                    spTv[speciei] = T;
                }
           // } 
        }
        else if (sphv[speciei] != 0.0)
        {
            spTv[speciei] = spTv[this->thermo().composition().vibTempAssociativity(speciei)]; // condition fulfilled for ions, electrons, and atoms with Eel on
        }
    }
    // END NEW VINCENT ********************************************************
    
    /* // This was commented originally 
    for (label i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = p;
    */ 
}


template<class Chemistry2Model>
void Foam::Euler2Implicit<Chemistry2Model>::solve
(
    scalarField& c,
    scalarField& cfwd,
    scalar& rriirN, // NEW VINCENT
    scalar& rriirO, // NEW VINCENT
    scalar& T,
    scalar& Tv, // NEW VINCENT
    List<scalar>& spTv, // NEW VINCENT
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    //Info<< "info statement in void Foam::Euler2Implicit<Chemistry2Model>::solve" << endl; // PRINTED = YES
    
    const label nSpecie = this->nSpecie();
    simpleMatrix<scalar> RR(nSpecie, 0, 0);
    simpleMatrix<scalar> RRfwd(nSpecie, 0, 0); // NEW VINCENT 25/03/2016

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }
    // Calculate the absolute enthalpy // DELETED VINCENT
    // Calculate the trans-rotational enthalpy // NEW VINCENT
    // Calculate the vibrational enthalpy // NEW VINCENT
    scalar cTot = sum(c);
    /*typename Chemistry2Model::thermoType mixture
    (
        (c[0]/cTot)*this->specieThermo_[0]
    );
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i]; // NOTE VINCENT: c[i]/cTot = Xi
    }
    // scalar ha = mixture.Ha(p, T); // DELETED VINCENT
    // NEW VINCENT ************************************************************
    Info << "1) Tt: " << T << " and ht: " << mixture.HEt(p, T) << endl;
    scalar ht = mixture.HEt(p, T);
    List<scalar> sphv(nSpecie);
    forAll(sphv, speciei)
    {
        //ht += (c[speciei]/cTot)*this->specieThermo_[speciei].HEt(p, T);
        sphv[speciei] = this->specieThermo_[speciei].HEvel(p, spTv[speciei]);
    }
    Info << "2) Tt: " << T << " and ht: " << mixture.HEt(p, T) << endl;*/
    // END NEW VINCENT ********************************************************

    scalar deltaTEst = min(deltaT, subDeltaT);

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pr, cr;
        label lRef, rRef;

        scalar omegai = this->omegaI(i, c, T, Tv, spTv, p, pf, cf, lRef, pr, cr, rRef); //MODIFIED VINCENT

        scalar corr = 1.0;
        if (eqRateLimiter_)
        {
            if (omegai < 0.0)
            {
                corr = 1.0/(1.0 + pr*deltaTEst);
            }
            else
            {
                corr = 1.0/(1.0 + pf*deltaTEst);
            }
        }

        updateRRInReactionI(i, pr, pf, corr, lRef, rRef, p, T, Tv, RR); // MODIFIED VINCENT (Tv)
        updateRRInReactionI(i, 0, pf, corr, lRef, rRef, p, T, Tv, RRfwd); // NEW VINCENT 25/03/2016
        
        // NEW VINCENT 30/03/2016 *************************************************
        if(this->reactions()[i].controlT() == impactIonisation)
        {
            if(this->reactions()[i].species().contains("N+"))
            {
                rriirN = RR[this->reactions()[i].rhs()[0].index][rRef] + RR[this->reactions()[i].rhs()[0].index][lRef];
            }
            else if(this->reactions()[i].species().contains("O+"))
            {
                rriirO = RR[this->reactions()[i].rhs()[0].index][rRef] + RR[this->reactions()[i].rhs()[0].index][lRef];
            }
            //Info << "There is at least one impactIonisation reaction considered." << endl;
            //Info << "reac_rate iir N: " << rriirN << endl;
            //Info << "reac_rate iir O: " << rriirO << endl;
        }
        // END NEW VINCENT 30/03/2016 *********************************************
    }
    
    // Calculate the stable/accurate time-step
    scalar tMin = GREAT;

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = 0;
        for (label j=0; j<nSpecie; j++)
        {
            d -= RR[i][j]*c[j];
        }

        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(cTot - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }

    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);

    // Add the diagonal and source contributions from the time-derivative
    for (label i=0; i<nSpecie; i++)
    {
        RR[i][i] += 1.0/deltaT;
        RR.source()[i] = c[i]/deltaT;
        
        RRfwd[i][i] += 1.0/deltaT; // NEW VINCENT 25/03/2016
        RRfwd.source()[i] = cfwd[i]/deltaT; // NEW VINCENT 25/03/2016
    }

    // Solve for the new composition
    c = RR.LUsolve();
    cfwd = RRfwd.LUsolve(); // NEW VINCENT 25/03/2016

    // Limit the composition
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
        cfwd[i] = max(0.0, cfwd[i]); // NEW VINCENT 25/03/2016
    }

    // Update the temperatures
    cTot = sum(c);
    /*Info << "3) Tt: " << T << " and ht: " << mixture.HEt(p, T) << endl;
    mixture = (c[0]/cTot)*this->specieThermo_[0];
    for (label i=1; i<nSpecie; i++)
    {
        mixture += (c[i]/cTot)*this->specieThermo_[i];
    }
    Info << "4) Tt: " << T << " and ht: " << mixture.HEt(p, T) << endl;
    // T = mixture.THa(ha, p, T); // DELETED VINCENT
    // NEW VINCENT ************************************************************
    T = mixture.TtHEt(ht, p, T);
    Info << "5) Tt: " << T << " and ht: " << mixture.HEt(p, T) << endl;
    
    scalar vibCutOffTemp = Foam::rho2ReactionThermo::vibrationalCutOffTemp;

    forAll(sphv, speciei)
    {
        if (this->specieThermo_[speciei].particleType() == 2) 
        {
            //if(composition().noVibrationalTemp(speciei) == 1 or downgradeSingleVibMode)
            //{
                if(T > vibCutOffTemp)
                {
                    spTv[speciei] = this->specieThermo_[speciei].TvelHEvel(sphv[speciei], p, spTv[speciei]);
                }
                else
                {
                    spTv[speciei] = T;
                }
           // } 
        }
        else if (sphv[speciei] != 0.0)
        {
            spTv[speciei] = spTv[this->thermo().composition().vibTempAssociativity(speciei)]; // condition fulfilled for ions, electrons, and atoms with Eel on
        }
    }*/
    // END NEW VINCENT ********************************************************
    
    /* // This was commented originally 
    for (label i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = p;
    */ 
}
// ************************************************************************* //
