/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "chemistry2Model.H"
#include "reacting2Mixture.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistry2Model<CompType, ThermoType>::chemistry2Model
(
    const fvMesh& mesh
)
:
    CompType(mesh),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reacting2Mixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reacting2Mixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),

    RR_(nSpecie_),
    RRf_(nSpecie_), // NEW VINCENT 25/03/2016
    
    preferentialFactor_(nSpecie_),
    simpleHarmonicOscillatorVibCutOff_(this->thermo().composition().solvedVibEqSpecies().size())
{
    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "RR_" + Y_[fieldI].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }
    
    // NEW VINCENT 25/03/2016 *************************************************
    // create the fields for the forward chemistry sources
    forAll(RRf_, fieldI)
    {
        RRf_.set
        (
            fieldI,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "RRfwd_" + Y_[fieldI].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }
    
    reactionRateNplusiir_ = 0.0; reactionRateOplusiir_ = 0.0;
    
    if(this->CVModel_ == "ParkTTv")
    {
        FixedList<scalar, 2> initList(0.0);
        if(this->preferentialModel_ == "constant")
        {
            initList[1] = readScalar(this->subDict("chemistryVibrationCoupling").subDict(this->CVModel_ + "Coeffs")
                .subDict("preferentialModel").lookup("constantFactor"));
        }    
        forAll(preferentialFactor_, speciei)
        {
            if(this->preferentialModel_ == "constant")
            {
                preferentialFactor_.set
                (
                    speciei, new FixedList<scalar, 2>(initList)
                );
            }
            else if(this->preferentialModel_ == "lineFitted")
            {
                preferentialFactor_.set
                (
                    speciei, new FixedList<scalar, 2>(this->subDict("chemistryVibrationCoupling").subDict(this->CVModel_ + "Coeffs").subDict("preferentialModel")
                        .subDict("lineFittedCoeffs").subDict("vibFactor").lookupOrDefault(Y_[speciei].name(), initList))
                );
            }
        }
    }
    else if(this->CVModel_ == "CVDV")
    {
        reciprocalUFactor_ = 1.0/readScalar(this->subDict("chemistryVibrationCoupling").subDict(this->CVModel_ + "Coeffs").lookup("reciprocalU"));
        
        forAll(simpleHarmonicOscillatorVibCutOff_, speciei) 
        {   
            simpleHarmonicOscillatorVibCutOff_.set
            (
                speciei, 
                new label
                (
                    this->subDict("chemistryVibrationCoupling").subDict(this->CVModel_ + "Coeffs")
                    .subDict("simpleHarmonicOscillatorVibCutOff").lookupOrDefault
                    (
                        this->thermo().composition().solvedVibEqSpecies()[speciei], 
                        1 + int(specieThermo_[speciei].dissociationPotential()/(specieThermo_[speciei].R()*specieThermo_[speciei].vibrationalList()[1])) // Thivet 1991, p. 2802
                    )
                )
            );
        }    
    }
    // END NEW VINCENT 25/03/2016 *********************************************   
    
    Info<< "chemistry2Model: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistry2Model<CompType, ThermoType>::~chemistry2Model()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::chemistry2Model<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar Tv, // NEW VINCENT 
    const scalar p
) const
{
    //Info<< "info statement in Foam::chemistry2Model<CompType, ThermoType>::omega(const scalarField& c,const scalar T,const scalar p) " << endl; //PRINTED = FALSE
    
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    tmp<scalarField> tom(new scalarField(nEqns(), 0.0));
    scalarField& om = tom();

    forAll(reactions_, i)
    {
        const Reaction2<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, Tv, p, pf, cf, lRef, pr, cr, rRef
        ); // MODIFIED VINCENT

        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return tom;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar Tv, // NEW VINCENT
    const List<scalar>& spTv, // NEW VINCENT 
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    // Info<< "info statement in Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::omegaI" << endl; //PRINTED = YES
    
    const Reaction2<ThermoType>& R = reactions_[index];
    scalar w = omega(R, c, T, Tv, spTv, p, pf, cf, lRef, pr, cr, rRef); // MODIFIED VINCENT
    return(w);
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::omega // NB VINCENT: overloaded for Tv[speciei]
(
    const Reaction2<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar Tv, // NEW VINCENT
    const List<scalar>& spTv, // NEW VINCENT 
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    // Info<< "info statement in overloaded Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::omega(const Reaction2<ThermoType>& R,const scalarField& c, ...)" << endl; PRINTED = YES
    
    scalarField c2(nSpecie_, 0.0);
    for (label i = 0; i < nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    // NEW VINCENT ******************************************************
    scalar kf = 0; // NEW VINCENT 07/02/2016
    scalar kr = 0; // NEW VINCENT 07/02/2016
    
    if (this->CVModel_ == "ParkTTv")
    {
        const Foam::ControllingTemperatureType controllingTemperature = R.controlT();
        scalar colTf = 0.0, colTr = 0.0;
        
        // The first species to appear on the left-hand side of the equation is the one chosen for spTv[]
        // for a dissociation reaction.
        const label labelSp = R.lhs()[0].index;
        
        switch(controllingTemperature)
        {
            case chargeExchange:          colTf = T;                colTr = T;                         break;
            case dissociation:            colTf = pow(T, this->exponentPark_)*pow(spTv[labelSp], 1.0-this->exponentPark_); colTr = T; break;
            case exchange:                colTf = T;                colTr = T;                         break;
            case impactDissociationOfN2:  colTf = T;                colTr = spTv[labelSp];             break;
            case impactIonisation:        colTf = spTv[R.rhs()[1].index];       colTr = spTv[R.rhs()[1].index];    break;
            case associativeIonisation:   colTf = T;                colTr = spTv[R.rhs()[0].index];    break;
            case transrotational:         colTf = T;                colTr = T;                         break;
            case vibrational:             colTf = spTv[labelSp];    colTr = spTv[labelSp];             break;
        }  
        
        kf = R.kf(p, colTf, c2);
        kr = R.kr(kf, p, colTr, c2);
    }
    else if (this->CVModel_ == "CVDV")
    {
        scalar colTf = T, colTr = T;
        
        // The first species to appear on the left-hand side of the equation is the one that undergoes dissociation
        const label labelSp = R.lhs()[0].index;

        const scalar reciprocalU = reciprocalUFactor_*specieThermo_[labelSp].R()/specieThermo_[labelSp].dissociationPotential();
        const scalar reciprocalTf = 1.0/spTv[labelSp] - 1.0/T - reciprocalU;
        
        scalar sumEnergyLevelsT = 0.0, sumEnergyLevelsTf = 0.0, sumEnergyLevelsTv = 0.0, sumEnergyLevelsU = 0.0;
        const scalar charVibTemp = specieThermo_[labelSp].vibrationalList()[1];
        
        for(label level=0 ; level<simpleHarmonicOscillatorVibCutOff_[labelSp]; level++)
        {
            sumEnergyLevelsT  += exp(-level*charVibTemp/T);
            sumEnergyLevelsTf += exp(-level*charVibTemp*reciprocalTf);
            sumEnergyLevelsTv += exp(-level*charVibTemp/spTv[labelSp]);
            sumEnergyLevelsU  += exp(level*charVibTemp*reciprocalU);
        }
        
        scalar CVDVCoeff = sumEnergyLevelsT*sumEnergyLevelsTf/(sumEnergyLevelsTv*sumEnergyLevelsU);
      //scalar CVDVCoeff = termQ(labelSp, T)*termQ(labelSp, 1.0/reciprocalTf)/(termQ(labelSp, spTv[labelSp])*termQ(labelSp, -1.0/reciprocalU)); TOO SLOW

        kf = R.kf(p, colTf, c2)*CVDVCoeff;
        kr = R.kr(kf, p, colTr, c2);
    }
    else
    {
        // Standard single-temperature formulation
        kf = R.kf(p, T, c2);
        kr = R.kr(kf, p, T, c2);
    }     
    // END NEW VINCENT **************************************************
    
    //const scalar kf = R.kf(p, T, c2); // DELETED VINCENT
    //const scalar kr = R.kr(kf, p, T, c2); // DELETED VINCENT

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }
    
    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::omega
(
    const Reaction2<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar Tv, // NEW VINCENT
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    //Info<< "info statement in Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::omega(const Reaction2<ThermoType>& R,const scalarField& c, ...)" << endl; //PRINTED = NO
    
    scalarField c2(nSpecie_, 0.0);
    for (label i = 0; i < nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    // NEW VINCENT ******************************************************
    scalar kf = 0; // NEW VINCENT 07/02/2016
    scalar kr = 0; // NEW VINCENT 07/02/2016
    
    if (this->CVModel_ == "ParkTTv")
    {
        const Foam::ControllingTemperatureType controllingTemperature = R.controlT();
        scalar colTf = 0.0, colTr = 0.0;
        
        switch(controllingTemperature)
        {
            case chargeExchange:      colTf = T;     colTr = T;  break;
            case dissociation:        colTf = pow(T, this->exponentPark_)*pow(Tv, 1.0-this->exponentPark_); colTr = T; break;
            case exchange:            colTf = T;     colTr = T;  break;
            case impactDissociationOfN2:  colTf = T; colTr = Tv; break;
            case impactIonisation:    colTf = Tv;    colTr = Tv; break;
            case associativeIonisation: colTf = T; colTr = Tv; break;
            case transrotational:     colTf = T;     colTr = T;  break;
            case vibrational:         colTf = Tv;    colTr = Tv; break;
        }  
        
        kf = R.kf(p, colTf, c2);
        kr = R.kr(kf, p, colTr, c2); 
    }
    else if (this->CVModel_ == "CVDV")
    {
        scalar colTf = T, colTr = T;
        
        const label labelSp = R.lhs()[0].index; // The first species to appear on the left-hand side of the equation is the one that undergoes dissociation
        
        const scalar reciprocalU = reciprocalUFactor_*specieThermo_[labelSp].R()/specieThermo_[labelSp].dissociationPotential();
        const scalar reciprocalTf = 1.0/Tv - 1.0/T - reciprocalU;
        scalar CVDVCoeff = termQ(labelSp, T)*termQ(labelSp, 1.0/reciprocalTf)/(termQ(labelSp, Tv)*termQ(labelSp, -1.0/reciprocalU));
        
        kf = R.kf(p, colTf, c2)*CVDVCoeff;
        kr = R.kr(kf, p, colTr, c2);
    }
    else
    {
        // Standard single-temperature formulation
        kf = R.kf(p, T, c2);
        kr = R.kr(kf, p, T, c2);
    }    
    // END BRAND NEW VINCENT **************************************************
    
    //const scalar kf = R.kf(p, T, c2); // DELETED VINCENT
    //const scalar kr = R.kr(kf, p, T, c2); // DELETED VINCENT

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }
    
    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::chemistry2Model<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    //Info<< "info statement in void Foam::chemistry2Model<CompType, ThermoType>::derivatives" << endl; // PRINTED = FALSE
    
    const scalar T = c[nSpecie_];
    const scalar Tv = c[nSpecie_ + 1]; // NEW VINCENT
    const scalar p = c[nSpecie_ + 2]; // MODIFIED VINCENT

    dcdt = omega(c, T, Tv, p); // MODIFIED VINCENT

    // constant pressure
    // dTtr/dt = ...
    scalar rho = 0.0;
    scalar cSum = 0.0;
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar W = specieThermo_[i].W();
        cSum += c[i];
        rho += W*c[i];
    }
    scalar cp = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        //cp += c[i]*specieThermo_[i].cp(p, T); // DELETED VINCENT
        cp += c[i]*specieThermo_[i].cp(p, T, 0); //TODO NEW VINCENT aiming for cp_t(p, Tt)
    }
    cp /= rho;

    scalar dT = 0.0;
    for (label i = 0; i < nSpecie_; i++)
    {
        //const scalar hi = specieThermo_[i].ha(p, T, T); //TODO il va falloir un dTtr et un dTv
        const scalar hi = specieThermo_[i].het(p, T); // NEW VINCENT
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[nSpecie_] = -dT;

    // dTv/dt = ...
    //TODO
    dcdt[nSpecie_ + 1] = 0.0; // TODO NEW VINCENT
    
    // dp/dt = ...
    dcdt[nSpecie_ + 2] = 0.0; // MODIFIED VINCENT
}


template<class CompType, class ThermoType>
void Foam::chemistry2Model<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    //Info<< "info statement in void Foam::chemistry2Model<CompType, ThermoType>::jacobian" << endl; // PRINTED = FALSE
    
    const scalar T = c[nSpecie_];
    const scalar Tv = c[nSpecie_ + 1]; // NEW VINCENT
    const scalar p = c[nSpecie_ + 2]; // MODIFIED VINCENT

    scalarField c2(nSpecie_, 0.0);
    forAll(c2, i)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    // length of the first argument must be nSpecie()
    dcdt = omega(c2, T, Tv, p); // MODIFIED VINCENT

    forAll(reactions_, ri)
    {
        const Reaction2<ThermoType>& R = reactions_[ri];

        const scalar kf0 = R.kf(p, T, c2); //TODO introduce the switch
        const scalar kr0 = R.kr(p, T, c2); //TODO

        forAll(R.lhs(), j)
        {
            const label sj = R.lhs()[j].index;
            scalar kf = kf0;
            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar el = R.lhs()[i].exponent;
                if (i == j)
                {
                    if (el < 1.0)
                    {
                        if (c2[si] > SMALL)
                        {
                            kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        kf *= el*pow(c2[si], el - 1.0);
                    }
                }
                else
                {
                    kf *= pow(c2[si], el);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] -= sl*kf;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] += sr*kf;
            }
        }

        forAll(R.rhs(), j)
        {
            const label sj = R.rhs()[j].index;
            scalar kr = kr0;
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar er = R.rhs()[i].exponent;
                if (i == j)
                {
                    if (er < 1.0)
                    {
                        if (c2[si] > SMALL)
                        {
                            kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                        }
                        else
                        {
                            kr = 0.0;
                        }
                    }
                    else
                    {
                        kr *= er*pow(c2[si], er - 1.0);
                    }
                }
                else
                {
                    kr *= pow(c2[si], er);
                }
            }

            forAll(R.lhs(), i)
            {
                const label si = R.lhs()[i].index;
                const scalar sl = R.lhs()[i].stoichCoeff;
                dfdc[si][sj] += sl*kr;
            }
            forAll(R.rhs(), i)
            {
                const label si = R.rhs()[i].index;
                const scalar sr = R.rhs()[i].stoichCoeff;
                dfdc[si][sj] -= sr*kr;
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1.0e-3;
    const scalarField dcdT0(omega(c2, T - delta, Tv - delta, p)); // MODIFIED VINCENT
    const scalarField dcdT1(omega(c2, T + delta, Tv - delta, p)); // MODIFIED VINCENT

    for (label i = 0; i < nEqns(); i++)
    {
        dfdc[i][nSpecie()] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }
    
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistry2Model<CompType, ThermoType>::tc() const
{
    Info<< "info statement in Foam::chemistry2Model<CompType, ThermoType>::tc()" << endl;
    
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, SMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc();
    const scalarField& T = this->thermo().Tt(); // MODIFIED VINCENT
    const scalarField& Tv = this->thermo().Tv(); // NEW VINCENT
    const scalarField& p = this->thermo().p();

    const label nReaction = reactions_.size();

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            scalar rhoi = rho[celli];
            scalar Ti = T[celli];
            scalar Tvi = Tv[celli]; // NEW VINCENT
            scalar pi = p[celli];
            scalarField c(nSpecie_);
            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                scalar Yi = Y_[i][celli];
                c[i] = rhoi*Yi/specieThermo_[i].W();
                cSum += c[i];
            }

            forAll(reactions_, i)
            {
                const Reaction2<ThermoType>& R = reactions_[i];

                omega(R, c, Ti, Tvi, pi, pf, cf, lRef, pr, cr, rRef); // NEW VINCENT

                forAll(R.rhs(), s)
                {
                    scalar sr = R.rhs()[s].stoichCoeff;
                    tc[celli] += sr*pf*cf;
                }
            }
            tc[celli] = nReaction*cSum/tc[celli];
        }
    }


    ttc().correctBoundaryConditions();
    
    return ttc;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistry2Model<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    if (this->chemistry_)
    {
        scalarField& Sh = tSh();
        
        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                const scalar hi = specieThermo_[i].Hc();
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }

    return tSh;
}


// NEW VINCENT ****************************************************************
template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistry2Model<CompType, ThermoType>::Scv(label i) const
{
    tmp<volScalarField> tScv
    (
        new volScalarField
        (
            IOobject
            (
                "Scv_" + Y_[i].name(),
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Scv = tScv();
        
        const scalarField& TtCells = this->thermo().Tt();
        const scalarField& TvCells = this->thermo().composition().Tv()[i];
        const scalarField& pCells = this->thermo().p();
            
        if (this->CVModel_ == "ParkTTv")
        {
            if (this->ScvModel_ == "preferential")
            {
                const scalar dissociationPotential = specieThermo_[i].dissociationPotential();
                forAll(Scv, cellI)
                {
                    Scv[cellI] = ((preferentialFactor_[i][0]*TtCells[cellI]+preferentialFactor_[i][1])*dissociationPotential 
                        + specieThermo_[i].HEel(pCells[cellI], TvCells[cellI]))*RR_[i][cellI];
                }
            }
            else if (this->ScvModel_ == "nonPreferential")
            {
                forAll(Scv, cellI)
                {
                    const scalar Dprime = specieThermo_[i].HEvel(pCells[cellI], TvCells[cellI]);
                    Scv[cellI] = Dprime*RR_[i][cellI];
                }
            } 
        }
        else if (this->CVModel_ == "CVDV")
        {
            const scalar speciesR = specieThermo_[i].R();
            const scalar charVibTemp = specieThermo_[i].vibrationalList()[1];
            const label cutoffValue = simpleHarmonicOscillatorVibCutOff_[i];
            const scalar reciprocalU = reciprocalUFactor_*speciesR/specieThermo_[i].dissociationPotential();
            
            forAll(Scv, cellI)
            {
                scalar sumEnergyEnergyLevelsTf = 0.0, sumEnergyEnergyLevelsU = 0.0,
                    sumEnergyLevelsTf = 0.0, sumEnergyLevelsU = 0.0;
                
                scalar reciprocalTf = 1.0/TvCells[cellI] - 1.0/TtCells[cellI] - reciprocalU;
                
                for(label level=0 ; level<cutoffValue; level++)
                {
                    sumEnergyEnergyLevelsTf += (level*speciesR*charVibTemp)
                        * exp(-level*charVibTemp*reciprocalTf);
                        
                    sumEnergyEnergyLevelsU += (level*speciesR*charVibTemp)
                        * exp(level*charVibTemp*reciprocalU);            
                
                    sumEnergyLevelsTf += exp(-level*charVibTemp*reciprocalTf);
                    sumEnergyLevelsU += exp(level*charVibTemp*reciprocalU);
                }
            
                /*Scv[cellI] = sumEnergyEnergyLevelsTf/termQ(i, 1.0/reciprocalTf)*RRf_[i][cellI] 
                    + sumEnergyEnergyLevelsU/termQ(i, -1.0/reciprocalU)*(RR_[i][cellI]-RRf_[i][cellI]);*/
                    
                Scv[cellI] = sumEnergyEnergyLevelsTf/sumEnergyLevelsTf*RRf_[i][cellI] 
                    + sumEnergyEnergyLevelsU/sumEnergyLevelsU*(RR_[i][cellI]-RRf_[i][cellI]);
            }
        }
    }

    return tScv;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistry2Model<CompType, ThermoType>::Seiir(label i) const
{
    tmp<volScalarField> tSeiir
    (
        new volScalarField
        (
            IOobject
            (
                "Seiir_" + Y_[i].name(),
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Seiir = tSeiir();
        label NplusLabel = -1; label OplusLabel = -1;
        forAll(Y_, speciei)
        {
            if(specieThermo_[speciei].nameSpecies() == "N+")
            {
                NplusLabel = speciei;
            }
            if(specieThermo_[speciei].nameSpecies() == "O+")
            {
                OplusLabel = speciei;
            }
        }
         
        if(NplusLabel != -1 or OplusLabel != -1)
        {
            Seiir = 0.0;
            if(specieThermo_[i].nameSpecies() == "N2") // N + e- <-> N+ + 2e-
            {
                forAll(Seiir, celli)
                {
                    Seiir[celli] = specieThermo_[i].iHat()*reactionRateNplusiir_*specieThermo_[NplusLabel].W()*1.0e-3;
                }
            }
            else if(specieThermo_[i].nameSpecies() == "O") // O + e- <-> O+ + 2e-
            {
                forAll(Seiir, celli)
                {
                    Seiir[celli] = specieThermo_[i].iHat()*reactionRateOplusiir_*specieThermo_[OplusLabel].W()*1.0e-3;
                }
            }
        }
    }

    return tSeiir;
}
// END NEW VINCENT ************************************************************


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistry2Model<CompType, ThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class ThermoType>
Foam::label Foam::chemistry2Model<CompType, ThermoType>::nEqns() const
{
    // nEqns = number of species + temperature_tr + temperature_v + pressure
    return nSpecie_ + 2 + 1;  // MODIFIED VINCENT
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::chemistry2Model<CompType, ThermoType>::calculateRR
(
    const label reactionI,
    const label specieI
) const
{
    Info<< "info statement in Foam::chemistry2Model<CompType, ThermoType>::calculateRR(...)" << endl; //PRINTED = FALSE
    
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    tmp<DimensionedField<scalar, volMesh> > tRR
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& RR = tRR();

    const scalarField& T = this->thermo().Tt(); // MODIFIED VINCENT
    const scalarField& Tv = this->thermo().Tv(); // NEW VINCENT
    List<scalarField> spTv(nSpecie_); // NEW VINCENT
    forAll (spTv, speciei)
    { 
        spTv[speciei] = this->thermo().composition().Tv()[speciei];
    } // NEW VINCENT
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar Tvi = Tv[celli]; // NEW VINCENT
        List<scalar> spTvi(nSpecie_); // NEW VINCENT
        forAll (spTv, speciei)
        {
            spTvi[speciei] = spTv[speciei][celli];
        } // NEW VINCENT   
        const scalar pi = p[celli];

        scalarField c(nSpecie_, 0.0);
        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c[i] = rhoi*Yi/specieThermo_[i].W();
        }

        const scalar w = omegaI
        (
            reactionI,
            c,
            Ti,
            Tvi, // MODIFIED VINCENT
            spTvi, // NEW VINCENT
            pi,
            pf,
            cf,
            lRef,
            pr,
            cr,
            rRef
        );

        RR[celli] = w*specieThermo_[specieI].W();

    }

    return tRR;
}


template<class CompType, class ThermoType>
void Foam::chemistry2Model<CompType, ThermoType>::calculate()
{
    //Info<< "info statement in Foam::chemistry2Model<CompType, ThermoType>::calculate()" << endl; //PRINTED = FALSE
    
    if (!this->chemistry_)
    {
        return;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().Tt(); // MODIFIED VINCENT
    const scalarField& Tv = this->thermo().Tv(); // NEW VINCENT
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar Tvi = Tv[celli]; // NEW VINCENT
        const scalar pi = p[celli];
        scalarField c(nSpecie_, 0.0);
        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c[i] = rhoi*Yi/specieThermo_[i].W();
        }

        const scalarField dcdt(omega(c, Ti, Tvi, pi)); // MODIFIED VINCENT

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dcdt[i]*specieThermo_[i].W();
        }
    }
}


template<class CompType, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    //Info<< "info statement in Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::solve(const DeltaTType& deltaT)" << endl; // PRINTED = TRUE, WITH OR WOUT CHEMISTRY
    
    CompType::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().Tt(); // MODIFIED VINCENT
    const scalarField& Tv = this->thermo().Tv(); // NEW VINCENT
    List<scalarField> spTv(nSpecie_); // NEW VINCENT
    forAll (spTv, speciei)
    { 
        spTv[speciei] = this->thermo().composition().Tv()[speciei];
    } // NEW VINCENT
    const scalarField& p = this->thermo().p();

    scalarField c(nSpecie_);
    scalarField cfwd(nSpecie_); // NEW VINCENT 25/03/2016
    scalarField c0(nSpecie_);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];
        scalar Tvi = Tv[celli]; // NEW VINCENT
        List<scalar> spTvi(nSpecie_); // NEW VINCENT
        forAll (spTv, speciei)
        {
            spTvi[speciei] = spTv[speciei][celli];
        } // NEW VINCENT   

        for (label i=0; i<nSpecie_; i++)
        {
            c[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
            c0[i] = c[i];
        }

        // Initialise time progress
        scalar timeLeft = deltaT[celli];
        

        // Calculate the chemical source terms
        cfwd = c; // NEW VINCENT 25/03/2016
        while (timeLeft > SMALL)
        {
            scalar dt = timeLeft;           
            //this->solve(c, Ti, Tvi, spTvi, pi, dt, this->deltaTChem_[celli]); // MODIFIED VINCENT
            this->solve(c, cfwd, reactionRateNplusiir_, reactionRateOplusiir_, Ti, Tvi, spTvi, pi, dt, this->deltaTChem_[celli]); // NEW VINCENT 25/03/2016
            
            timeLeft -= dt;
        }

        deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = (c[i] - c0[i])*specieThermo_[i].W()/deltaT[celli];
            RRf_[i][celli] = (cfwd[i] - c0[i])*specieThermo_[i].W()/deltaT[celli]; // NEW VINCENT 25/03/2016
        }
    }

    //Info<< "deltaTMin : " << deltaTMin << endl;
    return deltaTMin;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::solve
(
    const scalar deltaT
)
{
    //Info<< "info statement in Foam::chemistry2Model<CompType, ThermoType>::solve(const scalar deltaT)" << endl; this function is used at first
    
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar> >(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class CompType, class ThermoType>
Foam::scalar Foam::chemistry2Model<CompType, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    Info<< "info statement in Foam::chemistry2Model<CompType, ThermoType>::solve(const scalarField& deltaT)" << endl; // PRINTED = FALSE
    return this->solve<scalarField>(deltaT);
}


template<class CompType, class ThermoType>
void Foam::chemistry2Model<CompType, ThermoType>::solve
(
    scalarField &c,
    scalar& T,
    scalar& Tv, // MODIFIED VINCENT
    List<scalar>& spTv, // NEW VINCENT
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    notImplemented
    (
        "chemistry2Model::solve"
        "("
            "scalarField&, "
            "scalar&, "
            "scalar&, "
            "scalar&, " // MODIFIED VINCENT
            "List<scalar>&, " // NEW VINCENT
            "scalar&, "
            "scalar&"
        ") const"
    );
}


template<class CompType, class ThermoType>
void Foam::chemistry2Model<CompType, ThermoType>::solve
(
    scalarField& c,
    scalarField& cfwd,
    scalar& rriirN,
    scalar& rriirO,
    scalar& T,
    scalar& Tv, // MODIFIED VINCENT
    List<scalar>& spTv, // NEW VINCENT
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    notImplemented
    (
        "chemistry2Model::solve"
        "("
            "scalarField&, "
            "scalarField&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&, " // MODIFIED VINCENT
            "List<scalar>&, " // NEW VINCENT
            "scalar&, "
            "scalar&"
        ") const"
    );
}


// ************************************************************************* //
