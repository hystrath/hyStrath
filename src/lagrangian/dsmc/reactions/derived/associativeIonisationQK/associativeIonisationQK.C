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

Description

\*---------------------------------------------------------------------------*/

#include "associativeIonisationQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(associativeIonisationQK, 0);

addToRunTimeSelectionTable(dsmcReaction, associativeIonisationQK, dictionary);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void associativeIonisationQK::setProperties()
{
    dsmcReaction::setProperties();

    if (reactantIds_.size() != 2)
    {
        //- There must be exactly 2 reactants
        FatalErrorIn("associativeIonisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two reactants, instead of "
            << reactantIds_.size() << nl
            << exit(FatalError);
    }

    //- Reading in the intermediate molecule
    const word intermediateMolecule =
        propsDict_.lookupOrDefault<word>
        (
            "intermediateMolecule", word::null
        );

    if (intermediateMolecule == word::null)
    {
        FatalErrorIn("associativeIonisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "The intermediate molecule is not specified" << nl
            << exit(FatalError);
    }

    intermediateMoleculeId_ =
        findIndex(cloud_.typeIdList(), intermediateMolecule);

    //- Check that reactants belong to the typeIdList as defined in
    //  constant/dsmcProperties
    if (intermediateMoleculeId_ == -1)
    {
        FatalErrorIn("dsmcReaction::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "Cannot find type id of the intermediate molecule: "
            << intermediateMoleculeId_ << nl
            << exit(FatalError);
    }

    //- Reading in associative ionisation products
    const wordList productsAssociativeIonisation
    (
        propsDict_.lookup("associativeIonisationProducts")
    );

    if (productsAssociativeIonisation.size() != 2)
    {
        FatalErrorIn("associativeIonisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two products, instead of "
            << productsAssociativeIonisation.size() << nl
            << exit(FatalError);
    }

    productIdsAssociativeIonisation_.setSize
    (
        productsAssociativeIonisation.size()
    );

    labelList productsIds(label(2), -1);

    forAll(productsAssociativeIonisation, p)
    {
        productsIds[p] =
            findIndex
            (
                cloud_.typeIdList(),
                productsAssociativeIonisation[p]
            );

        if (productsIds[p] == -1)
        {
            FatalErrorIn("associativeIonisationQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Cannot find type id of product: "
                << productsAssociativeIonisation[p] << nl
                << exit(FatalError);
        }
    }


    //- Check that reactants belong to the typeIdList as defined in
    //  constant/dsmcProperties
    if (intermediateMoleculeId_ == -1)
    {
        FatalErrorIn("associativeIonisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "Cannot find type id of the intermediate molecule: "
            << intermediateMoleculeId_ << nl
            << exit(FatalError);
    }

    heatOfReactionAssociativeIonisationJoules_.setSize(2, 0.0);

    if (reactantTypes_[0] == 10 and reactantTypes_[1] == 10)
    {
        //- This reation is a forward associative ionisation reaction
        forwardAssociativeIonisation_ = true;

        //- The two atoms recombine to give the intermediate molecule
        //  with a heat of reaction of opposite sign that of dissociation
        heatOfReactionAssociativeIonisationJoules_[0] =
            -physicoChemical::k.value()
           *cloud_.constProps(intermediateMoleculeId_).thetaD();

        //- The intermediate molecule ionises
        heatOfReactionAssociativeIonisationJoules_[1] =
            physicoChemical::k.value()
           *cloud_.constProps(intermediateMoleculeId_).ionisationTemperature();

        //- Check that one product is an ionised molecule and one is an
        //  electron
        if
        (
            cloud_.constProps(productsIds[0]).charge()
           *cloud_.constProps(productsIds[1]).charge() != -1
        )
        {
            FatalErrorIn("associativeIonisationQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "For forward associative ionisation reactions, one product "
                << "must be an ionised molecule and the second an electron"
                << nl << exit(FatalError);
        }
        else
        {
            if
            (
                cloud_.constProps(productsIds[0]).charge()
             == 1
            )
            {
                posIonisedMol_ = 0;
                //- The first product listed is the ionised molecule
                productIdsAssociativeIonisation_[0] = productsIds[0];
                productIdsAssociativeIonisation_[1] = productsIds[1];

            }
            else
            {
                //- The first product listed is the ionised molecule
                productIdsAssociativeIonisation_[0] = productsIds[1];
                productIdsAssociativeIonisation_[1] = productsIds[0];
            }
        }

    }
    else if ((reactantTypes_[0] == 10) xor (reactantTypes_[1] == 10))
    {
        //- If one reactant is an atom, then the second must be an atom
        FatalErrorIn("associativeIonisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "One reactant is a neutral atom and the second is not." << nl
            << "For forward associative ionisation reactions, both reactants "
            << "are atoms" << nl
            << exit(FatalError);
    }
    else if (reactantTypes_[0] == 0 or reactantTypes_[0] == 21)
    {
        if
        (
            (reactantTypes_[0] == 0  and reactantTypes_[1] != 21)
         or (reactantTypes_[0] == 21 and reactantTypes_[1] != 0)
        )
        {
            //- If one reactant is an electron, then the second must be an
            //  ionised molecule
            FatalErrorIn("associativeIonisationQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "For reverse associative ionisation reactions, one reactant "
                << "must be an electron and the second an ionised molecule"
                << nl << exit(FatalError);
        }
        else
        {
            posIonisedMol_ = 0;
            if (reactantTypes_[1] == 21)
            {
                posIonisedMol_ = 1;
            }

            //- The ionised molecule and electron recombine to give the
            //  intermediate molecule with a heat of reaction of opposite sign
            //  that of ionisation
            heatOfReactionAssociativeIonisationJoules_[0] =
                -physicoChemical::k.value()
               *cloud_.constProps(intermediateMoleculeId_).ionisationTemperature();

            //- The intermediate molecule dissociates
            heatOfReactionAssociativeIonisationJoules_[1] =
                physicoChemical::k.value()
               *cloud_.constProps(intermediateMoleculeId_).thetaD();

            //- Check that both products are neutral atoms
            if
            (
                cloud_.constProps(productsIds[0]).type() != 10
             or cloud_.constProps(productsIds[1]).type() != 10
            )
            {
                FatalErrorIn("associativeIonisationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "For reverse associative ionisation reactions, both "
                    << "products must be neutral atoms"
                    << exit(FatalError);
            }
            else
            {
                productIdsAssociativeIonisation_[0] = productsIds[0];
                productIdsAssociativeIonisation_[1] = productsIds[1];
            }
        }
    }
    else
    {
        FatalErrorIn("associativeIonisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "For forward associative ionisation reactions, both reactants "
            << "are atoms" << nl
            << "For reverse associative ionisation reactions, one reactant "
            << "must be an electron and the second an ionised molecule" << nl
            << "In the present configuration, none of these conditions are "
            << "fulfilled" << nl
            << exit(FatalError);
    }
}


void associativeIonisationQK::testForwardAssociativeIonisation
(
    const dsmcParcel& p,
    const dsmcParcel& q,
    const scalar translationalEnergy,
    scalar& collisionEnergy,
    scalar& totalReactionProbability,
    scalar& reactionProbability
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    const scalar EEleP =
        cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
    const scalar EEleQ =
        cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

    const scalar omegaIntermediate = cloud_.constProps(intermediateMoleculeId_).omega();
    const scalar thetaVIntermediate = cloud_.constProps(intermediateMoleculeId_).thetaV_m(0);
    const scalarList& EElistIntermediate = cloud_.constProps(intermediateMoleculeId_).electronicEnergyList();
    const labelList& gListIntermediate = cloud_.constProps(intermediateMoleculeId_).electronicDegeneracyList();

    //- Collision energy is the sum of the relative kinetic energy of the two
    //  atoms, plus their electronic energies
    scalar Ec = translationalEnergy + EEleP + EEleQ;

    //- Trial L-B redistribution
    label iMax = Ec/(physicoChemical::k.value()*thetaVIntermediate);

    if (iMax > 0)
    {
        const label trialVibLevelIntermediate =
            cloud_.postCollisionVibrationalEnergyLevel
            (
                true,
                0.0,
                iMax,
                thetaVIntermediate,
                cloud_.constProps(intermediateMoleculeId_).thetaD(),
                cloud_.constProps(intermediateMoleculeId_).TrefZv()[0],
                omegaIntermediate,
                cloud_.constProps(intermediateMoleculeId_).Zref()[0],
                Ec // TODO VINCENT add last optional params to all postCollVib
            );

        //- Check if the intermediate excited molecule is in the ground
        //  vibrational state after the trial L-B redistribution
        if (trialVibLevelIntermediate == 0)
        {
            //- The intermediate molecule is 'formed' and is tested for
            //  ionisation
            Ec -= heatOfReactionAssociativeIonisationJoules_[0];

            const label ELevelIntermediate =
                cloud_.postCollisionElectronicEnergyLevel
                (
                    Ec,
                    cloud_.constProps(intermediateMoleculeId_).nElectronicLevels(),
                    omegaIntermediate,
                    EElistIntermediate,
                    gListIntermediate
                );

            //- Relative translational energy after electronic energy
            //  redistribution
            Ec -= EElistIntermediate[ELevelIntermediate];

            iMax = Ec/(physicoChemical::k.value()*thetaVIntermediate);

            if(iMax > 0)
            {
                label postCollisionVibLevel =
                    cloud_.postCollisionVibrationalEnergyLevel
                    (
                        true,
                        0.0,
                        iMax,
                        thetaVIntermediate,
                        cloud_.constProps(intermediateMoleculeId_).thetaD(),
                        cloud_.constProps(intermediateMoleculeId_).TrefZv()[0],
                        omegaIntermediate,
                        cloud_.constProps(intermediateMoleculeId_).Zref()[0],
                        Ec // TODO VINCENT add last optional params to all postCollVib
                    );

                //- Relative translational energy after vibrational energy
                //  redistribution
                Ec -= postCollisionVibLevel*thetaVIntermediate*physicoChemical::k.value();
            }

            const scalar energyRatio =
                cloud_.postCollisionRotationalEnergy
                (
                    cloud_.constProps(intermediateMoleculeId_).rotationalDegreesOfFreedom(),
                    2.5 - omegaIntermediate
                );

            scalar ERot = energyRatio*Ec;

            //- Relative translational energy after rotational energy
            //  redistribution
            Ec -= ERot;

            //- Redistribution finished, test it for ionisation
            if
            (
                Ec + EElistIntermediate[ELevelIntermediate]
              > heatOfReactionAssociativeIonisationJoules_[1]
            )
            {
                //- Add reaction to the list of competing reactions with
                //  probability reactionProbability
                reactionProbability = 1.0;
                totalReactionProbability += reactionProbability;
            }
        }
    }
}


void associativeIonisationQK::testReverseAssociativeIonisation
(
    const dsmcParcel& p,
    const dsmcParcel& q,
    const scalar translationalEnergy,
    scalar& collisionEnergy,
    scalar& totalReactionProbability,
    scalar& reactionProbability
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    scalar internalEnergyIonisedMolecule = 0.0;
    if (posIonisedMol_ == 0)
    {
        const scalar ERot = p.ERot();
        const scalar EVib = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
        const scalar EElec =
            cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];

        internalEnergyIonisedMolecule = ERot + EVib + EElec;
    }
    else
    {
        const scalar ERot = q.ERot();
        const scalar EVib = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());
        const scalar EElec =
            cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        internalEnergyIonisedMolecule = ERot + EVib + EElec;
    }

    const scalar omegaIntermediate = cloud_.constProps(intermediateMoleculeId_).omega();
    const scalar thetaVIntermediate = cloud_.constProps(intermediateMoleculeId_).thetaV_m(0);
    const scalarList& EElistIntermediate = cloud_.constProps(intermediateMoleculeId_).electronicEnergyList();
    const labelList& gListIntermediate = cloud_.constProps(intermediateMoleculeId_).electronicDegeneracyList();

    //- Collision energy is the sum of the relative kinetic energy and all
    //  internal energy modes of the ionised molecule
    scalar Ec = translationalEnergy + internalEnergyIonisedMolecule;

    //- Trial L-B redistribution
    const label trialELevelIntermediate =
        cloud_.postCollisionElectronicEnergyLevel
        (
            Ec,
            cloud_.constProps(intermediateMoleculeId_).nElectronicLevels(),
            omegaIntermediate,
            EElistIntermediate,
            gListIntermediate
        );

    //- Check if the intermediate excited molecule is in the ground
    //  electronic state after the trial L-B redistribution
    if (trialELevelIntermediate == 0)
    {
        //- The intermediate molecule is 'formed' and is tested for
        //  ionisation
        Ec -= heatOfReactionAssociativeIonisationJoules_[0];

        const label ELevelIntermediate =
            cloud_.postCollisionElectronicEnergyLevel
            (
                Ec,
                cloud_.constProps(intermediateMoleculeId_).nElectronicLevels(),
                omegaIntermediate,
                EElistIntermediate,
                gListIntermediate
            );

        //- Relative translational energy after electronic energy
        //  redistribution
        Ec -= EElistIntermediate[ELevelIntermediate];

        const scalar iMax = Ec/(physicoChemical::k.value()*thetaVIntermediate);

        label postCollisionVibLevel = 0;

        if(iMax > 0)
        {
            postCollisionVibLevel =
                cloud_.postCollisionVibrationalEnergyLevel
                (
                    true,
                    0.0,
                    iMax,
                    thetaVIntermediate,
                    cloud_.constProps(intermediateMoleculeId_).thetaD(),
                    cloud_.constProps(intermediateMoleculeId_).TrefZv()[0],
                    omegaIntermediate,
                    cloud_.constProps(intermediateMoleculeId_).Zref()[0],
                    Ec // TODO VINCENT add last optional params to all postCollVib
                );
        }

        //- Relative translational energy after vibrational energy
        //  redistribution
        const scalar EVibIntermediate = postCollisionVibLevel*thetaVIntermediate
           *physicoChemical::k.value();
        Ec -= EVibIntermediate;

        const scalar energyRatio =
            cloud_.postCollisionRotationalEnergy
            (
                cloud_.constProps(intermediateMoleculeId_).rotationalDegreesOfFreedom(),
                2.5 - omegaIntermediate
            );

        scalar ERot = energyRatio*Ec;

        //- Relative translational energy after rotational energy
        //  redistribution
        Ec -= ERot;

        //- Redistribution finished, test for dissociation
        const label imaxIntermediateMolecule =
            (Ec + EVibIntermediate)/(physicoChemical::k.value()*thetaVIntermediate);

        const label idIntermediateMolecule =
            cloud_.constProps(intermediateMoleculeId_).charDissQuantumLevel_m(0);

        if (imaxIntermediateMolecule > idIntermediateMolecule)
        {
            //- Add reaction to the list of competing reactions with
            //  probability reactionProbability
            reactionProbability = 1.0;
            totalReactionProbability += reactionProbability;
        }
    }
}


void associativeIonisationQK::forwardAssociativeIonisation
(
    dsmcParcel& p,
    dsmcParcel& q,
    scalar collisionEnergy
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    //- Forward associative ionisation
    nTotAssociativeIonisationReactions_++;
    nAssociativeIonisationReactionsPerTimeStep_++;

    if (allowSplitting_)
    {
        relax_ = false;

        vector UP = p.U();
        vector UQ = q.U();

        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        const scalar cRsqr = magSqr(UP - UQ);

        const scalar translationalEnergy = 0.5*mR*cRsqr;

        const scalar EEleP =
            cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        const scalar EEleQ =
            cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        const scalar Ecoll = translationalEnergy + EEleP + EEleQ
            - heatOfReactionAssociativeIonisationJoules_[0]
            - heatOfReactionAssociativeIonisationJoules_[1];

        //- Centre of mass velocity of molecules (pre-split)
        const vector& Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

        const scalar mProd1 = cloud_.constProps(productIdsAssociativeIonisation_[0]).mass();
        const scalar mProd2 = cloud_.constProps(productIdsAssociativeIonisation_[1]).mass();
        const scalar mRProducts = mProd1*mProd2/(mProd1 + mProd2);

        //- Assumption: no energy redistribution for both particles
        //  All the energy is stored in the translational mode
        const scalar relVel = sqrt(2.0*Ecoll/mRProducts);

        //- Variable Hard Sphere collision part for collision of molecules
        const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
        const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

        const vector& postCollisionRelU =
            relVel
           *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

        UP = Ucm + (postCollisionRelU*mProd2/(mProd1 + mProd2));
        UQ = Ucm - (postCollisionRelU*mProd1/(mProd1 + mProd2));

        //- p is originally a neutral atom and becomes an ionised molecule
        p.typeId() = productIdsAssociativeIonisation_[0];
        p.U() = UP;
        p.ERot() = 0.0;
        p.vibLevel().setSize
        (
            cloud_.constProps
            (
                intermediateMoleculeId_
            ).nVibrationalModes(),
            0
        );
        p.ELevel() = 0;

        //- q is originally a neutral atom and becomes an electron
        q.typeId() = productIdsAssociativeIonisation_[1];
        q.U() = UQ;
        q.ERot() = 0.0;
        q.vibLevel().setSize(0);
        q.ELevel() = 0;
    }
}


void associativeIonisationQK::reverseAssociativeIonisation
(
    dsmcParcel& p,
    dsmcParcel& q,
    scalar collisionEnergy
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    //- Reverse associative ionisation
    nTotAssociativeIonisationReactions_++;
    nAssociativeIonisationReactionsPerTimeStep_++;

    if (allowSplitting_)
    {
        relax_ = false;

        vector UP = p.U();
        vector UQ = q.U();

        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        const scalar cRsqr = magSqr(UP - UQ);

        const scalar translationalEnergy = 0.5*mR*cRsqr;

        scalar internalEnergyIonisedMolecule = 0.0;
        if (posIonisedMol_ == 0)
        {
            const scalar ERot = p.ERot();
            const scalar EVib = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
            const scalar EElec =
                cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];

            internalEnergyIonisedMolecule = ERot + EVib + EElec;
        }
        else
        {
            const scalar ERot = q.ERot();
            const scalar EVib = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());
            const scalar EElec =
                cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

            internalEnergyIonisedMolecule = ERot + EVib + EElec;
        }

        const scalar Ecoll = translationalEnergy + internalEnergyIonisedMolecule
            - heatOfReactionAssociativeIonisationJoules_[0]
            - heatOfReactionAssociativeIonisationJoules_[1];

        //- Centre of mass velocity of molecules (pre-split)
        const vector& Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

        const scalar mProd1 = cloud_.constProps(productIdsAssociativeIonisation_[0]).mass();
        const scalar mProd2 = cloud_.constProps(productIdsAssociativeIonisation_[1]).mass();
        const scalar mRProducts = mProd1*mProd2/(mProd1 + mProd2);

        //- Assumption: no energy redistribution for both particles
        //  All the energy is stored in the translational mode
        const scalar relVel = sqrt(2.0*Ecoll/mRProducts);

        //- Variable Hard Sphere collision part for collision of molecules
        const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
        const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

        const vector& postCollisionRelU =
            relVel
           *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );

        UP = Ucm + (postCollisionRelU*mProd2/(mProd1 + mProd2));
        UQ = Ucm - (postCollisionRelU*mProd1/(mProd1 + mProd2));

        //- p becomes a neutral atom
        p.typeId() = productIdsAssociativeIonisation_[0];
        p.U() = UP;
        p.ERot() = 0.0;
        p.vibLevel().setSize(0);
        p.ELevel() = 0;

        //- q becomes a neutral atom
        q.typeId() = productIdsAssociativeIonisation_[1];
        q.U() = UQ;
        q.ERot() = 0.0;
        q.vibLevel().setSize(0);
        q.ELevel() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
associativeIonisationQK::associativeIonisationQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    intermediateMoleculeId_(-1),
    productIdsAssociativeIonisation_(),
    forwardAssociativeIonisation_(false),
    posIonisedMol_(0),
    associativeIonisationStr_(word::null),
    nTotAssociativeIonisationReactions_(0),
    nAssociativeIonisationReactionsPerTimeStep_(0),
    heatOfReactionAssociativeIonisationJoules_(),
    volume_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

associativeIonisationQK::~associativeIonisationQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void associativeIonisationQK::initialConfiguration()
{
    setProperties();

    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA = cloud_.typeIdList()[productIdsAssociativeIonisation_[0]];
    const word& productB = cloud_.typeIdList()[productIdsAssociativeIonisation_[1]];

    associativeIonisationStr_ = "Associative ionisation reaction " + reactantA
        + " + " + reactantB + " --> " + productA + " + " + productB;
}


bool associativeIonisationQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
    //- Function used when setting the pair addressing matrix
    const label reactantPId = findIndex(reactantIds_, typeIdP);
    const label reactantQId = findIndex(reactantIds_, typeIdQ);

    //- If both indices were found in the list of reactants, there Ids will be
    //  different from -1
    if ((reactantPId != -1) && (reactantQId != -1))
    {
        //- Case of similar species colliding
        if((reactantPId == reactantQId) and (reactantIds_[0] == reactantIds_[1]))
        {
            return true;
        }

        //- Case of dissimilar species colliding
        if((reactantPId != reactantQId) and (reactantIds_[0] != reactantIds_[1]))
        {
            return true;
        }
    }

    return false;
}


void associativeIonisationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}


void associativeIonisationQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;

    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    //  If P is the first reactant A
    //  NB: Q is necessarily M otherwise this class would not have been selected
    if (typeIdP == reactantIds_[0])
    {
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);

        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;

        //- Possible reactions:
        // 1. Associative ionisation

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(1, 0.0);
        scalarList collisionEnergies(1, 0.0);

        if (forwardAssociativeIonisation_)
        {
            testForwardAssociativeIonisation
            (
                p,
                q,
                translationalEnergy,
                collisionEnergies[0],
                totalReactionProbability,
                reactionProbabilities[0]
            );
        }
        else
        {
            testReverseAssociativeIonisation
            (
                p,
                q,
                translationalEnergy,
                collisionEnergies[0],
                totalReactionProbability,
                reactionProbabilities[0]
            );
        }

        //- Decide if a reaction is to occur
        if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            //- Current reaction is to occur
            if (forwardAssociativeIonisation_)
            {
                forwardAssociativeIonisation(p, q, collisionEnergies[0]);
            }
            else
            {
                reverseAssociativeIonisation(p, q, collisionEnergies[0]);
            }
        }
    }
    else
    {
        //- If P is the second reactant M, then switch arguments in this
        //  function and P will be first
        associativeIonisationQK::reaction(q, p);
    }
}

void associativeIonisationQK::outputResults(const label& counterIndex)
{
    if (writeRatesToTerminal_)
    {
        //- measure density
        const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        labelList molsReactants(label(2), 0);

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                const label pos = findIndex(reactantIds_, p->typeId());

                if (pos != -1)
                {
                    molsReactants[pos]++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
        }

        scalarList numberDensities(2, cloud_.nParticle()/volume);
        numberDensities[0] *= molsReactants[0];
        numberDensities[1] *= molsReactants[1];

        label nTotAssociativeIonisationReactions =
            nTotAssociativeIonisationReactions_;
        label nAssociativeIonisationReactionsPerTimeStep =
            nAssociativeIonisationReactionsPerTimeStep_;

        const scalar deltaT = mesh_.time().deltaT().value();
        scalar factor = 0.0;

        if (reactantIds_[0] == reactantIds_[1] && numberDensities[0] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[0]
                   *volume
                );
        }
        else if (numberDensities[0] > 0.0 && numberDensities[1] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[1]
                   *volume
                );
        }

        if (Pstream::parRun())
        {
            //- Parallel communication
            reduce(molsReactants[0], sumOp<label>());
            reduce(molsReactants[1], sumOp<label>());
            reduce
            (
                nTotAssociativeIonisationReactions,
                sumOp<label>()
            );
            reduce
            (
                nAssociativeIonisationReactionsPerTimeStep,
                sumOp<label>()
            );
        }

        const scalar reactionRateAssociativeIonisation =
            factor*nTotAssociativeIonisationReactions;

        Info<< associativeIonisationStr_
            << ", reaction rate = " << reactionRateAssociativeIonisation
            << ", nReactions = "
            << nAssociativeIonisationReactionsPerTimeStep
            << endl;
    }
    else
    {
        label nTotAssociativeIonisationReactions =
            nTotAssociativeIonisationReactions_;
        label nAssociativeIonisationReactionsPerTimeStep =
            nAssociativeIonisationReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            //- Parallel communication
            reduce(nTotAssociativeIonisationReactions, sumOp<label>());
            reduce(nAssociativeIonisationReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotAssociativeIonisationReactions > 0)
        {
            Info<< associativeIonisationStr_
                << " is active, nReactions this time step = "
                << nAssociativeIonisationReactionsPerTimeStep
                << endl;
         }
    }

    nAssociativeIonisationReactionsPerTimeStep_ = 0;
}

}
// End namespace Foam

// ************************************************************************* //
