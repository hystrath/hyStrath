/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

Diffuse wall where ions are neutralised and current into wall used to calculate
electric potential boundary condition.

\*---------------------------------------------------------------------------*/

#include "pdDiffuseNeutralisingWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdDiffuseNeutralisingWallPatch, 0);

addToRunTimeSelectionTable(pdPatchBoundary, pdDiffuseNeutralisingWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdDiffuseNeutralisingWallPatch::pdDiffuseNeutralisingWallPatch
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    //nTotReactions_(0),
    //nReactionsPerTimeStep_(0),
    typeElec_(),
    ionIds_(),
    productsToNeut_()
    //writeRatesToTerminal_(false)

{
    writeInTimeDir_             = false;
    writeInCase_                = false;
    measurePropertiesAtWall_    = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdDiffuseNeutralisingWallPatch::~pdDiffuseNeutralisingWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void pdDiffuseNeutralisingWallPatch::initialConfiguration()
{

}

void pdDiffuseNeutralisingWallPatch::calculateProperties()
{
    if(typeElec_ == 1) //insulated
    {
        Info << "Insulated material not implimented!" << endl;
    }
    else if(typeElec_ == 2) // conductor
    {
        //const scalar e  = electromagnetic::e.value();   //elementary charge
        scalar totalQ   = 0.0;

        pdEmFields& eM = cloud_.emFields();

        forAll(eM.wallQ_.boundaryField()[patchId_],fI)
        {
            totalQ += eM.wallQ_.boundaryField()[patchId_][fI];
            //totalQ -= e*eM.efRhoBF()[patchId_][fI] * mag(patch.faceAreas()[fI]); //- need to check whether this will work with efmodel
        }

        if(Pstream::parRun())
        {
            reduce(totalQ, sumOp<scalar>());  //total charge on all processors
        }

        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        forAll(eM.wallQ_.boundaryField()[patchId_],fI)
        {
            eM.wallQ_.boundaryFieldRef()[patchId_][fI] = totalQ*mag(patch.faceAreas()[fI])/totalPatchSurfaceArea_;
        }
    }
    else
    {
        pdEmFields& eM = cloud_.emFields();
        eM.wallQ_.boundaryFieldRef()[patchId_] = 0.0;
    }
}

void pdDiffuseNeutralisingWallPatch::controlParticle(pdParcel& p, pdParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);

    label iD = findIndex(ionIds_, p.typeId());

    if(iD == -1)
    {
        diffuseReflect(p, td);
    }
    else
    {
        //neutralise the particle
        p.typeId() = productsToNeut_[iD][0];
        p.A() = vector::zero; //zero the acceleration vector

        diffuseReflect(p, td);

        /*
        forAll(productsToNeut_[iD],n)
        {
            //currently only supports single species products
        }
        */
    }

    measurePropertiesAfterControl(p);
}

void pdDiffuseNeutralisingWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void pdDiffuseNeutralisingWallPatch::diffuseReflect(pdParcel& p, pdParcel::trackingData& td)
{
    vector& U = p.U();

    scalar& ERot = p.ERot();

    scalar& EVib = p.EVib();

    label typeId = p.typeId();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.sample01<scalar>()),
            U.y()*(0.8 + 0.2*rndGen.sample01<scalar>()),
            U.z()*(0.8 + 0.2*rndGen.sample01<scalar>())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    const scalar& T = temperature_;

    scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();

    scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();

    U =
        sqrt(physicoChemical::k.value()*T/mass)
    *(
            rndGen.GaussNormal<scalar>()*tw1
        + rndGen.GaussNormal<scalar>()*tw2
        - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))*nw
        );

    U += velocity_;

    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);

    EVib = cloud_.equipartitionVibrationalEnergy(T, vibrationalDof, typeId);
}


void pdDiffuseNeutralisingWallPatch::updateProperties(const dictionary& newDict)
{
//     Info << "old dictionary: " << propsDict_  << endl;

    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void pdDiffuseNeutralisingWallPatch::setProperties()
{
    /** Electromagnetic properties **/
    typeElec_ = readScalar(propsDict_.lookup("typeElec"));

    /*if(typeElec_ == 1)
    {
        Info << "   Warning: " << patchName() << " has been set as insulated (typeElec = 2). Remember to select the appropriate fv boundary scheme." << endl;
    }else if(typeElec_ == 2)
    {
        Info << "   Warning: " << patchName() << " has been set as conducting (typeElec = 2). Remember to select the appropriate fv boundary scheme." << endl;
    }
    else
    {
        Info << "   Warning: " << patchName() << " has been set as grounded. Remember to select the appropriate fv boundary scheme or if this is incorrect, change typeElec (1. insulted, 2. conducting)."  << endl;
    }*/


    /** Diffuse wall properties **/
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));

    /** Reacton properties **/
    //lookup list of ions to neutralise and build list of their ids
    const List<word> ionSpecies (propsDict_.lookup("ionsToNeutralise"));

    forAll (ionSpecies, i)
    {
        forAll (ionSpecies, p)
        {
            if (i != p && ionSpecies[i] == ionSpecies[p])
            {
                FatalErrorIn("pdDiffuseNeutralisingWallPatch::setProperties()")
                    << "ions cannot be same species." << nl
                    << exit(FatalError);
            }
        }
    }

    ionIds_.setSize(ionSpecies.size(),-1);

    //allowAbsorption_        = Switch(propsDict_.lookup("allowAbsorption"));
    //writeRatesToTerminal_    = Switch(propsDict_.lookup("writeRatesToTerminal"))

    forAll (ionIds_, i)
    {
        ionIds_[i] = findIndex(cloud_.typeIdList(),ionSpecies[i]);

        // check that ions belong to the typeIdList (constant/pdProperties)
        if(ionIds_[i] == -1)
        {
            FatalErrorIn("pdDiffuseNeutralisingWallPatch::setProperties()")
                << "Cannot find type id: " << ionSpecies[i] << nl
                << exit(FatalError);
        }

    }

     //lookup list of neutralisation products and build list of their ids
    List< List<word> > productNeutrals (propsDict_.lookup("productsOfNeutralisation"));

    if(productNeutrals.size() != ionIds_.size())
    {
        FatalErrorIn("pdDiffuseNeutralisingWallPatch::setProperties()")
            << "number of ion to be neutralised = " << ionIds_.size()
            << " is not the same as the number of products = " << productNeutrals.size()
            << exit(FatalError);
    }

    productsToNeut_.setSize(productNeutrals.size());

    forAll(productNeutrals, n)
    {
        const List<word>& productsForNeut = productNeutrals[n];

        productsToNeut_[n].setSize(productsForNeut.size(), -1);

        forAll(productsToNeut_[n], p)
        {
            productsToNeut_[n][p] = findIndex(cloud_.typeIdList(), productsForNeut[p]);

            if(productsToNeut_[n][p] == -1)
            {
                FatalErrorIn("pdDiffuseNeutralisingWallPatch::setProperties()")
                    << "Cannot find type id: " << productsForNeut[p] << nl
                    << exit(FatalError);
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
