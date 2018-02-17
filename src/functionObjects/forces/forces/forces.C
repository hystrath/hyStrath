/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "forces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

#include "multi2Thermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forces, 0);

    addToRunTimeSelectionTable(functionObject, forces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObjects::forces::fieldName(const word& name) const
{
    return this->name() + ":" + name;
}


void Foam::functionObjects::forces::createFiles()
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !forceFilePtr_.valid())
    {
        forceFilePtr_ = createFile("force");
        writeIntegratedHeader("Force", forceFilePtr_());
        momentFilePtr_ = createFile("moment");
        writeIntegratedHeader("Moment", momentFilePtr_());

        if (nBin_ > 1)
        {
            forceBinFilePtr_ = createFile("forceBin");
            writeBinHeader("Force", forceBinFilePtr_());
            momentBinFilePtr_ = createFile("momentBin");
            writeBinHeader("Moment", momentBinFilePtr_());
        }

        if (localSystem_)
        {
            localForceFilePtr_ = createFile("localForce");
            writeIntegratedHeader("Force", localForceFilePtr_());
            localMomentFilePtr_ = createFile("localMoment");
            writeIntegratedHeader("Moment", localMomentFilePtr_());

            if (nBin_ > 1)
            {
                localForceBinFilePtr_ = createFile("localForceBin");
                writeBinHeader("Force", localForceBinFilePtr_());
                localMomentBinFilePtr_ = createFile("localMomentBin");
                writeBinHeader("Moment", localMomentBinFilePtr_());
            }
        }
    }
}


void Foam::functionObjects::forces::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(total_x total_y total_z)");
    writeTabbed(os, "(pressure_x pressure_y pressure_z)");
    writeTabbed(os, "(viscous_x viscous_y viscous_z)");

    if (porosity_)
    {
        writeTabbed(os, "(porous_x porous_y porous_z)");
    }

    os  << endl;
}


void Foam::functionObjects::forces::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header + " bins");
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointi)
    {
        binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
        os  << tab << binPoints[pointi].x();
    }
    os  << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].y();
    }
    os  << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].z();
    }
    os  << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    for (label j = 0; j < nBin_; j++)
    {
        const word jn(Foam::name(j) + ':');
        os  << tab << jn << "(total_x total_y total_z)"
            << tab << jn << "(pressure_x pressure_y pressure_z)"
            << tab << jn << "(viscous_x viscous_y viscous_z)";

        if (porosity_)
        {
            os  << tab << jn << "(porous_x porous_y porous_z)";
        }
    }

    os << endl;
}



void Foam::functionObjects::forces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_ << " or p:" << pName_
                << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_
                << exit(FatalError);
        }
    }

    initialiseBins();

    initialised_ = true;
}


void Foam::functionObjects::forces::initialiseBins()
{
    if (nBin_ > 1)
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Determine extents of patches
        binMin_ = GREAT;
        scalar binMax = -GREAT;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            const polyPatch& pp = pbm[patchi];
            scalarField d(pp.faceCentres() & binDir_);
            binMin_ = min(min(d), binMin_);
            binMax = max(max(d), binMax);
        }

        // Include porosity
        if (porosity_)
        {
            const HashTable<const porosityModel*> models =
                obr_.lookupClass<porosityModel>();

            const scalarField dd(mesh_.C() & binDir_);

            forAllConstIter(HashTable<const porosityModel*>, models, iter)
            {
                const porosityModel& pm = *iter();
                const labelList& cellZoneIDs = pm.cellZoneIDs();

                forAll(cellZoneIDs, i)
                {
                    label zonei = cellZoneIDs[i];
                    const cellZone& cZone = mesh_.cellZones()[zonei];
                    const scalarField d(dd, cZone);
                    binMin_ = min(min(d), binMin_);
                    binMax = max(max(d), binMax);
                }
            }
        }

        reduce(binMin_, minOp<scalar>());
        reduce(binMax, maxOp<scalar>());

        // Slightly boost binMax so that region of interest is fully
        // within bounds
        binMax = 1.0001*(binMax - binMin_) + binMin_;

        binDx_ = (binMax - binMin_)/scalar(nBin_);

        // Create the bin points used for writing
        binPoints_.setSize(nBin_);
        forAll(binPoints_, i)
        {
            binPoints_[i] = (i + 0.5)*binDir_*binDx_;
        }
    }

    // Allocate storage for forces and moments
    forAll(force_, i)
    {
        force_[i].setSize(nBin_, vector::zero);
        moment_[i].setSize(nBin_, vector::zero);
    }
}


void Foam::functionObjects::forces::resetFields()
{
    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    moment_[0] = Zero;
    moment_[1] = Zero;
    moment_[2] = Zero;

    if (writeFields_)
    {
        volVectorField& force =
            const_cast<volVectorField&>
            (
                lookupObject<volVectorField>(fieldName("force"))
            );

        force == dimensionedVector("0", force.dimensions(), Zero);

        volVectorField& moment =
            const_cast<volVectorField&>
            (
                lookupObject<volVectorField>(fieldName("moment"))
            );

        moment == dimensionedVector("0", moment.dimensions(), Zero);
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (foundObject<multi2Thermo>(multi2Thermo::dictName))
    {
        const multi2Thermo& thermo =
            lookupObject<multi2Thermo>(multi2Thermo::dictName);

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::mu() const
{
    if (foundObject<multi2Thermo>(basic2Thermo::dictName))
    {
        const multi2Thermo& thermo =
             lookupObject<multi2Thermo>(basic2Thermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const transportModel& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::rho() const
{
    /*if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(lookupObject<volScalarField>(rhoName_));
    }*/
    
    return(lookupObject<volScalarField>("rho")); // NEW VINCENT: compressible cases with forceCoeff active
}


Foam::scalar Foam::functionObjects::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


void Foam::functionObjects::forces::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        force_[0][0] += sum(fN);
        force_[1][0] += sum(fT);
        force_[2][0] += sum(fP);
        moment_[0][0] += sum(Md^fN);
        moment_[1][0] += sum(Md^fT);
        moment_[2][0] += sum(Md^fP);
    }
    else
    {
        scalarField dd((d & binDir_) - binMin_);

        forAll(dd, i)
        {
            label bini = min(max(floor(dd[i]/binDx_), 0), force_[0].size() - 1);

            force_[0][bini] += fN[i];
            force_[1][bini] += fT[i];
            force_[2][bini] += fP[i];
            moment_[0][bini] += Md[i]^fN[i];
            moment_[1][bini] += Md[i]^fT[i];
            moment_[2][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::functionObjects::forces::addToFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        const_cast<volVectorField&>
        (
            lookupObject<volVectorField>(fieldName("force"))
        );

    vectorField& pf = force.boundaryFieldRef()[patchi];
    pf += fN + fT + fP;

    volVectorField& moment =
        const_cast<volVectorField&>
        (
            lookupObject<volVectorField>(fieldName("moment"))
        );

    vectorField& pm = moment.boundaryFieldRef()[patchi];
    pm += Md;
}


void Foam::functionObjects::forces::addToFields
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    volVectorField& force =
        const_cast<volVectorField&>
        (
            lookupObject<volVectorField>(fieldName("force"))
        );

    volVectorField& moment =
        const_cast<volVectorField&>
        (
            lookupObject<volVectorField>(fieldName("moment"))
        );

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i] + fP[i];
        moment[celli] += Md[i];
    }
}


void Foam::functionObjects::forces::writeIntegratedForceMoment
(
    const string& descriptor,
    const vectorField& fm0,
    const vectorField& fm1,
    const vectorField& fm2,
    autoPtr<OFstream>& osPtr
) const
{
    vector pressure = sum(fm0);
    vector viscous = sum(fm1);
    vector porous = sum(fm2);
    vector total = pressure + viscous + porous;

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << total << nl
        << "        Pressure : " << pressure << nl
        << "        Viscous  : " << viscous << nl;

    if (porosity_)
    {
        Log << "        Porous   : " << porous << nl;
    }

    if (writeToFile())
    {
        Ostream& os = osPtr();

        os  << obr_.time().value()
            << tab << total
            << tab << pressure
            << tab << viscous;

        if (porosity_)
        {
            os  << tab << porous;
        }

        os  << endl;
    }
}


void Foam::functionObjects::forces::writeForces()
{
    Log << type() << " " << name() << " write:" << nl;

    writeIntegratedForceMoment
    (
        "forces",
        force_[0],
        force_[1],
        force_[2],
        forceFilePtr_
    );

    writeIntegratedForceMoment
    (
        "moments",
        moment_[0],
        moment_[1],
        moment_[2],
        momentFilePtr_
    );

    if (localSystem_)
    {
        writeIntegratedForceMoment
        (
            "local forces",
            coordSys_.localVector(force_[0]),
            coordSys_.localVector(force_[1]),
            coordSys_.localVector(force_[2]),
            localForceFilePtr_
        );

        writeIntegratedForceMoment
        (
            "local moments",
            coordSys_.localVector(moment_[0]),
            coordSys_.localVector(moment_[1]),
            coordSys_.localVector(moment_[2]),
            localMomentFilePtr_
        );
    }

    Log << endl;
}


void Foam::functionObjects::forces::writeBinnedForceMoment
(
    const List<Field<vector>>& fm,
    autoPtr<OFstream>& osPtr
) const
{
    if ((nBin_ == 1) || !writeToFile())
    {
        return;
    }

    List<Field<vector>> f(fm);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];
        }
    }

    Ostream& os = osPtr();

    writeTime(os);

    forAll(f[0], i)
    {
        vector total = f[0][i] + f[1][i] + f[2][i];

        os  << tab << total
            << tab << f[0][i]
            << tab << f[1][i];

        if (porosity_)
        {
            os  << tab << f[2][i];
        }
    }

    os  << nl;
}


void Foam::functionObjects::forces::writeBins()
{
    writeBinnedForceMoment(force_, forceBinFilePtr_);
    writeBinnedForceMoment(moment_, momentBinFilePtr_);

    if (localSystem_)
    {
        List<Field<vector>> lf(3);
        List<Field<vector>> lm(3);
        lf[0] = coordSys_.localVector(force_[0]);
        lf[1] = coordSys_.localVector(force_[1]);
        lf[2] = coordSys_.localVector(force_[2]);
        lm[0] = coordSys_.localVector(moment_[0]);
        lm[1] = coordSys_.localVector(moment_[1]);
        lm[2] = coordSys_.localVector(moment_[2]);

        writeBinnedForceMoment(lf, localForceBinFilePtr_);
        writeBinnedForceMoment(lm, localMomentBinFilePtr_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forces::forces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    force_(3),
    moment_(3),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    localForceFilePtr_(),
    localMomentFilePtr_(),
    localForceBinFilePtr_(),
    localMomentBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


Foam::functionObjects::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    force_(3),
    moment_(3),
    forceFilePtr_(),
    momentFilePtr_(),
    forceBinFilePtr_(),
    momentBinFilePtr_(),
    localForceFilePtr_(),
    localMomentFilePtr_(),
    localForceBinFilePtr_(),
    localMomentBinFilePtr_(),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }

/*
    // Turn off writing to file
    writeToFile_ = false;

    forAll(force_, i)
    {
        force_[i].setSize(nBin_, vector::zero);
        moment_[i].setSize(nBin_, vector::zero);
    }
*/
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forces::~forces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fD", "fD");
    }
    else
    {
        // Optional entries U and p
        pName_ = dict.lookupOrDefault<word>("p", "p");
        UName_ = dict.lookupOrDefault<word>("U", "U");
        rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = readScalar(dict.lookup("rhoInf"));
        }

        // Reference pressure, 0 by default
        pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
    }

    coordSys_.clear();

    // Centre of rotation for moment calculations
    // specified directly, from coordinate system, or implicitly (0 0 0)
    if (!dict.readIfPresent<point>("CofR", coordSys_.origin()))
    {
        coordSys_ = coordinateSystem(obr_, dict);
        localSystem_ = true;
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        binDict.lookup("nBin") >> nBin_;

        if (nBin_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Number of bins (nBin) must be zero or greater"
                << exit(FatalIOError);
        }
        else if (nBin_ == 0)
        {
            // Case of no bins equates to a single bin to collect all data
            nBin_ = 1;
        }
        else
        {
            binDict.lookup("cumulative") >> binCumulative_;
            binDict.lookup("direction") >> binDir_;
            binDir_ /= mag(binDir_);
        }
    }

    writeFields_ = dict.lookupOrDefault("writeFields", false);

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        volVectorField* forcePtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("force"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("0", dimForce, Zero)
            )
        );

        mesh_.objectRegistry::store(forcePtr);

        volVectorField* momentPtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("moment"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("0", dimForce*dimLength, Zero)
            )
        );

        mesh_.objectRegistry::store(momentPtr);
    }

    return true;
}


void Foam::functionObjects::forces::calcForcesMoment()
{
    initialise();

    resetFields();

    if (directForceDensity_)
    {
        const volVectorField& fD = lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - coordSys_.origin()
            );

            scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Tangential force (total force minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            // Porous force
            vectorField fP(Md.size(), Zero);

            addToFields(patchi, Md, fN, fT, fP);

            applyBins(Md, fN, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }
    else
    {
        const volScalarField& p = lookupObject<volScalarField>(pName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::Boundary& devRhoReffb
            = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();

            vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - coordSys_.origin()
            );

            vectorField fN
            (
                rho(p)*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
            );

            vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            vectorField fP(Md.size(), Zero);

            addToFields(patchi, Md, fN, fT, fP);

            applyBins(Md, fN, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }

    if (porosity_)
    {
        const volVectorField& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            // Non-const access required if mesh is changing
            porosityModel& pm = const_cast<porosityModel&>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zonei = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - coordSys_.origin());

                const vectorField fDummy(Md.size(), Zero);

                addToFields(cZone, Md, fDummy, fDummy, fP);

                applyBins(Md, fDummy, fDummy, fP, d);
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::functionObjects::forces::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]);
}


Foam::vector Foam::functionObjects::forces::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]);
}


bool Foam::functionObjects::forces::execute()
{
    calcForcesMoment();

    if (Pstream::master())
    {
        createFiles();

        writeForces();

        writeBins();

        Log << endl;
    }

    // Write state/results information
    setResult("normalForce", sum(force_[0]));
    setResult("tangentialForce", sum(force_[1]));
    setResult("porousForce", sum(force_[2]));

    setResult("normalMoment", sum(moment_[0]));
    setResult("tangentialMoment", sum(moment_[1]));
    setResult("porousMoment", sum(moment_[2]));

    return true;
}


bool Foam::functionObjects::forces::write()
{
    if (writeFields_)
    {
        lookupObject<volVectorField>(fieldName("force")).write();
        lookupObject<volVectorField>(fieldName("moment")).write();
    }

    return true;
}


// ************************************************************************* //
