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

#include "fieldAverage.H"
#include "volFields.H"
#include "fieldAverageItem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldAverage, 0);
    addToRunTimeSelectionTable(functionObject, fieldAverage, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldAverage::resetFields()
{
    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            if (obr().found(faItems_[i].meanFieldName()))
            {
                obr().checkOut(*obr()[faItems_[i].meanFieldName()]);
            }
        }

        if (faItems_[i].prime2Mean())
        {
            if (obr().found(faItems_[i].prime2MeanFieldName()))
            {
                obr().checkOut(*obr()[faItems_[i].prime2MeanFieldName()]);
            }
        }
    }
}


void Foam::functionObjects::fieldAverage::initialize()
{
    resetFields();

    Log << type() << " " << name() << ":" << nl;

    // Add mean fields to the field lists
    forAll(faItems_, fieldi)
    {
        addMeanField<scalar>(fieldi);
        addMeanField<vector>(fieldi);
        addMeanField<sphericalTensor>(fieldi);
        addMeanField<symmTensor>(fieldi);
        addMeanField<tensor>(fieldi);
    }

    // Add prime-squared mean fields to the field lists
    forAll(faItems_, fieldi)
    {
        addPrime2MeanField<scalar, scalar>(fieldi);
        addPrime2MeanField<vector, symmTensor>(fieldi);
    }

    forAll(faItems_, fieldi)
    {
        if (!faItems_[fieldi].active())
        {
            WarningInFunction
                << "Field " << faItems_[fieldi].fieldName()
                << " not found in database for averaging";
        }
    }

    // ensure first averaging works unconditionally
    prevTimeIndex_ = -1;

    Log << endl;
    initialised_ = true;
}


void Foam::functionObjects::fieldAverage::restart()
{
    Log << "    Restarting averaging at time " << obr().time().timeName()
        << nl << endl;

    totalIter_.clear();
    totalIter_.setSize(faItems_.size(), 1);

    totalTime_.clear();
    totalTime_.setSize(faItems_.size(), obr().time().deltaTValue());

    initialize();
}


void Foam::functionObjects::fieldAverage::calcAverages()
{
    if (!initialised_)
    {
        initialize();
    }

    const label currentTimeIndex = obr().time().timeIndex();
    const scalar currentTime = obr().time().value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }
    else
    {
        prevTimeIndex_ = currentTimeIndex;
    }

    bool doRestart = false;
    if (periodicRestart_ && currentTime > restartPeriod_*periodIndex_)
    {
        doRestart = true;
        periodIndex_++;
    }

    if (currentTime >= restartTime_)
    {
        doRestart = true;      // Restart is overdue.
        restartTime_ = GREAT;  // Avoid triggering again
    }

    if (doRestart)
    {
        restart();
    }

    Log << type() << " " << name() << " write:" << nl
        << "    Calculating averages" << nl;

    addMeanSqrToPrime2Mean<scalar, scalar>();
    addMeanSqrToPrime2Mean<vector, symmTensor>();

    calculateMeanFields<scalar>();
    calculateMeanFields<vector>();
    calculateMeanFields<sphericalTensor>();
    calculateMeanFields<symmTensor>();
    calculateMeanFields<tensor>();

    calculatePrime2MeanFields<scalar, scalar>();
    calculatePrime2MeanFields<vector, symmTensor>();

    forAll(faItems_, fieldi)
    {
        totalIter_[fieldi]++;
        totalTime_[fieldi] += obr().time().deltaTValue();
    }

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAverages() const
{
    Log << "    Writing average fields" << endl;

    writeFields<scalar>();
    writeFields<vector>();
    writeFields<sphericalTensor>();
    writeFields<symmTensor>();
    writeFields<tensor>();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAveragingProperties()
{
    forAll(faItems_, fieldi)
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        dictionary propsDict;
        propsDict.add("totalIter", totalIter_[fieldi]);
        propsDict.add("totalTime", totalTime_[fieldi]);
        setProperty(fieldName, propsDict);
    }
}


void Foam::functionObjects::fieldAverage::readAveragingProperties()
{
    totalIter_.clear();
    totalIter_.setSize(faItems_.size(), 1);

    totalTime_.clear();
    totalTime_.setSize(faItems_.size(), obr().time().deltaTValue());

    if (restartOnRestart_ || restartOnOutput_)
    {
        Info<< "    Starting averaging at time " << obr().time().timeName()
            << nl;
    }
    else
    {
        Info<< "    Restarting averaging for fields:" << nl;


        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            if (foundProperty(fieldName))
            {
                dictionary fieldDict;
                getDict(fieldName, fieldDict);

                totalIter_[fieldi] = readLabel(fieldDict.lookup("totalIter"));
                totalTime_[fieldi] = readScalar(fieldDict.lookup("totalTime"));

                Info<< "        " << fieldName
                    << " iters = " << totalIter_[fieldi]
                    << " time = " << totalTime_[fieldi] << nl;
            }
            else
            {
                Info<< "        " << fieldName
                    << ": starting averaging at time "
                    << obr().time().timeName() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::fieldAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    prevTimeIndex_(-1),
    restartOnRestart_(false),
    restartOnOutput_(false),
    periodicRestart_(false),
    restartPeriod_(GREAT),
    restartTime_(GREAT),
    initialised_(false),
    faItems_(),
    totalIter_(),
    totalTime_(),
    periodIndex_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::~fieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // Make certain that the values are consistent with the defaults:
    initialised_ = false;
    restartOnRestart_ = false;
    restartOnOutput_  = false;
    periodicRestart_  = false;
    restartPeriod_    = GREAT;
    restartTime_      = GREAT;

    Info<< type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    dict.readIfPresent("restartOnOutput",  restartOnOutput_);
    dict.readIfPresent("periodicRestart",  periodicRestart_);
    dict.lookup("fields") >> faItems_;

    const scalar currentTime = obr().time().value();

    if (periodicRestart_)
    {
        dict.lookup("restartPeriod") >> restartPeriod_;

        if (restartPeriod_ > 0)
        {
            // Determine the appropriate interval for the next restart
            periodIndex_ = 1;
            while (currentTime > restartPeriod_*periodIndex_)
            {
                ++periodIndex_;
            }

            Info<< "    Restart period " << restartPeriod_
                << " - next restart at " << (restartPeriod_*periodIndex_)
                << nl << endl;
        }
        else
        {
            periodicRestart_ = false;

            Info<< "    Restart period " << restartPeriod_
                << " - ignored"
                << nl << endl;
        }
    }

    if (dict.readIfPresent("restartTime", restartTime_))
    {
        if (currentTime > restartTime_)
        {
            // The restart time is already in the past - ignore
            restartTime_ = GREAT;
        }
        else
        {
            Info<< "    Restart scheduled at time " << restartTime_
                << nl << endl;
        }
    }

    readAveragingProperties();

    Info<< endl;

    return true;
}


bool Foam::functionObjects::fieldAverage::execute()
{
    calcAverages();

    return true;
}


bool Foam::functionObjects::fieldAverage::write()
{
    writeAverages();
    writeAveragingProperties();

    if (restartOnOutput_)
    {
        restart();
    }

    return true;
}


// ************************************************************************* //
