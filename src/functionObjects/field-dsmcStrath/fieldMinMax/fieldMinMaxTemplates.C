/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "fieldMinMax.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldMinMax::output
(
    const word& fieldName,
    const word& outputName,
    const vector& minC,
    const vector& maxC,
    const label minProci,
    const label maxProci,
    const Type& minValue,
    const Type& maxValue
)
{
    OFstream& file = this->file();

    if (location_)
    {
        writeTime(file());

        writeTabbed(file, fieldName);

        file<< token::TAB << minValue
            << token::TAB << minC;

        if (Pstream::parRun())
        {
            file<< token::TAB << minProci;
        }

        file<< token::TAB << maxValue
            << token::TAB << maxC;

        if (Pstream::parRun())
        {
            file<< token::TAB << maxProci;
        }

        file<< endl;

        Log << "    min(" << outputName << ") = " << minValue
            << " at location " << minC;

        if (Pstream::parRun())
        {
            Log << " on processor " << minProci;
        }

        Log << nl << "    max(" << outputName << ") = " << maxValue
            << " at location " << maxC;

        if (Pstream::parRun())
        {
            Log << " on processor " << maxProci;
        }
    }
    else
    {
        file<< token::TAB << minValue << token::TAB << maxValue;

        Log << "    min/max(" << outputName << ") = "
            << minValue << ' ' << maxValue;
    }

    Log << endl;

    // Write state/results information
    word nameStr('(' + outputName + ')');
    this->setResult("min" + nameStr, minValue);
    this->setResult("min" + nameStr + "_position", minC);
    this->setResult("min" + nameStr + "_processor", minProci);
    this->setResult("max" + nameStr, maxValue);
    this->setResult("max" + nameStr + "_position", maxC);
    this->setResult("max" + nameStr + "_processor", maxProci);
}


template<class Type>
void Foam::functionObjects::fieldMinMax::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const label proci = Pstream::myProcNo();

        const fieldType& field = lookupObject<fieldType>(fieldName);

        const volVectorField::Boundary& CfBoundary =
            mesh_.C().boundaryField();

        switch (mode)
        {
            case mdMag:
            {
                const volScalarField magField(mag(field));
                const volScalarField::Boundary& magFieldBoundary =
                    magField.boundaryField();

                scalarList minVs(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                label minProci = findMin(magField);
                minVs[proci] = magField[minProci];
                minCs[proci] = mesh_.C()[minProci];


                labelList maxIs(Pstream::nProcs());
                scalarList maxVs(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                label maxProci = findMax(magField);
                maxVs[proci] = magField[maxProci];
                maxCs[proci] = mesh_.C()[maxProci];

                forAll(magFieldBoundary, patchi)
                {
                    const scalarField& mfp = magFieldBoundary[patchi];
                    if (mfp.size())
                    {
                        const vectorField& Cfp = CfBoundary[patchi];

                        label minPi = findMin(mfp);
                        if (mfp[minPi] < minVs[proci])
                        {
                            minVs[proci] = mfp[minPi];
                            minCs[proci] = Cfp[minPi];
                        }

                        label maxPi = findMax(mfp);
                        if (mfp[maxPi] > maxVs[proci])
                        {
                            maxVs[proci] = mfp[maxPi];
                            maxCs[proci] = Cfp[maxPi];
                        }
                    }
                }

                Pstream::gatherList(minVs);
                Pstream::scatterList(minVs);
                Pstream::gatherList(minCs);
                Pstream::scatterList(minCs);

                Pstream::gatherList(maxVs);
                Pstream::scatterList(maxVs);
                Pstream::gatherList(maxCs);
                Pstream::scatterList(maxCs);

                label mini = findMin(minVs);
                scalar minValue = minVs[mini];
                const vector& minC = minCs[mini];

                label maxi = findMax(maxVs);
                scalar maxValue = maxVs[maxi];
                const vector& maxC = maxCs[maxi];

                output
                (
                    fieldName,
                    word("mag(" + fieldName + ")"),
                    minC,
                    maxC,
                    mini,
                    maxi,
                    minValue,
                    maxValue
                );
                break;
            }
            case mdCmpt:
            {
                const typename fieldType::Boundary&
                    fieldBoundary = field.boundaryField();

                List<Type> minVs(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                label minProci = findMin(field);
                minVs[proci] = field[minProci];
                minCs[proci] = mesh_.C()[minProci];

                Pstream::gatherList(minVs);
                Pstream::gatherList(minCs);

                List<Type> maxVs(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                label maxProci = findMax(field);
                maxVs[proci] = field[maxProci];
                maxCs[proci] = mesh_.C()[maxProci];

                forAll(fieldBoundary, patchi)
                {
                    const Field<Type>& fp = fieldBoundary[patchi];
                    if (fp.size())
                    {
                        const vectorField& Cfp = CfBoundary[patchi];

                        label minPi = findMin(fp);
                        if (fp[minPi] < minVs[proci])
                        {
                            minVs[proci] = fp[minPi];
                            minCs[proci] = Cfp[minPi];
                        }

                        label maxPi = findMax(fp);
                        if (fp[maxPi] > maxVs[proci])
                        {
                            maxVs[proci] = fp[maxPi];
                            maxCs[proci] = Cfp[maxPi];
                        }
                    }
                }

                Pstream::gatherList(minVs);
                Pstream::scatterList(minVs);
                Pstream::gatherList(minCs);
                Pstream::scatterList(minCs);

                Pstream::gatherList(maxVs);
                Pstream::scatterList(maxVs);
                Pstream::gatherList(maxCs);
                Pstream::scatterList(maxCs);

                label mini = findMin(minVs);
                Type minValue = minVs[mini];
                const vector& minC = minCs[mini];

                label maxi = findMax(maxVs);
                Type maxValue = maxVs[maxi];
                const vector& maxC = maxCs[maxi];

                output
                (
                    fieldName,
                    fieldName,
                    minC,
                    maxC,
                    mini,
                    maxi,
                    minValue,
                    maxValue
                );
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown min/max mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
