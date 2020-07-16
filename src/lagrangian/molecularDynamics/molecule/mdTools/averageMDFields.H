/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
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

\*---------------------------------------------------------------------------*/

if (runTime.outputTime())
{
    /*-----------------------------------------------------------------------*\
        Number density
    \*-----------------------------------------------------------------------*/

    scalarField totalRhoN_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesRhoN, rN)
    {
        allSpeciesRhoN[rN].internalField() =
            allSpeciesN_RU[rN]
            /mesh.cellVolumes()
            /nAveragingSteps;

        totalRhoN_sum += allSpeciesRhoN[rN].internalField();
    }

    totalRhoN.internalField() = totalRhoN_sum;


    /*-----------------------------------------------------------------------*\
        Mass density
    \*-----------------------------------------------------------------------*/

    scalarField totalRhoM_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesRhoM, rM)
    {
        allSpeciesRhoM[rM].internalField() =
            allSpeciesM_RU[rM]
            /mesh.cellVolumes()
            /nAveragingSteps;

        totalRhoM_sum += allSpeciesRhoM[rM].internalField();
    }

    totalRhoM.internalField() = totalRhoM_sum;

    /*-----------------------------------------------------------------------*\
        Bulk velocity
    \*-----------------------------------------------------------------------*/

    vectorField totalMomentum_sum(mesh.nCells(), vector::zero);

    scalarField totalMass_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesVelocity, v)
    {
        // A check for 1/0 molecules is required.

        vectorField& singleSpeciesVelocity
        (
            allSpeciesVelocity[v].internalField()
        );

        forAll(singleSpeciesVelocity, sSV)
        {
            if (allSpeciesN_RU[v][sSV])
            {
                singleSpeciesVelocity[sSV] =
                    allSpeciesVelocitySum_RU[v][sSV]
                    /allSpeciesN_RU[v][sSV];

                totalMomentum_sum[sSV] +=
                    allSpeciesM_RU[v][sSV]
                    /allSpeciesN_RU[v][sSV]
                    *allSpeciesVelocitySum_RU[v][sSV];

                totalMass_sum[sSV] += allSpeciesM_RU[v][sSV];
            }
            else
            {
                singleSpeciesVelocity[sSV] = vector::zero;
            }
        }
    }

    volVectorField::InternalField& itotalVelocity =
        totalVelocity.internalField();

    forAll(itotalVelocity, tV)
    {
        if (totalMass_sum[tV] > VSMALL)
        {
            itotalVelocity[tV] = totalMomentum_sum[tV]/totalMass_sum[tV];
        }
        else
        {
            itotalVelocity[tV] = vector::zero;
        }
    }

    /*-----------------------------------------------------------------------*\
        Kinetic temperature
    \*-----------------------------------------------------------------------*/

    scalarField totalTemperatureVTerms_sum(mesh.nCells(), 0.0);

    scalarField totalN_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesTemperature, t)
    {
        // A check for 1/0 molecules is required.

        scalarField& singleSpeciesTemp
        (
            allSpeciesTemperature[t].internalField()
        );

        forAll(singleSpeciesTemp, sST)
        {
            if (allSpeciesN_RU[t][sST])
            {
                singleSpeciesTemp[sST] =
                    allSpeciesM_RU[t][sST]
                    /allSpeciesN_RU[t][sST]
                    /(3.0 * moleculeCloud::kb * allSpeciesN_RU[t][sST])
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[t][sST]
                        -
                        (
                            allSpeciesVelocitySum_RU[t][sST]
                            &
                            allSpeciesVelocitySum_RU[t][sST]
                        )
                        /allSpeciesN_RU[t][sST]
                    );

                totalTemperatureVTerms_sum[sST] +=
                    allSpeciesM_RU[t][sST]
                   /allSpeciesN_RU[t][sST]
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[t][sST]
                      -
                        (
                            allSpeciesVelocitySum_RU[t][sST]
                            &
                            allSpeciesVelocitySum_RU[t][sST]
                        )
                        /allSpeciesN_RU[t][sST]
                    );

                totalN_sum[sST] += allSpeciesN_RU[t][sST];
            }
            else
            {
                singleSpeciesTemp[sST] = 0.0;
            }
        }
    }

    volScalarField::InternalField& itotalTemperature =
        totalTemperature.internalField();

    forAll(itotalTemperature, tT)
    {
        if (totalN_sum[tT] > 0)
        {
            itotalTemperature[tT] =
                totalTemperatureVTerms_sum[tT]
                /(3.0 * moleculeCloud::kb * totalN_sum[tT]);
        }
        else
        {
            itotalTemperature[tT] = 0.0;
        }
    }

    /*-----------------------------------------------------------------------*\
        Mean kinetic energy
    \*-----------------------------------------------------------------------*/

    scalarField totalKE_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesMeanKE, mKE)
    {
        // A check for 1/0 molecules is required.

        scalarField& singleSpeciesMeanKE
        (
            allSpeciesMeanKE[mKE].internalField()
        );

        forAll(singleSpeciesMeanKE, sSMKE)
        {
            if (allSpeciesN_RU[mKE][sSMKE])
            {
                singleSpeciesMeanKE[sSMKE] =
                    allSpeciesM_RU[mKE][sSMKE]
                   /allSpeciesN_RU[mKE][sSMKE]
                   /(2.0*allSpeciesN_RU[mKE][sSMKE])
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[mKE][sSMKE]
                    );

                totalKE_sum[sSMKE] +=
                    allSpeciesM_RU[mKE][sSMKE]
                    /allSpeciesN_RU[mKE][sSMKE]
                    /2.0
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[mKE][sSMKE]
                    );
            }
            else
            {
                singleSpeciesMeanKE[sSMKE] = 0.0;
            }
        }
    }

    volScalarField::InternalField& itotalMeanKE = totalMeanKE.internalField();

    forAll(itotalMeanKE, tMKE)
    {
        if (totalN_sum[tMKE] > 0)
        {
            itotalMeanKE[tMKE] =
                totalKE_sum[tMKE]
                /totalN_sum[tMKE];
        }
        else
        {
            itotalMeanKE[tMKE] = 0.0;
        }
    }
}
