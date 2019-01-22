/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "Boltzmann.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(Boltzmann, 0);

        addToRunTimeSelectionTable
        (
            conductivityModel,
            Boltzmann,
            mhdModel
        );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Boltzmann::Boltzmann
(
    const mhdModel& dict,
    const fvMesh& mesh
)
:
    conductivityModel(dict, mesh),
    crossSections_(dict.lookupOrDefault<word>("crossSection", "calculate"))
{
    Info << "Creating Boltzmann conductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Boltzmann::~Boltzmann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volScalarField Boltzmann::sigma() const
{
    volScalarField nDe = dict_.thermo().composition().nD("e-");
    Info << max(nDe) << endl;
    Info << min(nDe) << endl;
    volScalarField T = dict_.thermo().Tt();
    //Setting up physical constants
    const dimensionedScalar electron_mass = dimensionedScalar("electron_mass", dimensionSet(1, 0, 0, 0, 0, 0, 0), 9.10938356e-31);
	const dimensionedScalar k_Boltzmann = Foam::constant::physicoChemical::k;
	const dimensionedScalar electron_charge = dimensionedScalar("electron_mass", dimensionSet(0, 0, 1, 0, 0, 1, 0), 1.60217662e-19);
	const scalar NA = 6e+23;
	const dimensionedScalar Qen = dimensionedScalar("Qen", dimensionSet(0, 2, 0, 0, 0, 0, 0), 4e-20);

    Info << "Calculating electron velocity" << endl;
    volScalarField ve
    (
        IOobject
        (
            "ve",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(3*k_Boltzmann*T/electron_mass)
    );
    Info << max(ve) << endl;
    Info << "Calculating Coulomb logarithm" << endl;
    volScalarField Coulomb_log
    (
        IOobject
        (
            "Coulomb_log",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "Coulomb_log",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );
    forAll(Coulomb_log, cellI)
    {
        if (nDe[cellI] == 0.0)
        {
            Coulomb_log[cellI] = 0.0;
        }
        else
        {
            Coulomb_log[cellI] = log(1.059e+22*pow(T[cellI], 1.5)*pow(nDe[cellI], -0.5));
        }
    }

    Info << "Calculating number density of neutrals" << endl;
    volScalarField nn
    (
        IOobject
        (
            "nn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "nn",
            dimensionSet(0, -3, 0, 0, 0, 0, 0),
            0.0
        )
    );

    forAll(dict_.thermo().composition().species(), i)
    {
        word name = dict_.thermo().composition().species()[i];
        nn += dict_.thermo().composition().nD(name);
    }
    Info << max(nn) << endl;
	//Calculating each particular Qen

    scalar W;
	int i;
	scalar collisionIntegralNeutral;
	FixedList<scalar, 4> piOmegaCoeffs;
	scalar cfA, cfB, cfC, cfD;

// * * * * * * * * * * * * * * * * Reading species data from dictionaries * * * * * * * * * * * * * * * * 
	const dictionary thermoDEM =
        (
        	IFstream
        	(
        	    fileName(dict_.thermo().lookup("foamChemistryThermoFile")).expand()
       		 )()
    	);

	const scalar W_el = readScalar(thermoDEM.subDict("e-").subDict("specie").lookup("molWeight"))*1.0e-3;
		 Info << "DEBUG: W_el = "<< W_el << endl;
        const scalar Runi = Foam::constant::physicoChemical::R.value();
		 Info << "DEBUG: Runi = "<< Runi << endl;


		volScalarField collision_frequency
    		(
        		IOobject
        		(
        		    "collision_frequency",
                     mesh_.time().timeName(),
                     mesh_,
        		    IOobject::NO_READ,
        		    IOobject::AUTO_WRITE
        		),
			mesh_,
			dimensionedScalar("collission_frequency",dimensionSet(0,0,-1,0,0,0, 0),0)
    		);

		const dictionary transportDict =
        (
        	IFstream
        	(
        	    fileName(dict_.thermo().lookup("transportPropertiesFile")).expand()
       		 )()
    	);
// * * * * * * * * * * * * * * * * Reading collision data * * * * * * * * * * * * * * * *
// - Allocating memory for electron-species k collision frequency
    volScalarField nu_k //collision frequency for specie k
    (
        IOobject
        (
            "nu_k",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("nu_k", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    );
    //
    FixedList<scalar,4> defaultList;
    forAll(defaultList, i)
    {
        defaultList[i] = 0.0;
    }
    
    forAll(dict_.thermo().composition().species(), i)
    {
        word name = dict_.thermo().composition().species()[i];
        Info << name << endl;
        if (crossSections_ == "calculate")
        {
            if (name != "e-")
            {
                if (transportDict.subDict("collisionData").subDict("neutralNeutralInteractions").subDict("Gupta1990O")
                .subDict("Omega11").found("e-_"+name))
                {
                    piOmegaCoeffs = transportDict.subDict("collisionData").subDict("neutralNeutralInteractions")
                    .subDict("Gupta1990O").subDict("Omega11").lookupOrDefault<FixedList<scalar, 4>>("e-_"+name, defaultList);    
                }
                else if (transportDict.subDict("collisionData").subDict("neutralNeutralInteractions")
                .subDict("Gupta1990O").subDict("Omega11").found(name+"_e-"))
                {
                    piOmegaCoeffs = transportDict.subDict("collisionData").subDict("neutralNeutralInteractions")
                    .subDict("Gupta1990O").subDict("Omega11").lookupOrDefault<FixedList<scalar, 4>>(name+"_e-", defaultList);    
                }
                else
                {
                    FatalErrorIn("Electric conductivity model: ")
                    << "Collision integral data missing for electrons and" << name << "."
                    << exit(FatalError);
                }
            
                Info << "Calculating piOmegaCoeffs" << endl;
                cfA = piOmegaCoeffs[0];
                cfB = piOmegaCoeffs[1];
                cfC = piOmegaCoeffs[2];
                cfD = piOmegaCoeffs[3];
                Info << piOmegaCoeffs << endl;
                Info << "Reading molWeight of " << name << endl;
                W = readScalar(thermoDEM.subDict(name).subDict("specie").lookup("molWeight"));
                Info << W << endl;

                volScalarField piOmega
                (
                    IOobject
                    (
                        "piOmega",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "piOmega",
                        dimensionSet(0, 2, 0, 0, 0, 0, 0),
                        0.0
                    )
                );

                forAll(piOmega, cellI)
                {
                    piOmega[cellI] = 1e-20*Foam::exp(cfD)*pow(T[cellI], cfA*log(T[cellI])*log(T[cellI]) + cfB*log(T[cellI])+cfC);
                }
                Info << max(piOmega) << endl;
                Info << min(piOmega) << endl;
                nu_k = dict_.thermo().composition().nD(name)*piOmega*ve;
            }
            Info << "Adding to total collision frequency" << endl;
            collision_frequency = collision_frequency+nu_k; 
        }
        if (crossSections_ == "BityurinBocharov")
        {
            if (readLabel(thermoDEM.subDict(name).subDict("specie").lookup("charge")) == 1)
            {
                Info<< "Species" << name <<" is an ion" <<endl;
                volScalarField Qei //collision frequency for specie k
                (
                    IOobject
                    (
                        "Qei",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("Qei", dimensionSet(0, 2, 0, 0, 0, 0, 0), 0.0)
                );
                forAll(Qei, cellI)
                {
                    Qei[cellI] = 1.3962634016*pow(electron_charge.value(), 4)*Coulomb_log[cellI]/(k_Boltzmann.value()*k_Boltzmann.value()*
                    T[cellI]*T[cellI]);
                }
                Info << max(Qei) << endl;
                Info << min(Qei) << endl;
                nu_k = dict_.thermo().composition().nD(name)*Qei*ve;
            }
            if (readLabel(thermoDEM.subDict(name).subDict("specie").lookup("charge")) == 0)
            {
                Info<< "Species" << name <<" is a neutral" <<endl;
                nu_k = dict_.thermo().composition().nD(name)*Qen*ve;
            }
            Info << "Adding to total collision frequency" << endl;
            collision_frequency = collision_frequency+nu_k; 
        }


		Info << max(nu_k) << endl;
		Info << min(nu_k) << endl;
		Info << min(collision_frequency) << endl;
		Info << max(collision_frequency) << endl;		
	}
    
    volScalarField taue
    (
        IOobject
        (
            "taue",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("taue", dimensionSet(0, 0, 1, 0, 0, 0, 0), 0.0)
    );
    forAll(taue, cellI)
    {
        if(collision_frequency[cellI] == 0)
        {
            taue[cellI] = 0.0;
        }
        else
        {
            taue[cellI] = pow(collision_frequency[cellI], -1);
        } 
    }
    Info << min(taue) << endl;
    Info << max(taue) << endl;


// * * * * * * * * * * * * * * * * Allocating sigma * * * * * * * * * * * * * * * * 
    	volScalarField sigma
    	(
            IOobject
            (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
        	IOobject::NO_READ,
        	IOobject::AUTO_WRITE
            ),
	    electron_charge*electron_charge/electron_mass*taue*nDe
    	);
	Info << min(sigma) << endl;
	Info << max(sigma) << endl;
    return sigma;
}


}
}

// ************************************************************************* //
