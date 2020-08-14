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
            electricalConductivityModel,
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
    electricalConductivityModel(dict, mesh),
    crossSections_(dict.lookupOrDefault<word>("crossSection", "calculate")) // TODO this needs to disappear: redondant
{
    Info << "Loading the Boltzmann electricalConductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Boltzmann::~Boltzmann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Boltzmann::update()
{
    // NB: most of the operations in this function needs to be dealt with
    //     pre-computation in the constructor
    
    const volScalarField& nDe = dict_.thermo().composition().nD("e-");
    Info << max(nDe) << endl;
    Info << min(nDe) << endl;
    const volScalarField& T = dict_.thermo().T();
    
    const dimensionedScalar electron_mass =
        dimensionedScalar
        (
            "electron_mass",
            dimensionSet(1, 0, 0, 0, 0, 0, 0),
            9.10938356e-31
        );
	  const dimensionedScalar kB = constant::physicoChemical::k;
	  const dimensionedScalar electron_charge = constant::electromagnetic::e;
	  const scalar NA = constant::physicoChemical::NA.value();
	  const dimensionedScalar Qen =
	      dimensionedScalar
	      (
	          "Qen",
	          dimensionSet(0, 2, 0, 0, 0, 0, 0),
	          4e-20
        );
        
    //- Reading thermoDEM dictionary
    const dictionary thermoDEM =
    (
      	IFstream
      	(
      	    fileName(dict_.thermo().lookup("foamChemistryThermoFile")).expand()
     		 )()
  	);

    const scalar W_el =
        readScalar
        (
            thermoDEM.subDict("e-").subDict("specie").lookup("molWeight")
        )*1.0e-3;

    const scalar Runi = Foam::constant::physicoChemical::R.value();
    
    //- Reading transportProperties dictionary
    const dictionary transportDict =
    (
	      IFstream
	      (
	          fileName(dict_.thermo().lookup("transportPropertiesFile")).expand()
	       )()
    );
		 

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
        sqrt(3.0*kB*T/electron_mass)
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
            dimless,
            0.0
        )
    );
    
    forAll(Coulomb_log, cellI)
    {
        if (nDe[cellI] > SMALL)
        {
            Coulomb_log[cellI] =
                log(1.059e+22*pow(T[cellI], 1.5)*pow(nDe[cellI], -0.5));
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
        word speciesname = dict_.thermo().composition().species()[i];
        nn += dict_.thermo().composition().nD(speciesname);
    }
    Info << max(nn) << endl;
    
	  //- Calculating each Qen
    scalar W = 0.0;
	  label i = 0;
	  scalar collisionIntegralNeutral = 0.0;
	  FixedList<scalar, 4> piOmegaCoeffs;
	  scalar cfA, cfB, cfC, cfD;

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
			  dimensionedScalar("collision_frequency", dimless/dimTime, 0.0)
		);


    // Reading collision data
    // Allocating memory for electron-species k collision frequency
    
    //- TODO for specie k
    volScalarField Qei 
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
        dimensionedScalar
        (
            "Qei", dimensionSet(0, 2, 0, 0, 0, 0, 0), 0.0
        )
    );
                
    //- Collision frequency for specie k
    volScalarField nu_k
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
        dimensionedScalar("nu_k", dimless/dimTime, 0.0)
    );

    FixedList<scalar,4> defaultList;
    forAll(defaultList, i)
    {
        defaultList[i] = 0.0;
    }

    forAll(dict_.thermo().composition().species(), i)
    {
        word speciesname = dict_.thermo().composition().species()[i];
        Info << speciesname << endl;
        
        if (crossSections_ == "calculate")
        {
            if (speciesname != "e-")
            {
                if 
                (
                    transportDict.subDict("collisionData")
                        .subDict("tabulatedInteractions").subDict("Gupta1990O")
                        .subDict("Omega11").found("e-_" + speciesname)
                )
                {
                    piOmegaCoeffs =
                        transportDict.subDict("collisionData")
                        .subDict("tabulatedInteractions")
                        .subDict("Gupta1990O").subDict("Omega11")
                        .lookupOrDefault<FixedList<scalar, 4>>
                        (
                            "e-_" + speciesname, defaultList
                        );
                }
                else if
                (
                    transportDict.subDict("collisionData")
                        .subDict("tabulatedInteractions").subDict("Gupta1990O")
                        .subDict("Omega11").found(speciesname + "_e-"))
                {
                    piOmegaCoeffs =
                        transportDict.subDict("collisionData")
                            .subDict("tabulatedInteractions")
                            .subDict("Gupta1990O").subDict("Omega11")
                            .lookupOrDefault<FixedList<scalar, 4>>
                            (
                                speciesname + "_e-",
                                defaultList
                            );
                }
                else
                {
                    FatalErrorIn("Electric electricalConductivity model: ")
                        << "Collision integral data missing for electrons and"
                        << speciesname << "." << exit(FatalError);
                }

                Info << "Calculating piOmegaCoeffs" << endl;
                cfA = piOmegaCoeffs[0];
                cfB = piOmegaCoeffs[1];
                cfC = piOmegaCoeffs[2];
                cfD = piOmegaCoeffs[3];
                Info << piOmegaCoeffs << endl;
                W =
                    readScalar
                    (
                        thermoDEM.subDict(speciesname).subDict("specie")
                            .lookup("molWeight")
                    );

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
                    piOmega[cellI] = 1e-20*Foam::exp(cfD)
                       *pow
                        (
                            T[cellI],
                            cfA*log(T[cellI])*log(T[cellI]) 
                              + cfB*log(T[cellI])+cfC
                        );
                }
                Info << gMax(piOmega) << endl;
                Info << gMin(piOmega) << endl;
                nu_k = dict_.thermo().composition().nD(speciesname)*piOmega*ve;
            }
            
            collision_frequency = collision_frequency + nu_k;
        }
        
        if (crossSections_ == "BityurinBocharov")
        {
            const label speciescharge =
                readLabel
                (
                    thermoDEM.subDict(speciesname).subDict("specie")
                        .lookup("charge")
                );
            
            if (speciescharge == 1)
            {
                forAll(Qei, cellI)
                {
                    Qei[cellI] = 1.3962634016*pow(electron_charge.value(), 4)
                        *Coulomb_log[cellI]/pow(kB.value()*T[cellI], 2);
                }
                
                Info << gMax(Qei) << endl;
                Info << gMin(Qei) << endl;
                
                nu_k = dict_.thermo().composition().nD(speciesname)*Qei*ve;
            }
            else if (speciescharge == 0)
            {
                nu_k = dict_.thermo().composition().nD(speciesname)*Qen*ve;
            }
            
            collision_frequency = collision_frequency + nu_k;
        }

		    Info << gMax(nu_k) << endl;
		    Info << gMin(nu_k) << endl;
		    Info << gMin(collision_frequency) << endl;
		    Info << gMax(collision_frequency) << endl;
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
        dimensionedScalar("taue", dimTime, 0.0)
    );
    
    forAll(taue, cellI)
    {
        if(collision_frequency[cellI] != 0)
        {
            taue[cellI] = 1.0/collision_frequency[cellI];
        }
    }
    Info << gMin(taue) << endl;
    Info << gMax(taue) << endl;

    sigma_.primitiveFieldRef() = pow(electron_charge, 2)/electron_mass*taue*nDe;
    
    Info << gMin(sigma) << endl;
    Info << gMax(sigma) << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //End namespace mhd
} //End namespace Foam

// ************************************************************************* //
