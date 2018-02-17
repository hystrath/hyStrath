/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

Class
    Foam::DEShybrid

Description
    Hybrid convection scheme of Travin et al. for hybrid RAS/LES calculations.

    The scheme provides a blend between two convection schemes, based on local
    properties including the wall distance, velocity gradient and eddy
    viscosity.  The scheme was originally developed for DES calculations to
    blend a low-dissipative scheme, e.g. linear, in the vorticity-dominated,
    finely-resolved regions and a numerically more robust, e.g. upwind-biased,
    convection scheme in irrotational or coarsely-resolved regions.

    The routine calculates the blending factor denoted as "sigma" in the
    literature reference, where 0 <= sigma <= sigmaMax, which is then employed
    to set the weights:
    \f[
        weight = (1-sigma) w_1 + sigma w_2
    \f]

    where
    \vartable
        sigma | blending factor
        w_1   | scheme 1 weights
        w_2   | scheme 2 weights
    \endvartable

    First published in:
    \verbatim
        A. Travin, M. Shur, M. Strelets, P. Spalart (2000).
        Physical and numerical upgrades in the detached-eddy simulation of
        complex turbulent flows.
        In Proceedings of the 412th Euromech Colloquium on LES and Complex
        Transition and Turbulent Flows, Munich, Germany
    \endverbatim

    Original publication contained a typo for C_H3 constant. Corrected version
    with minor changes for 2 lower limiters published in:
    \verbatim
        P. Spalart, M. Shur, M. Strelets, A. Travin (2012).
        Sensitivity of Landing-Gear Noise Predictions by Large-Eddy
        Simulation to Numerics and Resolution.
        AIAA Paper 2012-1174, 50th AIAA Aerospace Sciences Meeting,
        Nashville / TN, Jan. 2012
    \endverbatim

    Example of the DEShybrid scheme specification using linear within the LES
    region and linearUpwind within the RAS region:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,U)      Gauss DEShybrid
            linear                        // scheme 1
            linearUpwind grad(U)          // scheme 2
            0.65                          // DES coefficient, typically = 0.65
            30                            // Reference velocity scale
            2                             // Reference length scale
            0                             // Minimum sigma limit (0-1)
            1                             // Maximum sigma limit (0-1)
            1.0e-03;                      // Limiter of B function, typically 1.0e-03
        .
        .
    }
    \endverbatim

Notes
  - Scheme 1 should be linear (or other low-dissipative schemes) which will
    be used in the detached/vortex shedding regions.
  - Scheme 2 should be an upwind/deferred correction/TVD scheme which will
    be used in the free-stream/Euler/boundary layer regions.
  - the scheme is compiled into a separate library, and not available to
    solvers by default.  In order to use the scheme, add the library as a
    run-time loaded library in the \$FOAM\_CASE/system/controlDict
    dictionary, e.g.:

        libs ("libturbulenceModelSchemes.so");

SourceFiles
    DEShybrid.C

\*---------------------------------------------------------------------------*/

#ifndef DEShybrid_H
#define DEShybrid_H

#include "surfaceInterpolationScheme.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "blendedSchemeBase.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          class DEShybrid Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DEShybrid
:
    public surfaceInterpolationScheme<Type>,
    public blendedSchemeBase<Type>
{
    // Private Data

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type>> tScheme2_;

        //- DES Coefficient
        scalar CDES_;

        //- Reference velocity scale [m/s]
        dimensionedScalar U0_;

        //- Reference length scale [m]
        dimensionedScalar L0_;

        //- Minimum bound for sigma (0 <= sigmaMin <= 1)
        scalar sigmaMin_;

        //- Maximum bound for sigma (0 <= sigmaMax <= 1)
        scalar sigmaMax_;

        //- Limiter of B function
        scalar OmegaLim_;

        //- Scheme constants
        scalar CH1_;
        scalar CH2_;
        scalar CH3_;

        //- Disallow default bitwise copy construct
        DEShybrid(const DEShybrid&);

        //- Disallow default bitwise assignment
        void operator=(const DEShybrid&);


    // Private Member Functions

        //- Calculate the blending factor
        tmp<surfaceScalarField> calcBlendingFactor
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const volScalarField& nuEff,
            const volVectorField& U,
            const volScalarField& delta
        ) const
        {
            tmp<volTensorField> gradU(fvc::grad(U));
            const volScalarField S(sqrt(2.0)*mag(symm(gradU())));
            const volScalarField Omega(sqrt(2.0)*mag(skew(gradU())));
            const dimensionedScalar tau0_ = L0_/U0_;

            const volScalarField B
            (
                CH3_*Omega*max(S, Omega)
               /max(0.5*(sqr(S) + sqr(Omega)), sqr(OmegaLim_/tau0_))
            );
            const volScalarField K
            (
                max(Foam::sqrt(0.5*(sqr(S) + sqr(Omega))), 0.1/tau0_)
            );

            const volScalarField lTurb(Foam::sqrt(nuEff/(pow(0.09, 1.5)*K)));
            const volScalarField g(tanh(pow4(B)));
            const volScalarField A
            (
                 CH2_*max(scalar(0), CDES_*delta/max(lTurb*g, 1.0e-15*L0_) - 0.5)
            );

            const volScalarField factor
            (
                IOobject
                (
                    typeName + ":Factor",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                max(sigmaMax_*tanh(pow(A, CH1_)), sigmaMin_)
            );

            if (blendedSchemeBaseName::debug)
            {
                factor.write();
            }

            return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    vf.name() + "BlendingFactor",
                    fvc::interpolate(factor)
                )
            );
        }


public:

    //- Runtime type information
    TypeName("DEShybrid");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        DEShybrid(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            CDES_(readScalar(is)),
            U0_("U0", dimLength/dimTime, readScalar(is)),
            L0_("L0", dimLength, readScalar(is)),
            sigmaMin_(readScalar(is)),
            sigmaMax_(readScalar(is)),
            OmegaLim_(readScalar(is)),
            CH1_(3.0),
            CH2_(1.0),
            CH3_(2.0)
        {
            if (U0_.value() <= 0)
            {
                FatalErrorInFunction
                    << "U0 coefficient must be greater than 0. "
                    << "Current value: " << U0_ << exit(FatalError);
            }
            if (L0_.value() <= 0)
            {
                FatalErrorInFunction
                    << "L0 coefficient must be greater than 0. "
                    << "Current value: " << L0_ << exit(FatalError);
            }
            if (sigmaMin_ < 0)
            {
                FatalErrorInFunction
                    << "sigmaMin coefficient must be greater than or equal to 0. "
                    << "Current value: " << sigmaMin_ << exit(FatalError);
            }
            if (sigmaMax_ < 0)
            {
                FatalErrorInFunction
                    << "sigmaMax coefficient must be greater than or equal to 0. "
                    << "Current value: " << sigmaMax_ << exit(FatalError);
            }
            if (sigmaMin_ > 1)
            {
                FatalErrorInFunction
                    << "sigmaMin coefficient must be less than or equal to 1. "
                    << "Current value: " << sigmaMin_ << exit(FatalError);
            }
            if (sigmaMax_ > 1)
            {
                FatalErrorInFunction
                    << "sigmaMax coefficient must be less than or equal to 1. "
                    << "Current value: " << sigmaMax_ << exit(FatalError);
            }
        }

        //- Construct from mesh, faceFlux and Istream
        DEShybrid
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            CDES_(readScalar(is)),
            U0_("U0", dimLength/dimTime, readScalar(is)),
            L0_("L0", dimLength, readScalar(is)),
            sigmaMin_(readScalar(is)),
            sigmaMax_(readScalar(is)),
            OmegaLim_(readScalar(is)),
            CH1_(3.0),
            CH2_(1.0),
            CH3_(2.0)
        {
            if (U0_.value() <= 0)
            {
                FatalErrorInFunction
                    << "U0 coefficient must be greater than 0. "
                    << "Current value: " << U0_ << exit(FatalError);
            }
            if (L0_.value() <= 0)
            {
                FatalErrorInFunction
                    << "L0 coefficient must be greater than 0. "
                    << "Current value: " << U0_ << exit(FatalError);
            }
            if (sigmaMin_ < 0)
            {
                FatalErrorInFunction
                    << "sigmaMin coefficient must be greater than or equal to 0. "
                    << "Current value: " << sigmaMin_ << exit(FatalError);
            }
            if (sigmaMax_ < 0)
            {
                FatalErrorInFunction
                    << "sigmaMax coefficient must be greater than or equal to 0. "
                    << "Current value: " << sigmaMax_ << exit(FatalError);
            }
            if (sigmaMin_ > 1)
            {
                FatalErrorInFunction
                    << "sigmaMin coefficient must be less than or equal to 1. "
                    << "Current value: " << sigmaMin_ << exit(FatalError);
            }
            if (sigmaMax_ > 1)
            {
                FatalErrorInFunction
                    << "sigmaMax coefficient must be less than or equal to 1. "
                    << "Current value: " << sigmaMax_ << exit(FatalError);
            }
        }


    // Member Functions

        //- Return the face-based blending factor
        virtual tmp<surfaceScalarField> blendingFactor
        (
             const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            typedef compressible::turbulenceModel cmpModel;
            typedef incompressible::turbulenceModel icoModel;

            // Assuming that LES delta field is called "delta"
            const volScalarField& delta = this->mesh().template
                 lookupObject<const volScalarField>("delta");

            // Could avoid the compressible/incompressible case by looking
            // up all fields from the database - but retrieving from model
            // ensures consistent fields are being employed e.g. for multiphase
            // where group name is used

            if (mesh.foundObject<icoModel>(icoModel::propertiesName))
            {
                const icoModel& model =
                    mesh.lookupObject<icoModel>(icoModel::propertiesName);

                return calcBlendingFactor(vf, model.nuEff(), model.U(), delta);
            }
            else if (mesh.foundObject<cmpModel>(cmpModel::propertiesName))
            {
                const cmpModel& model =
                    mesh.lookupObject<cmpModel>(cmpModel::propertiesName);

                return calcBlendingFactor
                (
                    vf, model.muEff()/model.rho(), model.U(), delta
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Scheme requires a turbulence model to be present. "
                    << "Unable to retrieve turbulence model from the mesh "
                    << "database" << exit(FatalError);

                return tmp<surfaceScalarField>(nullptr);
            }
        }


        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor(vf));

            return
                (scalar(1) - bf)*tScheme1_().weights(vf)
              + bf*tScheme2_().weights(vf);
        }


        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor(vf));

            return
                (scalar(1) - bf)*tScheme1_().interpolate(vf)
              + bf*tScheme2_().interpolate(vf);
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme1_().corrected() || tScheme2_().corrected();
        }


        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField bf(blendingFactor(vf));

            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        (scalar(1) - bf)
                      * tScheme1_().correction(vf)
                      + bf
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        (scalar(1) - bf)
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return (bf*tScheme2_().correction(vf));
            }
            else
            {
                return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
                (
                    nullptr
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
