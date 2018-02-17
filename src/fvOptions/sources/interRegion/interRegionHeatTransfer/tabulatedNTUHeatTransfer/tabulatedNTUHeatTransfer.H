/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

Class
    Foam::fv::tabulatedNTUHeatTransfer

Group
    grpFvOptionsSources

Description
    Tabulated heat transfer model.

    The heat flux is calculated based on the Number of Thermal Units (NTU).
    A 2-D table of NTU as functions of the primary and secondary mass flow
    rates is required.

    The exchanger geometry can be specified using either:
    - \c calculated:
        - inlet area of each region and core volume determined by interrogating
          mesh patches, and mesh-to-mesh interpolation volume
        - calculated core volume canbe partially blocked by specifying a
          \c coreBlockageRatio [0-1] entry
    - \c user:
        - inlet area of each region privided by the user
        - core volume automatically calculated by the mesh-to-mesh
          interpolation volume if not provided by user

    Heat transfer coefficient calculated by:

    \f[
        htc = C_{min} \frac{NTU}{V_{core}}
    \f]

    Where \f$ C_{min} \f$ is given by:

    \f[
        C_{min} = min \left(Cp_1 \dot{m}_1, Cp_2 \dot{m}_2 \right)
    \f]

    \heading Example usage

    \verbatim
    coolerToAir
    {
        type            tabulatedNTUHeatTransfer;
        active          yes;

        tabulatedNTUHeatTransferCoeffs
        {
            interpolationMethod cellVolumeWeight;
            nbrRegion       air;
            master          true;

            fields          (h);
            outOfBounds     clamp;
            fileName        "ntuTable";
            nbrModel        airToCooler;
            semiImplicit    no;


            geometryMode    user;
            Ain             0.01728;
            AinNbr          0.3456;
            Vcore           0.01244;  // Optional
    
            // geometryMode    calculated;
            // inletPatch      inlet_HWK;
            // inletPatchNbr   inlet_air;
            // inletBlockageRatio 0.10;
            // inletBlockageRatioNbr 0.04;
            // coreBlockageRatio 0;
        }
    }

SourceFiles
    tabulatedNTUHeatTransfer.C

SeeAlso
    interRegionHeatTransferModel.H

\*---------------------------------------------------------------------------*/

#ifndef tabulatedNTUHeatTransfer_H
#define tabulatedNTUHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"
#include "NamedEnum.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class tabulatedNTUHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedNTUHeatTransfer
:
    public interRegionHeatTransferModel
{
public:

    // Public enumerations

        enum geometryModeType
        {
            gmCalculated,
            gmUser
        };

        static const NamedEnum<geometryModeType, 2> geometryModelNames_;


private:

    // Private data

        //- Name of velocity field; default = U
        word UName_;

        //- Name of neighbour velocity field; default = U
        word UNbrName_;

        //- Name of density field; default = rho
        word rhoName_;

        //- Name of neighbour density field; default = rho
        word rhoNbrName_;

        //- Pointer to 2-D look-up table of NTU f(mDot1, mDot2)
        autoPtr<interpolation2DTable<scalar>> ntuTable_;

        //- Geometry input mode
        geometryModeType geometryMode_;

        //- Inlet area [m2]
        scalar Ain_;

        //- Neighbour region inlet area [m2]
        scalar AinNbr_;

        //- Heat exchanger core volume
        scalar Vcore_;


    // Private Member functions

        //- NTU table helper
        const interpolation2DTable<Foam::scalar>& ntuTable();

        //- Thermophysical properties helper
        const basicThermo& thermo(const fvMesh& mesh) const;

        //- Initialise geometry
        void initialiseGeometry();


public:

    //- Runtime type information
    TypeName("tabulatedNTUHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedNTUHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~tabulatedNTUHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc();


        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
