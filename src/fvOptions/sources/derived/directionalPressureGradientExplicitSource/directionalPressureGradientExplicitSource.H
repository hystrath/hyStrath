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
    Foam::fv::directionalPressureGradientExplicitSource

Group
    grpFvOptionsSources

Description
    Creates an explicit  pressure gradient source in such a way to deflect the
    flow towards an specific direction (flowDir). Alternatively add an extra
    pressure drop in the flowDir direction using a model.

    \heading Source usage
    Example usage:
    \verbatim
    airDeflection
    {
        type            directionalPressureGradientExplicitSource;
        active          true;

        directionalPressureGradientExplicitSourceCoeffs
        {
            selectionMode   cellZone;
            cellZone        cZone;

            fields      (U);            // Name of the field
            flowDir     (1 1 0);        // Desired flow direction
            faceZone    f0Zone;         // Face zone upstream cell zone
            relaxationFactor    0.3;    // Relaxation factor for flow
                                        // deflection (default 0.3)

            //Pressure drop model [Pa]
            model       volumetricFlowRateTable;//constant;//DarcyForchheimer;

            //DarcyForchheimer model
            // deltaP = (D*mu + 0.5*rho*magUn)*magUn*length_

            D           5e7;
            I           0;
            length      1e-3;

            //constant model
            pressureDrop    40;

            //volumetricFlowRateTable model
            outOfBounds     clamp;
            fileName        "volFlowRateTable";
        }
    }
    \endverbatim

    NOTE: In order to obtain the upwind velocities this function loops over
    the slaves cells of the faceZone specified in the dictionary, on the other
    hand, the cellZone to which this source term is applied should be composed
    of the master cells and they should be 'downwind' the faceZone.

SourceFiles
    directionalPressureGradientExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef directionalPressureGradientExplicitSource_H
#define directionalPressureGradientExplicitSource_H

#include "cellSetOption.H"
#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvOption.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
          Class directionalPressureGradientExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class directionalPressureGradientExplicitSource
:
    public cellSetOption
{
public:

    //  Public enumeration

        //- Modes of pressure drop
        enum pressureDropModel
        {
            pVolumetricFlowRateTable,
            pConstant,
            pDarcyForchheimer
        };


private:

    // Private data

        static const NamedEnum<pressureDropModel, 3> PressureDropModelNames_;

         //- Pressure drop model
        pressureDropModel model_;

        //- Pressure gradient before correction
        vectorField gradP0_;

        //- Change in pressure gradient
        vectorField dGradP_;

        //- Pressure drop due to porous media
        vectorField gradPporous_;

        //- Flow direction
        vector flowDir_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> invAPtr_;

        //- Darcy pressure loss coefficient
        scalar D_;

        //- Inertia pressure lost coefficient
        scalar I_;

        //- Porous media length
        scalar length_;

        //- Constant pressure drop
        scalar pressureDrop_;

        //- Volumetric flow rate vs pressure drop table
        interpolationTable<scalar> flowRate_;

        //- Name of the faceZone at the heat exchange inlet
        word faceZoneName_;

        //- Id for the face zone
        label zoneID_;

         //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- Relaxation factor
        scalar relaxationFactor_;

        //- Cells faces mapping
        labelList cellFaceMap_;


    // Private Member Functions

        //- Init
        void initialise();

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const vectorField& gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);

        //- Disallow default bitwise copy construct
        directionalPressureGradientExplicitSource
        (
            const directionalPressureGradientExplicitSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const directionalPressureGradientExplicitSource&);


public:

    //- Runtime type information
    TypeName("directionalPressureGradientExplicitSource");


    // Constructors

        //- Construct from explicit source name and mesh
        directionalPressureGradientExplicitSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
