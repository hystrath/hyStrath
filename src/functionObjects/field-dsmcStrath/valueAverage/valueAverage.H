/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::valueAverage

Group
    grpFieldFunctionObjects

Description
    This function object calculates the average value from the output of
    function objects that generate singular values.

Usage
    Example of function object specification:
    \verbatim
    valueAverage1
    {
        type        valueAverage;
        libs        ("libfieldFunctionObjects.so");
        ...
        writeToFile yes;
        log         yes;
        functionObject forceCoeffs1;
        fields      (Cm Cd Cl);
        window      0.5;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: valueAverage | yes         |
        writeToFile  | write average data to file |  no      | yes
        log          | write average data to standard output | no | yes
        fields       | list of fields to process | yes       |
        functionObject | Name of function object to retrieve data from | yes |
        resetOnRestart | Reset the averaging on restart | no | no |
    \endtable

    Output data is written to the file \<timeDir\>/valueAverage.dat

See also
    Foam::functionObject
    Foam::functionObjects::stateFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    valueAverage.C
    valueAverageTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_valueAverage_H
#define functionObjects_valueAverage_H

#include "regionFunctionObject.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class valueAverage Declaration
\*---------------------------------------------------------------------------*/

class valueAverage
:
    public regionFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of fields on which to operate
        wordList fieldNames_;

        //- Averaging window
        scalar window_;

        //- Average time per field
        List<scalar> totalTime_;

        //- Reset the averaging process on restart
        Switch resetOnRestart_;


    // Protected Member Functions

        //- Templated function to calculate the average
        template<class Type>
        void calc
        (
            const word& fieldName,
            const word& meanName,
            const scalar alpha,
            const scalar beta,
            bool& processed
        );

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        valueAverage(const valueAverage&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const valueAverage&) = delete;


public:

    //- Runtime type information
    TypeName("valueAverage");

    //- Constructor
    valueAverage
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~valueAverage();


    // Public Member Functions

        //- Read the field value average data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write the value average
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valueAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
