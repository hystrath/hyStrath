/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
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

Description

\*---------------------------------------------------------------------------*/

#include "integral.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //




//- scalar field, scalar field
integral::integral
(
    const scalar& binWidth,
    const scalarField& xData,
    const scalarField& yData,
    const word& integrationType
)
:
    area_(0.0)
{
    word type = integrationType;

    label noPanels = xData.size() -1;

    if(type == "simpsons")
    {
        // simpsons 1/3 rule - even no of panels
        if((noPanels % 2) == 0)
        {
            area_ += yData[0];
            area_ += yData[noPanels];
            
            for (label i=1; i<noPanels; i++)
            {
                if((i % 2) == 0) // even
                {
                    area_ += 2.0*yData[i];
                }
                else // odd
                {
                    area_ += 4.0*yData[i];
                }
            }
            
            area_ *= binWidth/3.0;

            Info << "Integration (simpsons 1/3) = " << area_  << endl;
        }

        else if((noPanels % 3) == 0) // simpsons 3/8 rule
        { 
            area_ += yData[0];
            area_ += yData[noPanels];

            label counter = 0;
            
            for (label i=1; i<noPanels; i++)
            {
                if(counter < 2) // even
                {
                    counter++;
                    area_ += 3.0*yData[i];
                }
                else // odd
                {
                    area_ += 2.0*yData[i];
                    counter = 0;
                }
            }
            
            area_ *= binWidth*3.0/8.0;

            Info << "Integration (simpsons 3/8) = " << area_  << endl;
        }

        else
        {

            type = "default";

            Info<< "Simpsons rule works on number of panels "
                <<  "that are even (Simpsons 1/3 rule), or divisible by 3 "
                <<  "(Simpsons 3/8 rule). "  << nl
                <<  "Number of panels in this example: " << noPanels << nl
                <<  "WARNING: Changing to the default type..." 
                << endl;

//             FatalErrorIn("integral::integral()")
//                 << "Integration type selected: " << type
//                 << " does not work on the no. of panels provided: " <<noPanels << nl
//                 << "Panels have to be even for 1/3 simpsons"  
//                 << " rule or divisible by 3 for simpsons 3/8 rule." << nl
//                 << exit(FatalError);
        }
    }

    if((type == "default") || (type == "trapezium"))
    {
        area_ += yData[0];
        area_ += yData[noPanels];

        for (label i=1; i<noPanels; i++)
        {
            area_ += 2.0*yData[i];
        }

        area_ *= 0.5*binWidth;

        Info << "Integration (trapezium/default) =  " << area_ << endl;
    }

//     else
//     {
//         FatalErrorIn("integral::integral()")
//             << "Cannot find integration type: " << integrationType << nl
//             << " Options are = simpsons, trapezium."
//             << exit(FatalError);
//     }
}




integral::~integral()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const scalar& integral::area() const
{
    return area_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
