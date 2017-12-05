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

#include "polyARCHER.H"
#include "graph.H"
#include "OFstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

void polyARCHER::outputInitialisation()
{
    
    {
        // deletes current content of file
        OFstream file(pathName_/nameOfFile_);
    
        if(file.good())
        {
            file << endl;
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }

    // Initialisation
    
    fileName fName(pathName_/nameOfFile_);

    std::ofstream file(fName.c_str(),ios_base::app);
    file.precision(11);

    if(file.is_open())
    {
        file << "*--------------------------------*- C++ -*----------------------------------*" << nl;
        file << "| =========                |                                                 | " << nl;
        file << "| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | " << nl;
        file << "|  \\    /   O peration     | Version:  2.3.x                                 | " << nl;
        file << "|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | " << nl;
        file << "|    \\/     M anipulation  |                                                 | " << nl;
        file << "*---------------------------------------------------------------------------*" << nl;
        
        file << "OpenFOAM running on ARCHER - temporary log file" <<  nl;
        file << "number of procs = " << Pstream::nProcs() << nl;
        file << nl;
    }
    else
    {
        FatalErrorIn("void polyARCHER()")
            << "Cannot open file " << fName
            << abort(FatalError);
    }

    file.close();
}

void polyARCHER::outputTime()
{
    fileName fName(pathName_/nameOfFile_);

    scalar TE = getTotalEnergy();
    
    std::ofstream file(fName.c_str(),ios_base::app);
    file.precision(11);

    if(file.is_open())
    {
        file<< "Time = " << time_.timeOutputValue() << nl;
        file<< "total energy = " << TE << nl;
        file<< "Duration: " << molCloud_.clock().instantDuration()
         << " s   av. write int. = " << molCloud_.clock().averageTimeWriteInterval()
         << " s   av. sim. = " << molCloud_.clock().averageTime()
         << " s   tot. = " << molCloud_.clock().totalDuration() << " s"
         << nl;
        
        file<< "ExecutionTime = " << time_.elapsedCpuTime() << " s"
            << "  ClockTime = " << time_.elapsedClockTime() << " s"
            << nl;
        file<< nl;            
    }
    else
    {
        FatalErrorIn("void polyARCHER()")
            << "Cannot open file " << fName
            << abort(FatalError);
    }

    file.close();    
    
}

} // End namespace Foam

// ************************************************************************* //
