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

#include "writeTimeData.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& fileName,
    const List< Pair<scalar> >& data
)
{
    OFstream timeFile(pathName/fileName+".raw");

    if (timeFile.good())
    {
        timeFile << data << endl;
    }
    else
    {
        FatalErrorIn("writeTimeData::writeTimeData()")
            << "Cannot open file " << timeFile.name()
            << abort(FatalError);
    }
}

writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const List< Pair<scalar> >& data,
    const label& dummy
)
{
    fileName writeFile(pathName/nameFile);
    
    scalarField xData(data.size());
    scalarField yData(data.size());

    forAll(data, d)
    {
        xData[d] = data[d].first();
        yData[d] = data[d].second();
    }

    graph outputGraph("title", "x", "y", xData, yData);

    outputGraph.write(writeFile, "raw");
}

//- scalar field
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData
)
{
    OFstream file(pathName/nameFile);

    if(file.good())
    {
        forAll(xData, n)
        {
            file 
                << xData[n]
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

//- vector field
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const vectorField& xData
)
{
    OFstream file(pathName/nameFile);

    if(file.good())
    {
        forAll(xData, n)
        {
            file 
                << xData[n].x() << "\t" 
                << xData[n].y() << "\t" 
                << xData[n].z() << "\t"
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}



//- scalar field, scalar field
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const scalarField& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(xData, n)
            {
                file 
                    << xData[n] << "\t" 
                    << yData[n]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


//- vector field, scalar field, 
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const vectorField& xData,
    const scalarField& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(yData, n)
            {
                file 
                    << xData[n].x() << "\t" 
                    << xData[n].y() << "\t" 
                    << xData[n].z() << "\t"
                    << yData[n]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


//- scalar field, vector field
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(yData, n)
            {
                file 
                    << xData[n] << "\t" 
                    << yData[n].x() << "\t" 
                    << yData[n].y() << "\t"
                    << yData[n].z() 
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


//- scalar field, complex field
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const List<complex>& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(yData, n)
            {
                file 
                    << xData[n] << "\t" 
                    << yData[n].Re() << "\t" 
                    << yData[n].Im() << "\t"
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


//- vector field, vector field, 
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const vectorField& xData,
    const vectorField& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(yData, n)
            {
                file 
                    << xData[n].x() << "\t" 
                    << xData[n].y() << "\t" 
                    << xData[n].z() << "\t"
                    << yData[n].x() << "\t" 
                    << yData[n].y() << "\t"
                    << yData[n].z() 
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


//- scalar field, tensor field, 
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const tensorField& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(yData, n)
            {
                file
                    << xData[n] << "\t"
                    << yData[n].xx() << "\t" << yData[n].xy() << "\t" << yData[n].xz() << "\t"
                    << yData[n].yx() << "\t" << yData[n].yy() << "\t" << yData[n].yz() << "\t"
                    << yData[n].zx() << "\t" << yData[n].zy() << "\t" << yData[n].zz()
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    } 
}

//- vector field, tensor field, 
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const vectorField& xData,
    const tensorField& yData
)
{
    if(xData.size() == yData.size())
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(yData, n)
            {
                file
                    << xData[n].x() << "\t" 
                    << xData[n].y() << "\t" 
                    << xData[n].z() << "\t"
                    << yData[n].xx() << "\t" << yData[n].xy() << "\t" << yData[n].xz() << "\t"
                    << yData[n].yx() << "\t" << yData[n].yy() << "\t" << yData[n].yz() << "\t"
                    << yData[n].zx() << "\t" << yData[n].zy() << "\t" << yData[n].zz()
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    } 
}





// two scalar fields (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const scalarField& yData,
    const bool& dummy
)
{
    if(xData.size() == yData.size())
    {
        fileName fName(pathName/nameFile);
    
        std::ofstream file(fName.c_str(),ios_base::app);
        file.precision(11);

        if(file.is_open())
        {
            forAll(xData, n)
            {
                file
                    << xData[n] << "\t"
                    << yData[n] << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }
    
        file.close();
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


// one scalar field (with append possible) [OLD]
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const bool& dummy
)
{
    if(dummy)
    {
        fileName fName(pathName/nameFile);

        std::ofstream file(fName.c_str(),ios_base::app);
        file.precision(11);

        if(file.is_open())
        {
            forAll(xData, n)
            {
                file << xData[n] << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }

        file.close();
    }
    else
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(xData, n)
            {
                file 
                    << xData[n]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }    
    }
}

// one scalar field - sideways (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const word& option, 
    const bool& dummy     

)
{
    if(option == "once")
    {
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            forAll(xData, n)
            {
                file 
                    << xData[n]
                    << endl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }    
    }
    if(option == "append")
    {
        fileName fName(pathName/nameFile);

        std::ofstream file(fName.c_str(),ios_base::app);
        file.precision(11);

        if(file.is_open())
        {
            forAll(xData, n)
            {
                file << xData[n] << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }

        file.close();
    }    
    if(option == "sidewaysAppend")
    {
//         Pout <<"xData = " << xData << endl;
        
        fileName fName(pathName/nameFile);

        std::ofstream file(fName.c_str(),ios_base::app);
        file.precision(11);

        if(file.is_open())
        {
            forAll(xData, n)
            {
                file << xData[n] << " ";
            }
            
            file << nl;
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }

        file.close();
    }
}


// one scalar field one VECTOR field (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData,
    const bool& dummy
)
{
    if(xData.size() == yData.size())
    {
        fileName fName(pathName/nameFile);
    
        std::ofstream file(fName.c_str(),ios_base::app);
        file.precision(11);
        if(file.is_open())
        {
            forAll(xData, n)
            {
                file
                    << xData[n] << "\t" 
                    << yData[n].x() << "\t" 
                    << yData[n].y() << "\t" 
                    << yData[n].z() << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }
    
        file.close();
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


// one VECTOR field (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const vectorField& yData,
    const bool& dummy
)
{
    fileName fName(pathName/nameFile);

    std::ofstream file(fName.c_str(),ios_base::app);
    file.precision(11);
    if(file.is_open())
    {
        forAll(yData, n)
        {
            file 
                << yData[n].x() << "\t" 
                << yData[n].y() << "\t" 
                << yData[n].z() << nl;
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << fName
            << abort(FatalError);
    }

    file.close();
}




// one scalar field one TENSOR field (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const tensorField& yData,
    const bool& dummy
)
{
    if(xData.size() == yData.size())
    {
        fileName fName(pathName/nameFile);
    
        std::ofstream file(fName.c_str(),ios_base::app);

        file.precision(11);
    
        if(file.is_open())
        {
            forAll(xData, n)
            {
                file
                    << xData[n] << "\t" 
                    << yData[n].xx() << "\t" << yData[n].xy() << "\t" << yData[n].xz() << "\t"
                    << yData[n].yx() << "\t" << yData[n].yy() << "\t" << yData[n].yz() << "\t"
                    << yData[n].zx() << "\t" << yData[n].zy() << "\t" << yData[n].zz()
                    << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }
    
        file.close();
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}

// one scalar field and one List<scalarField> (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const List<scalarField>& yData,
    const bool& dummy
)
{
    if(xData.size() == yData.size())
    {
        fileName fName(pathName/nameFile);
    
        std::ofstream file(fName.c_str(),ios_base::app);
    
        if(file.is_open())
        {
            forAll(xData, n)
            {
                label ySize = yData[n].size();
                
                file 
                    << xData[n] << "\t";

                    
                    for(label i = 0; i < ySize; i++)
                    {
                        file
                            << yData[n][i] << "\t";
                    }
                    
                file
                    << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }
    
        file.close();
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}

// one scalar field and one complex field (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const List<complex>& yData,
    const bool& dummy
)
{
    if(xData.size() == yData.size())
    {
        fileName fName(pathName/nameFile);
    
        std::ofstream file(fName.c_str(),ios_base::app);
    
        if(file.is_open())
        {
            forAll(yData, n)
            {
                file 
                    << xData[n] << "\t" 
                    << yData[n].Re() << "\t" 
                    << yData[n].Im() << nl;
            }
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }
        
        file.close();
    }
    else
    {
        Info << "WARNING: size of two fields for output are not equal: " 
             << xData.size() << " and " << yData.size()
             << nl << " in writeTimeData."
             << endl;
    }
}


// one TENSOR field (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const tensorField& yData,
    const bool& dummy
)
{
    fileName fName(pathName/nameFile);

    std::ofstream file(fName.c_str(),ios_base::app);
    file.precision(11);

    if(file.is_open())
    {
        forAll(yData, n)
        {
            file 
                << yData[n].xx() << "\t" << yData[n].xy() << "\t" << yData[n].xz() << "\t"
                << yData[n].yx() << "\t" << yData[n].yy() << "\t" << yData[n].yz() << "\t"
                << yData[n].zx() << "\t" << yData[n].zy() << "\t" << yData[n].zz()
                << nl;
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << fName
            << abort(FatalError);
    }

    file.close();
}


        
// write out List<scalarField>     component only    
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const List<scalarField>& data
)
{

    OFstream file(pathName/nameFile);

    if(file.good())
    {
        forAll(data, nX)
        {
            forAll(data[nX], nY)
            {
                file << data[nX][nY] << "\t";
            }
            
            file << endl;
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
} 

// write out List<vectorField>        component only
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const List<vectorField>& data,
    const word& option
)
{

    OFstream file(pathName/nameFile);

    if(file.good())
    {
        if(option == "x")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].x() << "\t";
                }
                
                file << endl;
            }
        }
        
        if(option == "y")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].y() << "\t";
                }
                
                file << endl;
            }
        }   
        
        if(option == "z")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].z() << "\t";
                }
                
                file << endl;
            }
        }       
        
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
} 

// List<vectorField> component only (with append possible)
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const List<vectorField>& data,
    const word& option,
    const bool& dummy
)
{
    fileName fName(pathName/nameFile);

    std::ofstream file(fName.c_str(),ios_base::app);
    file.precision(11);

    if(file.is_open())
    {
        forAll(data, i)
        {
            forAll(data[i], j)
            {
            
                if(option == "x")
                {
                    file << data[i][j].x() << " ";
                }
                if(option == "y")
                {
                    file << data[i][j].y() << " ";
                }
                if(option == "z")
                {
                    file << data[i][j].z() << " ";
                }
            }
            
            file << nl;
                
        }
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << fName
            << abort(FatalError);
    }

    file.close();
}

// write out List<tensorField>  component only
writeTimeData::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const List<tensorField>& data,
    const word& option
)
{

    OFstream file(pathName/nameFile);

    if(file.good())
    {
        if(option == "xx")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    file << data[nX][nY].xx() << "\t";
                }
                
                file << endl;
            }
        }
        
        if(option == "xy")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].xy() << "\t";
                }
                
                file << endl;
            }
        }   
        
        if(option == "xz")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].xz() << "\t";
                }
                
                file << endl;
            }
        }
        
        
        if(option == "yx")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    file << data[nX][nY].yx() << "\t";
                }
                
                file << endl;
            }
        }
        
        if(option == "yy")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].yy() << "\t";
                }
                
                file << endl;
            }
        }   
        
        if(option == "yz")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].yz() << "\t";
                }
                
                file << endl;
            }
        }  
        
        
        if(option == "zx")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    file << data[nX][nY].zx() << "\t";
                }
                
                file << endl;
            }
        }
        
        if(option == "zy")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].zy() << "\t";
                }
                
                file << endl;
            }
        }   
        
        if(option == "zz")
        {
            forAll(data, nX)
            {
                forAll(data[nX], nY)
                {
                    
                    file << data[nX][nY].zz() << "\t";
                }
                
                file << endl;
            }
        }        
        
    }
    else
    {
        FatalErrorIn("void writeTimeData::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
} 

writeTimeData::~writeTimeData()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
