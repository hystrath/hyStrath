/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.0.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])

define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

   convertToMeters 1;

	define(lambda_De,0.0021822)

   	define(Lx, calc(20*lambda_De))
	define(Lz, calc(5*lambda_De))
   	define(Ly, calc(lambda_De))
        
   	define(N_Lx,calc(128))
   	define(N_Lz,calc(1))
   	define(N_Ly,calc(128)) 

   	vertices
   	(
    		( 0.0 0.0 0.0) 	vlabel(a0)
		( Lx 0.0 0.0) vlabel(a1)
		( Lx Ly 0.0) vlabel(a2)
		( 0.0 Ly 0.0) vlabel(a3)
 		
    		( 0.0 0.0 -Lz) 	vlabel(b0)
		( Lx 0.0 -Lz) vlabel(b1)
		( Lx Ly -Lz) vlabel(b2)
		( 0.0 Ly -Lz) vlabel(b3)
	);

   	blocks
   	(
				
		hex (a0 a1 b1 b0 a3 a2 b2 b3) (N_Lx N_Ly N_Lz) simpleGrading (1 1 1) //0
	);

   edges
   (
   );

boundary
   (
    	Inlet
        {
            type patch;
            faces
            (
		(a0 b0 b3 a3)
	    );
	}

    	periodic_1
        {
            type cyclic;
	    neighbourPatch periodic_2;
            faces
            (
		(a0 a3 a2 a1)
	    );
	}

    	periodic_2
        {
            type cyclic;
	    neighbourPatch periodic_1;
            faces
            (
		(b0 b1 b2 b3)
	    );
	}

    	Wall
        {
            type wall;
            faces
            (
		(a1 b1 b2 a2)
	    );
	}

    	front
        {
            type empty;
            faces
            (
	
		(a0 a1 b1 b0)
	    );
	}
    	back
        {
            type empty;
            faces
            (
		(a3 b3 b2 a2)
	    );
	}
);

mergePatchPairs
(
);
