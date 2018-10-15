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

   	define(Lx, 6.283)
	define(Lz, 0.01)
   	define(Ly, 0.01)
        
   	define(N_Lx,calc(256))
   	define(N_Lz,calc(1))
   	define(N_Ly,calc(1)) 

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
    	periodic_1
        {
	    type cyclic;
	    neighbourPatch periodic_2;
            faces
            (
		(a0 b0 b3 a3)
	    );
	}

    	bottom
        {
            type empty;
            faces
            (
		(a0 a3 a2 a1)
	    );
	}

    	top
        {
            type empty;
            faces
            (
		(b0 b1 b2 b3)
	    );
	}

    	periodic_2
        {
	    type cyclic;
	    neighbourPatch periodic_1;
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
