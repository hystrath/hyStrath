/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2005 OpenCFD Ltd.
    \\/      M anipulation   |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    makeAxialMesh

Description
    Creates a axial-symmetric grid from a two-dimensional grid

    Needs two inputs from the command line:
    - name of the boundary that is the symmetry axis
    - name of the boundary that is to be splitted into two wedge-boundaries (front and back)
    - an option that moves the real axis away (offset)

    If no options are specified on the command line, it looks for a dictionary
    rotationDict. In this dictionary the rotational axis is specified. In addition the 
    patch that is on that axis may be specified.

    Assumes the following about the grid:
    - 2D-Grid (just one cell thick)
    - parallel to the XY-plane
    - the boundary that is the symmetry axis is a straight line
    - the "front" and "back" are one patch

    Afterwards: use the collapseEdges utility to remove the faces at the symmetry axis
 
\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "Time.H"
#include "repatchPolyTopoChanger.H"

#include "plane.H"
#include "line.H"
#include "faceSet.H"
#include "SortableList.H"
#include "SubField.H"

#include "mathematicalConstants.H"

using namespace Foam;

typedef line<point,point> linie;

// get the plane of the grid
// right now it assumes parallel to XY

plane getPlane(pointField& points) {
  scalar x=0,y=0,z=0;
  label  n=0;

  forAll(points,pointI) {
    point &pt=points[pointI];
    x += pt.x();    y += pt.y();    z += pt.z();
    n ++;
  }

  plane result(point(x/n,y/n,z/n),vector(0,0,1));

  return result;
}

// calculate the symmetry axis from the "symmetry patch"
// by finding the two points that are farthest apart

linie getAxis(const polyPatch &axis,plane &pl) {
  scalar dist=0;
  int cnt=0;
  point pt1(0,0,0),pt2(0,0,0);
  const pointField &pts=axis.localPoints();

  forAll(pts,pointI) {
    switch(cnt) {
    case 0:
      pt1=pl.nearestPoint(pts[pointI]);
      break;
    case 1:
      pt2=pl.nearestPoint(pts[pointI]);
      dist=mag(pt1-pt2);
      break;
    default:
      point pt=pl.nearestPoint(pts[pointI]);
      if( mag(pt-pt1)>dist || mag(pt-pt2)>dist) {
	if( mag(pt-pt1) > mag(pt-pt2)) {
	  pt2=pt;
	} else {
	  pt1=pt;
	}
	dist=mag(pt1-pt2);
      }
    }
    cnt++;
  }

  return linie(pt1,pt2);
}

scalar getDistance(const polyPatch &axisPatch,linie &axisLine) {
  const pointField &pts=axisPatch.localPoints();
  vector axisDir = axisLine.vec()/axisLine.mag();      

  scalar distance=0;

  forAll(pts,pointI) {
    const point &pt=pts[pointI];

    point axisPoint = axisLine.start() + axisDir * (axisDir & (pt-axisLine.start()));
    scalar radius=mag( axisPoint - pt );
    if(radius>distance) {
        distance=radius;
    }
  }
  
  return distance;
}

// angle (as suggested in the OpenFOAM-Userguide)
const scalar defaultAngle=5;   //DPS the plane is rotated +2.5 degrees and -2.5 degrees to make a 5 deg. wedge

// changes the coordinates of the points

/*
    This docu provided by E David Huckaby
    (plus a fix to a bug)

    Inputs:
        cutPlane - INFINITE plane (getPlane)
            - defined by position vector and normal vector
        axisLine - FINITE line (getAxis)
            - defined by two points
            - axis of rotation
        mesh - 
    

    Variables:
        "basePoint" is the projection of the cartesian point, "oldPoint"
            onto the cutting plane "cutPlane". This is equivalent to
            the nearest point from "oldPoint" to "cutPlane".
        "axisPoint" is the projection of "basePoint" along the axis, "axisLine".
            If "axisLine" was an INFINITE line then this would be equilivalent to the
            the neareast point from "axisLine" to "basePoint".  This is not the case,
            so the projection is done with a dot product ("&") with the unit vector in the
            direction of "axisLine", "axisDir".    
        "radius" is the radial coordinate of the cartesian point, in the plane
            perpendicular to the rotation axis.
	"radiusBasedPoint" is the basePoint coordinates relative to the axisPoint.
	    It is used if revolve = true to intruduce the factorRadius in order to
	    revolve the mesh instead of projecting it on wedges

*/

void changeCoordinates(
    polyMesh &mesh,
    plane cutPlane,
    linie &axisLine,
    scalar offset,
    const scalar wedgeAngle,
    bool revolve) 
{
  if (revolve)    Info << "Revolving nodes" << endl;
  else 		  Info << "Projecting nodes" << endl;
  const scalar angle=wedgeAngle/2;

  repatchPolyTopoChanger topo(mesh);

  pointField oldPoints=mesh.points();                //DPS the old points will be rotated too
  pointField newPoints(oldPoints.size());

  const scalar factor=std::sin(angle/180.*constant::mathematical::pi);
  const scalar factorRadius=std::cos(angle/180.*constant::mathematical::pi);
  vector axisDir = axisLine.vec()/axisLine.mag();      

  scalar minRadius=1e10,maxRadius=-1e10;

  scalar dist1=mag(axisLine.start()-cutPlane.nearestPoint(axisLine.start()));
  scalar dist2=mag(axisLine.end()-cutPlane.nearestPoint(axisLine.end()));

  if(dist1>SMALL || dist2>SMALL) {
      Warning << " End points of axis " << axisLine << " are " << dist1 
          << " and " << dist2 << " away from plane "  << cutPlane << endl;
  }

  forAll(oldPoints,pointI) {
      //    point oldPoint=oldPoints[pointI];
    const point &oldPoint=oldPoints[pointI];
    point basePoint=cutPlane.nearestPoint(oldPoint);

    scalar radius=GREAT;

    point axisPoint = axisLine.start() + axisDir * (axisDir & (basePoint-axisLine.start()));
    radius=mag( axisPoint - basePoint )+offset;
    point radiusBasedPoint = factorRadius * (basePoint - axisPoint);

    if(radius>maxRadius) {
        maxRadius=radius;
    }
    if(radius<minRadius) {
        minRadius=radius;
    }
    vector dir(oldPoint-basePoint);
    if(revolve) {
        newPoints[pointI] = axisPoint + radiusBasedPoint + factor * radius * dir/mag(dir);
    }
    else newPoints[pointI] = basePoint + factor * radius * dir/mag(dir);
    //    oldPoints[pointI] = basePoint - factor *radius * dir/mag(dir);  //DPS rotating original plane
  }

  Info << "Radius to axis: min = " << minRadius << " max = " << maxRadius << endl;

  //  mesh.movePoints(oldPoints);
  mesh.movePoints(newPoints);
}

// Split the wedge-patch into two patches
 
void splitWedge(polyMesh &mesh,word wname,plane pl) {
  repatchPolyTopoChanger topo(mesh);

  const polyBoundaryMesh& patches = mesh.boundaryMesh();
  const polyPatch &wedge=patches[patches.findPatchID(wname)];
  
  const vectorField::subField 	& fcs=wedge.faceCentres ();

  faceSet facesPos(mesh,"set1_"+wname,fcs.size()/2,IOobject::NO_WRITE);
  faceSet facesNeg(mesh,"set2_"+wname,fcs.size()/2,IOobject::NO_WRITE);

  // assign faces to the two patches
  forAll(fcs,faceI) {
    vector center=fcs[faceI];

    scalar dir=pl.normal() & (center-pl.nearestPoint(center));
    if( dir >0 ) {
      // This offset is dirty, there must be a better way
      facesPos.insert(faceI-wedge.whichFace(0));
    } else {
      facesNeg.insert(faceI-wedge.whichFace(0));
    }
  }

  Info << " Copying patches " << endl;

  // the rest is stolen from createPatch (others might say inspired by)

  List<polyPatch*> newPatches(patches.size() + 2);
  
  forAll(patches, patchI)
    {
      const polyPatch& pp = patches[patchI];
      
      newPatches[patchI] =
	pp.clone
	(
	 patches,
	 patchI,
	 pp.size(),
	 pp.start()
	 ).ptr();
      
    }
  
  label patchPos = newPatches.size() - 2;
  label patchNeg = newPatches.size() - 1;
  
  Info << " Creating Patches" << endl;

  newPatches[patchPos] =
    polyPatch::New
    (
     "empty",
     wname+"_pos",
     0,
     mesh.nFaces(),
     patchPos,
     patches
     ).ptr();

   newPatches[patchNeg] =
     polyPatch::New
     (
      "empty",
      wname+"_neg",
      0,
      mesh.nFaces(),
      patchNeg,
      patches
      ).ptr();
  

  // Actually add new list of patches
  topo.changePatches(newPatches);
  
  Info << " Creating Pos-patch " << endl;

  labelList faceLabelsPos(facesPos.toc());

  SortableList<label> patchFacesPos(faceLabelsPos);

  forAll(patchFacesPos, i)
    {
      label faceI = patchFacesPos[i];

      if (mesh.isInternalFace(faceI))
	{
	  FatalErrorIn("SplitWedge")
	    << "Face " << faceI << " specified in set " << facesPos.name()
	    << " is not an external face of the mesh." << endl
	    << "This application can only repatch existing boundary"
	    << " faces."
	    << exit(FatalError);
	}

      topo.changePatchID(patchFacesPos[i], patchPos);
    }

  Info << " Creating Neg-patch " << endl;

  labelList faceLabelsNeg(facesNeg.toc());

  SortableList<label> patchFacesNeg(faceLabelsNeg);

  forAll(patchFacesNeg, i)
    {
      label faceI = patchFacesNeg[i];

      if (mesh.isInternalFace(faceI))
	{
	  FatalErrorIn("SplitWedge")
	    << "Face " << faceI << " specified in set " << facesNeg.name()
	    << " is not an external face of the mesh." << endl
	    << "This application can only repatch existing boundary"
	    << " faces."
	    << exit(FatalError);
	}
      
      topo.changePatchID(patchFacesNeg[i], patchNeg);
    }

  Info << " Changing patches\n" << endl;

  topo.repatch();
}

void changeTypes(polyMesh &mesh,word wedge,word axis,bool hasOffset) {
  repatchPolyTopoChanger topo(mesh);
  const polyBoundaryMesh& patches = mesh.boundaryMesh();

  List<polyPatch*> newPatches(patches.size());
  
  forAll(patches, patchI)
    {
      const polyPatch& pp = patches[patchI];

      if(pp.name()==axis && !hasOffset) {
          Info << " Changing " << axis << " to symmetryPlane" << endl;
          newPatches[patchI] =
              polyPatch::New(
                  "symmetryPlane",
                  pp.name(),
                  pp.size(),
                  pp.start(),
                  patchI,
                  patches
              ).ptr();
      } else if(pp.name()==wedge+"_pos" || pp.name()==wedge+"_neg"){
          Info << " Changing " << pp.name() << " to wedge" << endl;
          newPatches[patchI] =
              polyPatch::New(
                  "wedge",
                  pp.name(),
                  pp.size(),
                  pp.start(),
                  patchI,
                  patches
              ).ptr();
      } else {
          newPatches[patchI] =
              pp.clone
              (
                  patches,
                  patchI,
                  pp.size(),
                  pp.start()
              ).ptr();
      }
    }

  Info << endl;

  topo.changePatches(newPatches);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    argList::validOptions.insert("axis","<axis face name>");
    argList::validOptions.insert("wedge","<wedge face  name>");
    argList::validOptions.insert("offset","<additional offset from axis>");
    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("wedgeAngle","<degrees>");
    argList::validOptions.insert("revolve","");
  
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Create mesh for time = " << runTime.value() << nl << endl;

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    const word oldInstance = mesh.pointsInstance();

    word axisName;
    word wedgeName;

    bool oldMode=false;
    bool revolve = args.options().found("revolve");

    scalar offset=0;
    scalar wedgeAngle=defaultAngle;

    vector rotation(0,0,0),origin(0,0,0);

    bool overwrite = args.options().found("overwrite");

    if(args.options().found("axis") && args.options().found("wedge")) {
      oldMode=true;
      Info << " Using old mode. Dictionary not used" << endl;
      axisName=args.options()["axis"];
      wedgeName=args.options()["wedge"];

      if(args.options().found("offset")) {
	offset=readScalar(IStringStream(args.options()["offset"])());
      }
    } else if(args.options().found("axis") || args.options().found("wedge")) {
        FatalErrorIn(args.executable())
            << " axis and wedge options have to be specified for legacy mode "
            << exit(FatalError);
    } else {
        IOdictionary rotationalDict
	  (
	   IOobject
	   (
	    "rotationDict",
	    runTime.system(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	    )
	   );
        
	revolve = readBool(rotationalDict.lookup("revolve"));

	if(rotationalDict.found("makeAxialOldMode") &&
	   readBool(rotationalDict.lookup("makeAxialOldMode"))) {
	  oldMode=true;
	  Info << "Using old mode" << endl;

	  if(rotationalDict.found("makeAxialOffset")) {
	    offset=readScalar(rotationalDict["makeAxialOffset"]);
	  }
	} else {
	  rotation=rotationalDict.lookup("rotationVector");
	  origin=rotationalDict.lookup("originVector");
	}
        axisName=word(rotationalDict.lookup("makeAxialAxisPatch"));

	wedgeName=word(rotationalDict["makeAxialWedgePatch"]);

        if( axisName == wedgeName) {
            FatalErrorIn(args.executable())
                << "Patch " << axisName << " can't be used as axis and wedge patch" << endl
                    << exit(FatalError);	
        }
    
        if(rotationalDict.found("wedgeAngle")) {
	    wedgeAngle=readScalar(rotationalDict["wedgeAngle"]);
        }
    }

    if(args.options().found("wedgeAngle")) {
	wedgeAngle=readScalar(IStringStream(args.options()["wedgeAngle"])());
    }

    if(offset<0) {
      FatalErrorIn(args.executable())
	<< "Offset " << offset << " smaller than 0" << endl
	<< exit(FatalError);	
    }
    if(offset>0) {
	Info << " Adding offset " << offset << endl;
    }

    pointField& points = const_cast<pointField&>(mesh.points());

    plane approx=getPlane(points);
    linie theAxis(point(0,0,0),point(0,0,0));

    Info << "Plane of the grid: " << approx << "\n" << endl;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (patches.findPatchID(wedgeName) == -1)
      {
	FatalErrorIn(args.executable())
	  << "Patch " << wedgeName << " does not exists in mesh" << endl
	  << "Patches are " << patches.names()
	  << exit(FatalError);
      }
    
    if(oldMode) {
      if (patches.findPatchID(axisName) == -1)
	{
	  FatalErrorIn(args.executable())
            << "Patch " << axisName << " does not exists in mesh" << endl
            << "Patches are " << patches.names()
            << exit(FatalError);
	}

      const polyPatch &axisPatch=patches[patches.findPatchID(axisName)];
      
      theAxis=getAxis(axisPatch,approx);
    } else {
      theAxis=linie(origin,origin+rotation);
      offset=0;
    }
    
    Info << "The rotation-axis: " << theAxis << "\n" << endl;

    Info << "Creating wedge with an opening angle of " << wedgeAngle << " degrees\n" << endl;

    changeCoordinates(
        mesh,
        approx,
        theAxis,
        offset,
        wedgeAngle,
	revolve
    );

    scalar distance=offset;

    if(!oldMode) {
        const polyPatch &axisPatch=patches[patches.findPatchID(axisName)];

        distance=getDistance(axisPatch,theAxis);

        Info << " Calculated Distance of " << axisName << " to axis is " << distance << endl;
    }

    if(distance>SMALL) {
        Info << "Distance of axis to patch is " << distance 
            << " -> not changing the patch type" << endl;
    }

    Info << "Splitting patch " << wedgeName << "\n" << endl;

    splitWedge(mesh,wedgeName,approx);

    Info << "Changing patch types " << endl;

    changeTypes(mesh,wedgeName,axisName,distance>SMALL);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    Info << "Writing mesh to time " << runTime.value() << endl;
    mesh.write();

    Info << "End\n" << endl;

    Info << "Now use collapseEdges to clean grid\n";

    return 0;
}


// ************************************************************************* //
