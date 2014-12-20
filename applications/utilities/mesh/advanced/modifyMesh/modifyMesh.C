/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Manipulates mesh boundary elements.

    Actions are:
        (boundary)points:
            - move

        (boundary)edges:
            - split and move introduced point

        (boundary)faces:
            - split(triangulate) and move introduced point

    For later:
    - collapse edges
    - split cells into pyramids around point
    - split cells with plane

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "morphMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "boundaryCutter.H"
#include "meshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label findPoint(const primitivePatch& pp, const point& nearPoint)
{
    const pointField& points = pp.points();
    const labelList& meshPoints = pp.meshPoints();

    // Find nearest and next nearest
    scalar minDistSqr = GREAT;
    label minI = -1;

    scalar almostMinDistSqr = GREAT;
    label almostMinI = -1;

    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        scalar distSqr = magSqr(nearPoint - points[pointI]);

        if (distSqr < minDistSqr)
        {
            almostMinDistSqr = minDistSqr;
            almostMinI = minI;

            minDistSqr = distSqr;
            minI = pointI;
        }
        else if (distSqr < almostMinDistSqr)
        {
            almostMinDistSqr = distSqr;
            almostMinI = pointI;
        }
    }


    // Decide if nearPoint unique enough.
    Info<< "Found to point " << nearPoint << nl
        << "    nearest point      : " << minI
        << " distance " <<  Foam::sqrt(minDistSqr)
        << " at " << points[minI] << nl
        << "    next nearest point : " << almostMinI
        << " distance " <<  Foam::sqrt(almostMinDistSqr)
        << " at " << points[almostMinI] << endl;

    if (almostMinDistSqr < 4*minDistSqr)
    {
        Info<< "Next nearest too close to nearest. Aborting" << endl;

        return -1;
    }
    else
    {
        return minI;
    }
}


label findEdge
(
    const primitiveMesh& mesh,
    const primitivePatch& pp,
    const point& nearPoint
)
{
    const pointField& localPoints = pp.localPoints();
    const pointField& points = pp.points();
    const labelList& meshPoints = pp.meshPoints();
    const edgeList& edges = pp.edges();

    // Find nearest and next nearest
    scalar minDist = GREAT;
    label minI = -1;

    scalar almostMinDist = GREAT;
    label almostMinI = -1;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        pointHit pHit(e.line(localPoints).nearestDist(nearPoint));

        if (pHit.hit())
        {
            if (pHit.distance() < minDist)
            {
                almostMinDist = minDist;
                almostMinI = minI;

                minDist = pHit.distance();
                minI = meshTools::findEdge
                (
                    mesh,
                    meshPoints[e[0]],
                    meshPoints[e[1]]
                );
            }
            else if (pHit.distance() < almostMinDist)
            {
                almostMinDist = pHit.distance();
                almostMinI = meshTools::findEdge
                (
                    mesh,
                    meshPoints[e[0]],
                    meshPoints[e[1]]
                );
            }
        }
    }

    if (minI == -1)
    {
        Info<< "Did not find edge close to point " << nearPoint << endl;

        return -1;
    }


    // Decide if nearPoint unique enough.
    Info<< "Found to point " << nearPoint << nl
        << "    nearest edge      : " << minI
        << " distance " << minDist << " endpoints "
        << mesh.edges()[minI].line(points) << nl
        << "    next nearest edge : " << almostMinI
        << " distance " << almostMinDist << " endpoints "
        << mesh.edges()[almostMinI].line(points) << nl
        << endl;

    if (almostMinDist < 2*minDist)
    {
        Info<< "Next nearest too close to nearest. Aborting" << endl;

        return -1;
    }
    else
    {
        return minI;
    }
}


label findFace
(
    const primitiveMesh& mesh,
    const primitivePatch& pp,
    const point& nearPoint
)
{
    const pointField& points = pp.points();

    // Find nearest and next nearest
    scalar minDist = GREAT;
    label minI = -1;

    scalar almostMinDist = GREAT;
    label almostMinI = -1;

    forAll(pp, patchFaceI)
    {
        pointHit pHit(pp[patchFaceI].nearestPoint(nearPoint, points));

        if (pHit.hit())
        {
            if (pHit.distance() < minDist)
            {
                almostMinDist = minDist;
                almostMinI = minI;

                minDist = pHit.distance();
                minI = patchFaceI + mesh.nInternalFaces();
            }
            else if (pHit.distance() < almostMinDist)
            {
                almostMinDist = pHit.distance();
                almostMinI = patchFaceI + mesh.nInternalFaces();
            }
        }
    }

    if (minI == -1)
    {
        Info<< "Did not find face close to point " << nearPoint << endl;

        return -1;
    }


    // Decide if nearPoint unique enough.
    Info<< "Found to point " << nearPoint << nl
        << "    nearest face      : " << minI
        << " distance " << minDist
        << " with face centre " << mesh.faceCentres()[minI] << nl
        << "    next nearest face : " << almostMinI
        << " distance " << almostMinDist
        << " with face centre " << mesh.faceCentres()[almostMinI] << nl
        << endl;

    if (almostMinDist < 2*minDist)
    {
        Info<< "Next nearest too close to nearest. Aborting" << endl;

        return -1;
    }
    else
    {
        return minI;
    }
}



// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Reading mesh for time = " << runTime.value() << endl;

    Info<< "Create mesh\n" << endl;
    morphMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    Info<< "Reading modifyMeshDict\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "modifyMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<FixedList<point, 2> > pointsToMove(dict.lookup("pointsToMove"));
    List<FixedList<point, 2> > edgesToSplit(dict.lookup("edgesToSplit"));
    List<FixedList<point, 2> > facesToTriangulate
    (
        dict.lookup("facesToTriangulate")
    );


    // Get calculating engine for all of outside
    const SubList<face> outsideFaces
    (
        mesh.faces(),
        mesh.nFaces() - mesh.nInternalFaces(),
        mesh.nInternalFaces()
    );

    primitivePatch allBoundary(outsideFaces, mesh.points());


    // Look up mesh labels and convert to input for boundaryCutter.
    Map<point> pointToPos(pointsToMove.size());
    Map<List<point> > edgeToCuts(edgesToSplit.size());
    Map<point> faceToSplit(facesToTriangulate.size());


    bool validInputs = true;

    Info<< nl << "Looking up points to move ..." << nl << endl;
    forAll(pointsToMove, i)
    {
        const FixedList<point, 2>& pts = pointsToMove[i];

        label pointI = findPoint(allBoundary, pts[0]);

        if (pointI == -1 || !pointToPos.insert(pointI, pts[1]))
        {
            Info<< "Could not insert mesh point " << pointI
                << " for input point " << pts[0] << nl
                << "Perhaps the point is already marked for moving?" << endl;
            validInputs = false;
        }
    }


    Info<< nl << "Looking up edges to split ..." << nl << endl;
    forAll(edgesToSplit, i)
    {
        const FixedList<point, 2>& pts = edgesToSplit[i];

        label edgeI = findEdge(mesh, allBoundary, pts[0]);

        if (edgeI == -1 || !edgeToCuts.insert(edgeI, List<point>(1, pts[1])))
        {
            Info<< "Could not insert mesh edge " << edgeI
                << " for input point " << pts[0] << nl
                << "Perhaps the edge is already marked for cutting?" << endl;

            validInputs = false;
        }
    }


    Info<< nl << "Looking up faces to triangulate ..." << nl << endl;
    forAll(facesToTriangulate, i)
    {
        const FixedList<point, 2>& pts = facesToTriangulate[i];

        label faceI = findFace(mesh, allBoundary, pts[0]);

        if (faceI == -1 || !faceToSplit.insert(faceI, pts[1]))
        {
            Info<< "Could not insert mesh face " << faceI
                << " for input point " << pts[0] << nl
                << "Perhaps the face is already marked for splitting?" << endl;

            validInputs = false;
        }
    }

    if (!validInputs)
    {
        Info<< nl << "There was a problem in one of the inputs in the"
            << " dictionary. Not modifying mesh." << endl;
    }
    else
    {
        Info<< nl << "All input points located. Modifying mesh." << endl;


        // Mesh change engine
        boundaryCutter cutter(mesh);

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Insert commands into meshMod
        cutter.setRefinement
        (
            pointToPos,
            edgeToCuts,
            faceToSplit,
            meshMod
        );

        mesh.setMorphTimeIndex(runTime.timeIndex());
        mesh.updateTopology(meshMod);
        cutter.updateTopology(mesh.morphMap());
        mesh.movePoints(mesh.morphMap().preMotionPoints());

        runTime++;

        // Write resulting mesh
        Info << "Writing modified mesh to time " << runTime.value() << endl;
        mesh.write();
    }


    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
