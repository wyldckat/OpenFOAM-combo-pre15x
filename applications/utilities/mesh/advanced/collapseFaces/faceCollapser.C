/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

#include "faceCollapser.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "ListOps.H"
#include "polyAddFace.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Return point as factor of edge between start and end.
// 0: pt == start
// 1: pt == end
Foam::scalar Foam::faceCollapser::edgeWeight
(
    const point& start,
    const point& end,
    const point& pt
)
{
    return mag(pt - start)/mag(end - start);
}


// Insert labelList into labelHashSet. Optional excluded element.
void Foam::faceCollapser::insert
(
    const labelList& elems,
    const label excludeElem,
    labelHashSet& set
)
{
    forAll(elems, i)
    {
        if (elems[i] != excludeElem)
        {
            set.insert(elems[i]);
        }
    }
}


// Find edge amongst candidate edges. FatalError if none.
Foam::label Foam::faceCollapser::findEdge
(
    const edgeList& edges,
    const labelList& edgeLabels,
    const label v0,
    const label v1
)
{
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];

        const edge& e = edges[edgeI];
    
        if
        (
            (e[0] == v0 && e[1] == v1)
         || (e[0] == v1 && e[1] == v0)
        )
        {
            return edgeI;
        }
    }

    FatalErrorIn("findEdge") << "Cannot find edge between vertices " << v0
        << " and " << v1 << " in edge labels " << edgeLabels
        << abort(FatalError);

    return -1;
}


// Gets list of new vertices + position and sorts them.
void Foam::faceCollapser::insertSorted
(
    const DynamicList<scalarLabel>& edgeCuts,
    DynamicList<label>& newFace
)
{
    scalar startWeight = -GREAT;

    forAll(edgeCuts, i)
    {
        scalar minWeight = GREAT;
        label minPointI = -1;

        forAll(edgeCuts, i)
        {
            scalar weight = edgeCuts[i].val();

            if (weight > startWeight && weight < minWeight)
            {
                minWeight = weight;
                minPointI = edgeCuts[i].index();
            }
        }

        // Use minPointI at minWeightI
        newFace.append(minPointI);

        startWeight = minWeight;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Find edge in edgeSet that is nearest to pointI. Return edgeI and nearest
// point on edge.
void Foam::faceCollapser::findNearestEdgePoint
(
    const labelHashSet& edgeSet,
    const label faceI,
    const label pointI,
    label& minEdgeI,
    point& minEdgePt
) const
{
    const point& pt = mesh_.points()[pointI];

    const labelList& fEdges = mesh_.faceEdges()[faceI];

    minEdgeI = -1;
    scalar minDist = GREAT;

    forAll(fEdges, i)
    {
        label edgeI = fEdges[i];

        if (edgeSet.found(edgeI))
        {
            const edge& e = mesh_.edges()[edgeI];

            point edgePt(e.line(mesh_.points()).nearestDist(pt).rawPoint());

            scalar dist = magSqr(edgePt - pt);

            if (dist < minDist)
            {
                minDist = dist;
                minEdgeI = edgeI;
                minEdgePt = edgePt;
            }
        }
    }
}


// Replace vertices in face
void Foam::faceCollapser::filterFace
(
    const Map<DynamicList<scalarLabel> >& splitEdges,
    const Map<label>& pointMap,
    const label faceI,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[faceI];
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    // Space for replaced vertices and split edges.
    DynamicList<label> newFace(10 * f.size());

    forAll(f, fp)
    {
        label v0 = f[fp];

        Map<label>::const_iterator iter = pointMap.find(v0);

        if (iter != pointMap.end())
        {
            newFace.append(iter());
        }
        else
        {
            newFace.append(v0);
        }

        // Look ahead one to get edge.
        label fp1 = (fp + 1) % f.size();

        label v1 = f[fp1];

        // Get split on edge if any.
        label edgeI = findEdge(mesh_.edges(), fEdges, v0, v1);

        Map<DynamicList<scalarLabel> >::const_iterator edgeFnd =
            splitEdges.find(edgeI);

        if (edgeFnd != splitEdges.end())
        {
            // edgeI has been split (by introducing new vertices).
            // Insert new vertices in face.
            insertSorted(edgeFnd(), newFace);
        }
    }
    face newF(newFace.shrink());

    //Info<< "Modifying face:" << faceI << " from " << f << " to " << newFace
    //    << endl;

    label nei = -1;

    label patchI = -1;

    if (mesh_.isInternalFace(faceI))
    {
        nei = mesh_.faceNeighbour()[faceI];
    }
    else
    {
        patchI = mesh_.boundaryMesh().whichPatch(faceI);
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            newF,                       // modified face
            faceI,                      // label of face being modified
            mesh_.faceOwner()[faceI],   // owner
            nei,                        // neighbour
            false,                      // face flip
            patchI,                     // patch for face
            false,                      // remove from zone
            -1,                         // zone for face
            false                       // face flip in zone
        )
    );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::faceCollapser::faceCollapser(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceCollapser::setRefinement
(
    const labelList& faceLabels,
    const labelList& fpStart,
    const labelList& fpEnd,
    polyTopoChange& meshMod
) const
{
    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();
    const faceList& faces = mesh_.faces();
    const labelListList& faceEdges = mesh_.faceEdges();
    const labelListList& edgeFaces = mesh_.edgeFaces();
    const labelListList& pointFaces = mesh_.pointFaces();


    labelHashSet keepPoints(4*faceLabels.size());
    labelHashSet keepEdges(4*faceLabels.size());

    forAll(faceLabels, i)
    {
        const label faceI = faceLabels[i];

        const face& f = faces[faceI];

        const label fpA = fpStart[i];
        const label fpB = fpEnd[i];

        // Mark all vertices from fpA..fpB as being kept.
        label fp = fpA;

        while (true)
        {
            keepPoints.insert(f[fp]);

            if (fp == fpB)
            {
                break;
            }

            fp = (fp + 1) % f.size();
        }

        // Check all other vertices (fpB..fpA) that they are not marked as being
        // removed (by another face). Would be conflict and is callers
        // responsability.
        do
        {
            fp = (fp + 1) % f.size();

            if (fp == fpA)
            {
                break;
            }

            if (keepPoints.found(f[fp]))
            {
                FatalErrorIn
                (
                    "faceCollapser::setRefinement(const labelList&"
                    ", const labelList&, const labelList&, polyTopoChange&"
                )   << "Face:" << faceI << " vertices:" << f
                    << " has conflict in that one of its vertices to be"
                    << " removed is specified by another face to be kept"
                    << endl
                    << "Vertices to be kept between indices " << fpA
                    << " and " << fpB << endl
                    << "Conflict at vertex to be removed " << f[fp]
                    << " at index " << fp << abort(FatalError);
            }
        } while(true);


        // Determine edges beteen fpA and fpB to keep
        const labelList& fEdges = faceEdges[faceI];

        forAll(fEdges, i)
        {
            const edge& e = edges[fEdges[i]];

            if (keepPoints.found(e[0]) && keepPoints.found(e[1]))
            {
                keepEdges.insert(fEdges[i]);
            }
        }
    }

    keepPoints.clear();


    // Determine projected non-keep points. Can either project onto new
    // vertex on edge or existing vertex.

    // From removed point to kept point. Point to be mapped to can be
    // either existing point or newly introduced point on (now split) edge.
    Map<label> removedPoints(faceLabels.size());

    // From split edge to newly introduced point(s). Can be more than one per
    // edge!
    Map<DynamicList<scalarLabel> > splitEdges(faceLabels.size());

    // Mark faces in any way affect by modifying points/edges. Used later on
    // to prevent having to redo all faces.
    labelHashSet affectedFaces(4*faceLabels.size());


    //
    // Add/remove vertices and construct mapping
    //

    forAll(faceLabels, i)
    {
        const label faceI = faceLabels[i];

        const face& f = faces[faceI];

        const label fpA = fpStart[i];
        const label fpB = fpEnd[i];


        label fp = (fpB + 1) % f.size();

        while (fp != fpA)
        {
            label v = f[fp];

            // Remove point
            meshMod.setAction(polyRemovePoint(v));

            // Mark all faces affected
            insert(pointFaces[v], faceI, affectedFaces);

            // Find nearest point on a keepEdge as:
            //    - edgeI
            //    - point on edge
            label edgeI;
            point edgePt;
            findNearestEdgePoint(keepEdges, faceI, v, edgeI, edgePt);

            // Snap to edge endpoint if too close
            const edge& e = edges[edgeI];

            scalar weight = edgeWeight(points[e[0]], points[e[1]], edgePt);

            if (weight < 0.05)
            {
                // Replace v by e[0]
                removedPoints.insert(v, e[0]);

                Info<< "Replacing point "
                    << v << " at " << points[v]
                    << " with existing vertex " << e[0] << " at " 
                    << points[e[0]] << endl;
            }
            else if (weight > 0.95)
            {
                // Replace v by e[1]
                removedPoints.insert(v, e[1]);

                Info<< "Replacing point "
                    << v << " at " << points[v]
                    << " with existing vertex " << e[1] << " at "
                    << points[e[1]] << endl;
            }
            else
            {
                // - Create point
                // - Split edgeI
                // - Replace v by new point
                label addedVertI = 
                    meshMod.setAction
                    (
                        polyAddPoint
                        (
                            edgePt,     // point
                            v,          // master point
                            -1,         // zone for point
                            true        // supports a cell
                        )
                    );

                Info<< "Added new point " << addedVertI << " on edge " << e
                    << " at " << edgePt
                    << " to replace " << v << " at " << points[v] << endl;

                removedPoints.insert(v, addedVertI);

                // Mark all faces affected
                insert(edgeFaces[edgeI], faceI, affectedFaces);

                // Store weighting factor and new vertex label on the edge.
                scalarLabel edgeCut(weight, addedVertI);

                Map<DynamicList<scalarLabel> >::iterator edgeFnd =
                    splitEdges.find(edgeI);

                if (edgeFnd != splitEdges.end())
                {
                    edgeFnd().append(edgeCut);
                }
                else
                {
                    DynamicList<scalarLabel> edgeData;
                    edgeData.append(edgeCut);

                    splitEdges.insert(edgeI, edgeData);
                }
            }

            fp = (fp + 1) % f.size();
        }
    }


    // Shrink splitEdges.
    for
    (
        Map<DynamicList<scalarLabel> >::iterator iter =
            splitEdges.begin();
        iter != splitEdges.end();
        ++iter
    )
    {
        iter().shrink();
    }


    //
    // Remove faces.
    //

    forAll(faceLabels, i)
    {
        const label faceI = faceLabels[i];

        //Info<< "Removing face:" << faceI << endl;
        meshMod.setAction(polyRemoveFace(faceI));

        // Update list of faces we still have to modify
        affectedFaces.erase(faceI);
    }

    //
    // Modify faces affected
    //

    for
    (
        labelHashSet::const_iterator iter = affectedFaces.begin();
        iter != affectedFaces.end();
        ++iter
    )
    {
        filterFace(splitEdges, removedPoints, iter.key(), meshMod);
    }
}


// ************************************************************************* //
