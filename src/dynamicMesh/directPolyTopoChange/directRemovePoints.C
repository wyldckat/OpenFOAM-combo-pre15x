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

\*---------------------------------------------------------------------------*/

#include "directRemovePoints.H"
#include "PstreamReduceOps.H"
#include "polyMesh.H"
#include "directPolyTopoChange.H"
#include "polyRemovePoint.H"
#include "polyModifyFace.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(directRemovePoints, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Change the vertices of the face whilst keeping everything else the same.
void Foam::directRemovePoints::modifyFace
(
    const label faceI,
    const face& newFace,
    directPolyTopoChange& meshMod
) const
{
    // Get other face data.
    label patchI = -1;
    label owner = mesh_.faceOwner()[faceI];
    label neighbour = -1;

    if (mesh_.isInternalFace(faceI))
    {
        neighbour = mesh_.faceNeighbour()[faceI];
    }
    else
    {
        patchI = mesh_.boundaryMesh().whichPatch(faceI);

        // Check that we're not modifying coupled faces.
        if (mesh_.boundaryMesh()[patchI].coupled())
        {
            FatalErrorIn("mergeEdges") << "Problem for face:" << faceI
                << abort(FatalError);
        }
    }

    label zoneID = mesh_.faceZones().whichZone(faceI);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            newFace,        // modified face
            faceI,          // label of face being modified
            owner,          // owner
            neighbour,      // neighbour
            false,          // face flip
            patchI,         // patch for face
            false,          // remove from zone
            zoneID,         // zone for face
            zoneFlip        // face flip in zone
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::directRemovePoints::directRemovePoints(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::directRemovePoints::countPointUsage
(
    const scalar minCos,
    boolList& pointCanBeDeleted
) const
{
    // Containers to store two edges per point:
    // -1   : not filled
    // >= 0 : edge label
    // -2   : more than two edges for point
    labelList edge0(mesh_.nPoints(), -1);
    labelList edge1(mesh_.nPoints(), -1);

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        forAll(e, eI)
        {
            label pointI = e[eI];

            if (edge0[pointI] == -2)
            {
                // Already too many edges
            }
            else if (edge0[pointI] == -1)
            {
                // Store first edge using point
                edge0[pointI] = edgeI;
            }
            else
            {
                // Already one edge using point. Check second container.

                if (edge1[pointI] == -1)
                {
                    // Store second edge using point
                    edge1[pointI] = edgeI;
                }
                else
                {
                    // Third edge using point. Mark.
                    edge0[pointI] = -2;
                    edge1[pointI] = -2;
                }
            }
        }
    }


    // Check the ones used by only 2 edges that these are sufficiently in line.
    const pointField& points = mesh_.points();

    pointCanBeDeleted.setSize(mesh_.nPoints());
    pointCanBeDeleted = false;
    label nDeleted = 0;

    forAll(edge0, pointI)
    {
        if (edge0[pointI] >= 0 && edge1[pointI] >= 0)
        {
            // Point used by two edges exactly

            const edge& e0 = edges[edge0[pointI]];
            const edge& e1 = edges[edge1[pointI]];

            label common = e0.commonVertex(e1);
            label vLeft = e0.otherVertex(common);
            label vRight = e1.otherVertex(common);

            vector e0Vec = points[common] - points[vLeft];
            e0Vec /= mag(e0Vec) + VSMALL;

            vector e1Vec = points[vRight] - points[common];
            e1Vec /= mag(e1Vec) + VSMALL;

            if ((e0Vec & e1Vec) > minCos)
            {
                pointCanBeDeleted[pointI] = true;
                nDeleted++;
            }
        }
        else if (edge0[pointI] == -1)
        {
            // point not used at all
            pointCanBeDeleted[pointI] = true;
            nDeleted++;
        }

    }
    edge0.clear();
    edge1.clear();


    // Point can be deleted only if all processors want to delete it
    syncTools::syncPointList
    (
        mesh_,
        pointCanBeDeleted,
        andEqOp<bool>(),
        true,               // null value
        false               // no separation
    );

    return returnReduce(nDeleted, sumOp<label>());
}


void Foam::directRemovePoints::setRefinement
(
    const boolList& pointCanBeDeleted,
    directPolyTopoChange& meshMod
)
{
    // Remove points
    // ~~~~~~~~~~~~~

    // Faces (in mesh face labels) affected by points removed. Will hopefully
    // be only a few.
    labelHashSet facesAffected;

    forAll(pointCanBeDeleted, pointI)
    {
        if (pointCanBeDeleted[pointI])
        {
            meshMod.setAction(polyRemovePoint(pointI));

            // Store faces affected
            const labelList& pFaces = mesh_.pointFaces()[pointI];

            forAll(pFaces, i)
            {
                facesAffected.insert(pFaces[i]);
            }
        }
    }


    // Update faces
    // ~~~~~~~~~~~~

    forAllConstIter(labelHashSet, facesAffected, iter)
    {
        label faceI = iter.key();

        const face& f = mesh_.faces()[faceI];

        face newFace(f.size());

        label newI = 0;

        forAll(f, fp)
        {
            label pointI = f[fp];

            if (!pointCanBeDeleted[pointI])
            {
                newFace[newI++] = pointI;
            }
        }
        newFace.setSize(newI);

        modifyFace(faceI, newFace, meshMod);
    }
}


// ************************************************************************* //
