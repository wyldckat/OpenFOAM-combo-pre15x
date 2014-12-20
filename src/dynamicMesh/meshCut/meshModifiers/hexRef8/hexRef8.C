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

#include "hexRef8.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "ListOps.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexRef8, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Make sure if one side of coupled patch gets cut the other side gets cut
// as well.
void Foam::hexRef8::syncCoupledCutFaces(labelHashSet& cutFaces) const
{
    if (debug)
    {
        Pout<< "hexRef8::syncCoupledCutFaces : cutFaces before sync:"
            << cutFaces.size() << endl;
    }

    // There can never be a situation that we cannot add to cut faces
    // (a cell can never have too many of them, even if not cut itself) so
    // take the union of the set on both sides of the processor faces.
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Parallel: send all cut patch faces on processor patches
    if (Pstream::parRun())
    {
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (typeid(pp) == typeid(processorPolyPatch))
            {
                // Collect all my cut patch faces
                labelList cutPatchFaces(pp.size());
                label cutI = 0;
                forAll(pp, patchFaceI)
                {
                    if (cutFaces.found(pp.start() + patchFaceI))
                    {
                        cutPatchFaces[cutI++] = patchFaceI;
                    }
                }
                cutPatchFaces.setSize(cutI);

                if (debug)
                {
                    Pout<< "hexRef8::syncCoupledCutFaces : Sending "
                        << cutPatchFaces.size() << " cut faces on "
                        << pp.name() << endl;
                }

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);
                {
                    OPstream toNeighb(procPatch.neighbProcNo());
                    toNeighb << cutPatchFaces;
                }
            }
        }
    }

    // Receive all cut patch faces and merge.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (typeid(pp) == typeid(cyclicPolyPatch))
        {
            label half = pp.size()/2;

            // Half0: merge in data from half1
            for (label patchFaceI = 0; patchFaceI < half; patchFaceI++)
            {
                if (cutFaces.found(pp.start() + patchFaceI))
                {
                    cutFaces.insert(pp.start() + patchFaceI + half);
                }
            }
            // Half1: merge in data from half0
            for (label patchFaceI = half; patchFaceI < pp.size(); patchFaceI++)
            {
                if (cutFaces.found(pp.start() + patchFaceI))
                {
                    cutFaces.insert(pp.start() + patchFaceI - half);
                }
            }
        }
        else if (Pstream::parRun() && typeid(pp) == typeid(processorPolyPatch))
        {
            // Parallel: receive data from neighbour and merge in
            labelList cutPatchFaces;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);
            {
                IPstream fromNeighb(procPatch.neighbProcNo());
                fromNeighb >> cutPatchFaces;
            }

            forAll(cutPatchFaces, i)
            {
                cutFaces.insert(pp.start() + cutPatchFaces[i]);
            }
        }
        else if (pp.coupled())
        {
            label hasCutFaces = false;

            forAll(pp, patchFaceI)
            {
                if (cutFaces.found(pp.start() + patchFaceI))
                {
                    hasCutFaces = true;

                    break;
                }
            }

            if (hasCutFaces)
            {
                FatalErrorIn
                (
                    "hexRef8::syncCoupledCutFaces(labelHashSet& cutFaces) const"
                )   << "Don't know how to handle refining across coupled patch "
                    << pp.name() << " of type " << pp.type()
                    << exit(FatalError);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef8::syncCoupledCutFaces : cutFaces after sync:"
            << cutFaces.size() << endl;
    }
}


// Collect all boundary edges of pp that are in cutEdges. Express edges in
// terms of patchface and starting index in face.
void Foam::hexRef8::collectCutEdges
(
    const labelHashSet& cutEdges,
    const polyPatch& pp,
    DynamicList<label>& cutFaces,
    DynamicList<label>& cutIndex
) const
{
    const edgeList& edges = pp.edges();
    const labelList& meshPoints = pp.meshPoints();

    for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
    {
        const edge& e = edges[edgeI];

        label v0 = meshPoints[e[0]];
        label v1 = meshPoints[e[1]];

        label meshEdgeI = pp.meshEdges()[edgeI];

        if (cutEdges.found(meshEdgeI))
        {
            // Convert edge to patchFace and index in patchFace
            label patchFaceI = pp.edgeFaces()[edgeI][0];

            cutFaces.append(patchFaceI);

            const face& f = pp[patchFaceI];

            // Get index of v0 in face. Might be start or end
            // point of edge so check.
            label fp = findIndex(f, v0);

            if (fp == -1)
            {
                FatalErrorIn("collectCutEdges") << "Problem fp:" << fp
                    << abort(FatalError);
            }


            label fpMin1 = (fp == 0 ? f.size()-1 : fp-1);

            if (f[fpMin1] == v1)
            {
                cutIndex.append(fpMin1);
            }
            else
            {
                cutIndex.append(fp);
            }
        }
    }

    cutFaces.shrink();
    cutIndex.shrink();
}


bool Foam::hexRef8::syncCoupledCutEdges(labelHashSet& cutEdges) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (debug)
    {
        Pout<< "hexRef8::syncCoupledCutEdges : cutEdges before sync:"
            << cutEdges.size() << endl;
    }

    // Parallel: send all cut patch edges on processor patches
    if (Pstream::parRun())
    {
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (typeid(pp) == typeid(processorPolyPatch))
            {
                DynamicList<label> cutFaces(pp.size()/4);
                DynamicList<label> cutIndex(pp.size()/4);

                collectCutEdges(cutEdges, pp, cutFaces, cutIndex);

                if (debug)
                {
                    Pout<< "hexRef8::syncCoupledCutEdges : Sending "
                        << cutFaces.size() << " cut edges on "
                        << pp.name() << endl;
                }

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);
                {
                    OPstream toNeighb(procPatch.neighbProcNo());
                    toNeighb << cutFaces << cutIndex;
                }
            }
        }
    }

    label nOldCutEdges = cutEdges.size();

    // Receive all cut patch edges and merge.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (typeid(pp) == typeid(cyclicPolyPatch))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            DynamicList<label> cutFaces(pp.size()/4);
            DynamicList<label> cutIndex(pp.size()/4);

            collectCutEdges(cutEdges, pp, cutFaces, cutIndex);

            // Transfer patchfaces to other side and convert patchFace+index
            // back into mesh edge.

            forAll(cutFaces, i)
            {
                label patchFaceI = cycPatch.transformLocalFace(cutFaces[i]);

                const face& f = pp[patchFaceI];

                label fp = (f.size() - cutIndex[i]) % f.size();

                // Starting mesh point of edge on this side
                label v0 = f[fp];

                // Next mesh point of edge. Is fp+1 on other side so fp-1 on
                // this side.
                label v1 = f.prevLabel(fp);
                
                // Find edge
                label edgeI = meshTools::findEdge
                (
                    mesh_.edges(),
                    mesh_.faceEdges()[patchFaceI + pp.start()],
                    v0,
                    v1
                );

                if (edgeI == -1)
                {
                    FatalErrorIn("syncCoupledCutEdges") << "Problem"
                        << abort(FatalError);
                }
                cutEdges.insert(edgeI);
            }
        }
        else if (Pstream::parRun() && typeid(pp) == typeid(processorPolyPatch))
        {
            labelList cutFaces;
            labelList cutIndex;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);
            {
                IPstream fromNeighb(procPatch.neighbProcNo());
                fromNeighb >> cutFaces >> cutIndex;
            }

            forAll(cutFaces, i)
            {
                label patchFaceI = cutFaces[i];

                const face& f = pp[patchFaceI];

                label fp = (f.size() - cutIndex[i]) % f.size();

                // Starting mesh point of edge on this side
                label v0 = f[fp];

                // Next mesh point of edge. Is fp+1 on other side so fp-1 on
                // this side.
                label v1 = f.prevLabel(fp);
                
                // Find edge
                label edgeI = meshTools::findEdge
                (
                    mesh_.edges(),
                    mesh_.faceEdges()[patchFaceI + pp.start()],
                    v0,
                    v1
                );

                if (edgeI == -1)
                {
                    FatalErrorIn("syncCoupledCutEdges") << "Problem"
                        << abort(FatalError);
                }
                cutEdges.insert(edgeI);
            }
        }
    }

    if (debug)
    {
        Pout<< "hexRef8::syncCoupledCutEdges : cutEdges after sync:"
            << cutEdges.size() << endl;
    }

    return cutEdges.size() > nOldCutEdges;
}


// Calculates points of cell. Sorts cellPoints in ascending order.
// (Note: only real requirement is that lowest numbered is first but sorting
//  won't harm)
void Foam::hexRef8::calcSortedCellPoints(labelListList& cellPoints) const
{
    cellPoints.setSize(mesh_.nCells());

    labelList nPoints(mesh_.nCells(), 0);

    forAll(mesh_.pointCells(), pointI)
    {
        const labelList& pCells = mesh_.pointCells()[pointI];

        forAll(pCells, pI)
        {
            label cellI = pCells[pI];

            nPoints[cellI]++;
        }
    }

    forAll(cellPoints, cellI)
    {
        cellPoints[cellI].setSize(nPoints[cellI]);
    }

    nPoints = 0;

    forAll(mesh_.pointCells(), pointI)
    {
        const labelList& pCells = mesh_.pointCells()[pointI];

        forAll(pCells, pI)
        {
            label cellI = pCells[pI];

            cellPoints[cellI][nPoints[cellI]++] = pointI;
        }
    }

    // Sort (only so master is first)
    forAll(cellPoints, cellI)
    {
        sort(cellPoints[cellI]);
    }
}


// Get two faces connected to edgeI on cellI. Get them in such a way that
// walking from edge to f0-centre to cell-centre to f1-centre
// produces a normal in the direction of edge[0]
void Foam::hexRef8::getEdgeFaces
(
    const label cellI,
    const label edgeI,
    label& f0,
    label& f1
) const
{
    label fA, fB;
    meshTools::getEdgeFaces(mesh_, cellI, edgeI, fA, fB);

    const edge& e = mesh_.edges()[edgeI];

    label pointI = e[0];

    label baseFaceI = -1;

    const labelList& pFaces = mesh_.pointFaces()[pointI];

    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        if
        (
            faceI != fA
         && faceI != fB
         && meshTools::faceOnCell(mesh_, cellI, faceI)
        )
        {
            baseFaceI = faceI;
            break;
        }
    }

    if (baseFaceI == -1)
    {
        FatalErrorIn("hexRef8::getEdgeFaces")
            << "Cell:" << cellI
            << " edge:" << edgeI
            << " fA:" << fA << " fB:" << fB << abort(FatalError);
    }

    // Now we have
    //   - fA,fB: faces on edgeI
    //   - baseFace: face on e[0]
    // Orient fA,fB such that consistent with outwards pointing normal on
    // pointI

    const face& f = mesh_.faces()[baseFaceI];
    const face& faceA = mesh_.faces()[fA];


    //
    // See if fA is oriented in same way as baseFaceI.
    //

    label baseIndex = findIndex(f, pointI);

    label nextPointI = f[(baseIndex+1) % f.size()];

    // Check if fA uses nextPointI as well.
    label indexA = findIndex(faceA, pointI);

    bool sameOrder = false;

    if
    (
        faceA.nextLabel(indexA) == nextPointI
     || faceA.prevLabel(indexA) == nextPointI
    )
    {
        // fA uses two consecutive points of baseFace so ordering is the
        // same as baseFace.
        sameOrder = true;
    }

    if (mesh_.faceOwner()[baseFaceI] == cellI)
    {
        if (sameOrder)
        {
            f0 = fA;
            f1 = fB;
        }
        else
        {
            f0 = fB;
            f1 = fA;
        }
    }
    else
    {
        if (sameOrder)
        {
            f0 = fB;
            f1 = fA;
        }
        else
        {
            f0 = fA;
            f1 = fB;
        }
    }
}


// Get the two edges using pointI on faceI.
// Such that e0 uses pointI and nextPointI
void Foam::hexRef8::getFaceEdges
(
    const label faceI,
    const label pointI,
    const label nextPointI,
    label& e0,
    label& e1
) const
{
    const labelList& fEdges = mesh_.faceEdges()[faceI];

    e0 = -1;
    e1 = -1;

    forAll(fEdges, i)
    {
        label edgeI = fEdges[i];

        const edge& e = mesh_.edges()[edgeI];

        if (e[0] == pointI)
        {
            if (e[1] == nextPointI)
            {
                e0 = edgeI;
            }
            else
            {
                e1 = edgeI;
            }
        }
        else if (e[1] == pointI)
        {
            if (e[0] == nextPointI)
            {
                e0 = edgeI;
            }
            else
            {
                e1 = edgeI;
            }
        }

        if (e0 != -1 && e1 != -1)
        {
            return;
        }
    }

    FatalErrorIn("hexRef8::getFaceEdges") << "Cannot find two edges on face "
        << faceI << " verts:" << mesh_.faces()[faceI] << " that use point "
        << pointI << abort(FatalError);
}
            

void Foam::hexRef8::getFaceInfo
(
    const label faceI,
    label& patchID,
    label& zoneID,
    label& zoneFlip
) const
{
    patchID = -1;

    if (!mesh_.isInternalFace(faceI))
    {
        patchID = mesh_.boundaryMesh().whichPatch(faceI);
    }

    zoneID = mesh_.faceZones().whichZone(faceI);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
    }
}


// Adds a face on top of existing faceI.
void Foam::hexRef8::addFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if ((nei == -1) || (own < nei))
    {
        // Ordering ok.
        meshMod.setAction
        (
            polyAddFace
            (
                newFace,                    // face
                own,                        // owner
                nei,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
    else
    {
        // Reverse owner/neighbour
        meshMod.setAction
        (
            polyAddFace
            (
                newFace.reverseFace(),      // face
                nei,                        // owner
                own,                        // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceI,                      // master face for addition
                false,                      // flux flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face zone flip
            )
        );
    }
}


// Adds an internal face from an edge. Assumes orientation correct.
// Problem is that the face is between four new vertices. So what do we provide
// as master? The only existing mesh item we have is the edge we have split.
// Have to be careful in only using it if it has internal faces since otherwise 
// polyMeshMorph will complain (because it cannot generate a sensible mapping
// for the face)
void Foam::hexRef8::addInternalFace
(
    polyTopoChange& meshMod,
    const label meshEdgeI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    // Check if edge has any internal faces we can use.
    label masterEdgeI = -1;

    const labelList& eFaces = mesh_.edgeFaces()[meshEdgeI];

    forAll(eFaces, i)
    {
        if (mesh_.isInternalFace(eFaces[i]))
        {
            // meshEdge uses internal faces so ok to inflate from it
            masterEdgeI = meshEdgeI;

            break;
        }
    }
    

    meshMod.setAction
    (
        polyAddFace
        (
            newFace,                    // face
            own,                        // owner
            nei,                        // neighbour
            -1,                         // master point
            masterEdgeI,                // master edge
            -1,                         // master face for addition
            false,                      // flux flip
            -1,                         // patch for face
            -1,                         // zone for face
            false                       // face zone flip
        )
    );
}


// Modifies existing faceI for either new owner/neighbour or new face points.
void Foam::hexRef8::modFace
(
    polyTopoChange& meshMod,
    const label faceI,
    const face& newFace,
    const label own,
    const label nei
) const
{
    label patchID, zoneID, zoneFlip;

    getFaceInfo(faceI, patchID, zoneID, zoneFlip);

    if
    (
        (own != mesh_.faceOwner()[faceI])
     || (
            mesh_.isInternalFace(faceI)
         && (nei != mesh_.faceNeighbour()[faceI])
        )
     || (newFace != mesh_.faces()[faceI])
    )
    {
        if ((nei == -1) || (own < nei))
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,            // modified face
                    faceI,              // label of face being modified
                    own,                // owner
                    nei,                // neighbour
                    false,              // face flip
                    patchID,            // patch for face
                    false,              // remove from zone
                    zoneID,             // zone for face
                    zoneFlip            // face flip in zone
                )
            );
        }
        else
        {
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    faceI,                  // label of face being modified
                    nei,                    // owner
                    own,                    // neighbour
                    false,                  // face flip
                    patchID,                // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }
}


// Given cell and point on cell determine new neighbour.
Foam::label Foam::hexRef8::newCell
(
    const labelListList& cellPoints,
    const label cellI,
    const label pointI
) const
{
    if (addedCells_[cellI].size() == 0)
    {
        // Unsplit cell.
        return cellI;
    }

    // Get local vertex numbering of pointI in cellPoints (master first)
    label index = findIndex(cellPoints[cellI], pointI);

    if (index == -1 || index >= addedCells_[cellI].size())
    {
        FatalErrorIn("hexRef8::newCell")
            << "Cannot find point " << pointI << " on cell " << cellI << nl
            << "Cell has vertices " << cellPoints[cellI] << abort(FatalError);
    }

    // See what new cell is on that cell/point combination.
    return addedCells_[cellI][index];
}


// Get new owner and neighbour (in unspecified order) of pointI on faceI.
void Foam::hexRef8::getFaceNeighbours
(
    const labelListList& cellPoints,
    const label faceI,
    const label pointI,

    label& own,
    label& nei
) const
{
    label oldOwn = mesh_.faceOwner()[faceI];

    own = newCell(cellPoints, oldOwn, pointI);

    nei = -1;

    if (mesh_.isInternalFace(faceI))
    {
        label oldNei = mesh_.faceNeighbour()[faceI];

        nei = newCell(cellPoints, oldNei, pointI);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::hexRef8::hexRef8(const polyMesh& mesh)
:
    mesh_(mesh),
    addedCells_(mesh.nCells())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hexRef8::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        Pout<< "hexRef8:" << endl;
    }

    //
    // Simple check on the mesh
    //

    const cellModel& hex = *(cellModeller::lookup("hex"));

    const cellShapeList& cellShapes = mesh_.cellShapes();

    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        const cellModel& model = cellShapes[cellI].model();
        
        if (model != hex)
        {
            FatalErrorIn
            (
                "hexRef8::setRefinement(const labelList&, polyTopoChange&)"
            )   << "cell to refine:" << cellI
                << " is not a hex" << nl
                << "cellShape:" << cellShapes[cellI] << abort(FatalError);
        }
    }


    // Determine whether mesh has coupled patches and needs to sync

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    bool hasCouples = false;

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            hasCouples = true;

            break;
        }
    }


    //
    // Get faces that will be cut
    //

    labelHashSet cutFaces(mesh_.nFaces()/10 + 10);

    forAll(cellLabels, i)
    {
        const labelList& cFaces = mesh_.cells()[cellLabels[i]];

        forAll(cFaces, i)
        {
            cutFaces.insert(cFaces[i]);
        }
    }

    // Synchronize cut faces on both sides of coupled patches
    if (hasCouples)
    {
        syncCoupledCutFaces(cutFaces);
    }


    //
    // Get edges that will be cut
    //

    labelHashSet cutEdges(mesh_.nEdges()/10 + 10);

    for
    (
        labelHashSet::const_iterator iter = cutFaces.begin();
        iter != cutFaces.end();
        ++iter
    )
    {
        const labelList& fEdges = mesh_.faceEdges()[iter.key()];

        forAll(fEdges, i)
        {
            cutEdges.insert(fEdges[i]);
        }
    }

    // Synchronize cut edges on both sides of coupled patches
    if (hasCouples)
    {
        while (true)
        {
            bool changed = syncCoupledCutEdges(cutEdges);

            reduce(changed, orOp<bool>());

            if (!changed)
            {
                break;
            }
        }
    }


    //
    // Get faces that will be affected but not cut
    //

    labelHashSet affectedFaces(mesh_.nFaces()/10 + 10);

    for
    (
        labelHashSet::const_iterator iter = cutEdges.begin();
        iter != cutEdges.end();
        ++iter
    )
    {
        const labelList& eFaces = mesh_.edgeFaces()[iter.key()];

        forAll(eFaces, i)
        {
            label faceI = eFaces[i];

            if (!cutFaces.found(faceI))
            {
                affectedFaces.insert(faceI);
            }
        }
    }

    if (debug)
    {
        Pout<< "    Cells cut     : " << cellLabels.size() << endl
            << "    Faces cut     : " << cutFaces.size() << endl
            << "    Edges cut     : " << cutEdges.size() << endl
            << "    Faces affected: " << affectedFaces.size() << endl;
    }

    // Sorted list of vertices per cell. Sorted so master point (lowest number)
    // is first. Determine also for unrefined cells
    labelListList cellPoints(mesh_.nCells());
    calcSortedCellPoints(cellPoints);


    // Added points: cell centre
    labelList cellPoint(mesh_.nCells(), -1);
    // Added points: face centre
    labelList facePoint(mesh_.nFaces(), -1);
    // Added points: edge mid
    labelList edgePoint(mesh_.nEdges(), -1);

    // Add cell centre points
    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        label masterPointI = cellPoints[cellI][0];

        cellPoint[cellI] =
            meshMod.setAction
            (
                polyAddPoint
                (
                    mesh_.cellCentres()[cellI],  // point
                    masterPointI,               // master point
                    -1,                         // zone for point
                    true                        // supports a cell
                )
            );
    }

    if (debug)
    {
        Pout<< "    Points added: " << cellLabels.size() << " cell centres"
            << endl;
    }

    // Add face centre points
    for
    (
        labelHashSet::const_iterator iter = cutFaces.begin();
        iter != cutFaces.end();
        ++iter
    )
    {
        label faceI = iter.key();

        label masterPointI = min(mesh_.faces()[faceI]);

        facePoint[faceI] =
            meshMod.setAction
            (
                polyAddPoint
                (
                    mesh_.faceCentres()[faceI],  // point
                    masterPointI,               // master point
                    -1,                         // zone for point
                    true                        // supports a cell
                )
            );
    }

    if (debug)
    {
        Pout<< "    Points added: " << cutFaces.size() << " face centres"
            << endl;
    }


    // Add edge-mid points
    for
    (
        labelHashSet::const_iterator iter = cutEdges.begin();
        iter != cutEdges.end();
        ++iter
    )
    {
        label edgeI = iter.key();

        const edge& e = mesh_.edges()[edgeI];

        label masterPointI = min(e.start(), e.end());

        const point& v0 = mesh_.points()[e.start()];
        const point& v1 = mesh_.points()[e.end()];

        point newPt = 0.5*(v0 + v1);

        edgePoint[edgeI] =
            meshMod.setAction
            (
                polyAddPoint
                (
                    newPt,              // point
                    masterPointI,       // master point
                    -1,                 // zone for point
                    true                // supports a cell
                )
            );
    }

    if (debug)
    {
        Pout<< "    Points added: " << cutEdges.size() << " edge mids"
            << endl;
    }


    // Add cells and store in addedCells.
    label nAdded = 0;
    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        label nPoints = cellPoints[cellI].size();

        labelList& myCells = addedCells_[cellI];

        myCells.setSize(nPoints);

        // Master
        myCells[0] = cellI;

        for (label i = 1; i < nPoints; i++)
        {
            myCells[i] =
                meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,                 // master point
                        -1,                 // master edge
                        -1,                 // master face
                        cellI,              // master cell
                        -1                  // zone for cell
                    )
                );
            nAdded++;
        }
    }

    if (debug)
    {
        Pout<< "    Cells added: " << nAdded << endl;
    }


    // Split original faces into 1 original and 3 new ones.
    nAdded = 0;
    for
    (
        labelHashSet::const_iterator iter = cutFaces.begin();
        iter != cutFaces.end();
        ++iter
    )
    {
        label faceI = iter.key();

        const face& f = mesh_.faces()[faceI];

        label masterPointI = min(f);

        forAll(f, fp)
        {
            label pointI = f[fp];
            label nextI = f[(fp + 1) % f.size()];

            // Get edges on face using pointI. (note: in face order:
            // e0=between pointI,nextI)
            label e0, e1;
            getFaceEdges(faceI, pointI, nextI, e0, e1);

            // Construct new face (in same order as original face!)
            face newFace(4);

            newFace[0] = pointI;
            newFace[1] = edgePoint[e0];
            newFace[2] = facePoint[faceI];
            newFace[3] = edgePoint[e1];

            label own, nei;
            getFaceNeighbours
            (
                cellPoints,
                faceI,
                pointI,

                own,
                nei
            );

            if (pointI == masterPointI)
            {
                modFace(meshMod, faceI, newFace, own, nei);
            }
            else
            {
                addFace(meshMod, faceI, newFace, own, nei);
                nAdded++;
            }
        }
    }

    if (debug)
    {
        Pout<< "    Faces split: " << cutFaces.size()
            << "   modified:" << cutFaces.size()
            << " added:" << nAdded << endl;
    }


    // Add all internal faces. These go from cell-centre to cut edge.
    nAdded = 0;
    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        const labelList& cEdges = mesh_.cellEdges()[cellI];

        forAll(cEdges, i)
        {
            // Get faces using edgeI in order consistent with e[0]
            label edgeI = cEdges[i];
            label f0, f1;
            getEdgeFaces(cellI, edgeI, f0, f1);

            const edge& e = mesh_.edges()[edgeI];

            // Owner and neighbour determined by sides of edge.
            label c0 = newCell(cellPoints, cellI, e[0]);
            label c1 = newCell(cellPoints, cellI, e[1]);

            if (c0 < c1)
            {
                // c0 is master.

                // Face ordered to point to e[0]
                face newF(4);
                newF[0] = cellPoint[cellI];
                newF[1] = facePoint[f0];
                newF[2] = edgePoint[edgeI];
                newF[3] = facePoint[f1];

                addInternalFace
                (
                    meshMod,
                    edgeI,
                    newF,
                    c0,
                    c1
                );
            }
            else
            {
                // c1 is master.

                // Face ordered to point to e[1]
                face newF(4);
                newF[0] = cellPoint[cellI];
                newF[1] = facePoint[f1];
                newF[2] = edgePoint[edgeI];
                newF[3] = facePoint[f0];

                addInternalFace
                (
                    meshMod,
                    edgeI,
                    newF,
                    c1,
                    c0
                );
            }
            nAdded++;
        }
    }

    if (debug)
    {
        Pout<< "    Faces added: " << nAdded
            << " faces internal to original cells" << endl;
    }


    // Faces that only get one or more edges split.
    for
    (
        labelHashSet::const_iterator iter = affectedFaces.begin();
        iter != affectedFaces.end();
        ++iter
    )
    {
        label faceI = iter.key();

        const face& f = mesh_.faces()[faceI];

        DynamicList<label> verts(2*f.size());

        forAll(f, fp)
        {
            label v0 = f[fp];

            verts.append(f[fp]);

            label v1 = f.nextLabel(fp);

            label edgeI = meshTools::findEdge(mesh_, v0, v1);

            if (edgePoint[edgeI] != -1)
            {
                verts.append(edgePoint[edgeI]);
            }
        }
        face newF(verts.shrink());

        label own = mesh_.faceOwner()[faceI];

        label nei = -1;

        if (mesh_.isInternalFace(faceI))
        {
            nei = mesh_.faceNeighbour()[faceI];
        }

        modFace(meshMod, faceI, newF, own, nei);
    }

    if (debug)
    {
        Pout<< "    Faces modified: " << affectedFaces.size() << " for"
            << " split edges" << endl;
    }
}


// ************************************************************************* //
