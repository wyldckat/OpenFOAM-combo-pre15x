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

#include "removeFaces.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "ListOps.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(removeFaces, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::removeFaces::nexti(const label size, const label index)
{
    return (index + 1) % size;
}


Foam::label Foam::removeFaces::previ(const label size, const label index)
{
    return (index == 0) ? size - 1: index - 1;
}


bool Foam::removeFaces::sameFaceOrdering
(
    const edge& e,
    const face& f0,
    const label f0Start,    // index of edge start on f0
    const face& f1,
    const label f1Start     // index of edge start on f1
)
{
    if (f0.nextLabel(f0Start) == e.end())
    {
        if (f1.nextLabel(f1Start) != e.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if (f1.nextLabel(f1Start) == e.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}    


// Get angle between two faces shared by edge.
Foam::scalar Foam::removeFaces::faceCos
(
    const label face0I,
    const label face1I,
    const bool checkOrder,
    const label sharedEdgeI
) const
{
    vector n0 = mesh_.faceAreas()[face0I];
    n0 /= mag(n0) + VSMALL;

    vector n1 = mesh_.faceAreas()[face1I];
    n1 /= mag(n1) + VSMALL;

    scalar cosAngle = n0 & n1;

    if (checkOrder)
    {
        const face& f0 = mesh_.faces()[face0I];
        const face& f1 = mesh_.faces()[face1I];

        const edge& e = mesh_.edges()[sharedEdgeI];

        if
        (
           !sameFaceOrdering
            (
                e,
                f0,
                findIndex(f0, e.start()),
                f1,
                findIndex(f1, e.start())
            )
        )
        {
            cosAngle = -cosAngle;
        }
    }

    return cosAngle;
}


// Get patch, zone info for faceI
void Foam::removeFaces::getFaceInfo
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


// Count number of remaining faces using edgeI. Returns true if only two
// remaining and sets face0I, face1I to those two.
bool Foam::removeFaces::getUnremovedFaces
(
    const labelHashSet& facesToRemove,
    const label edgeI,

    label& face0I,
    label& face1I
) const
{
    const labelList& eFaces = mesh_.edgeFaces()[edgeI];

    face0I = -1;
    face1I = -1;

    forAll(eFaces, eFaceI)
    {
        label faceI = eFaces[eFaceI];

        if (!facesToRemove.found(faceI))
        {
            if (face0I == -1)
            {
                face0I = faceI;
            }
            else if (face1I == -1)
            {
                face1I = faceI;
            }
            else
            {
                // Already two unremoved faces stored, this is third one.
                return false;
            }
        }
    }

    if (face1I != -1)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Count number of remaining edges using vertI. Returns true if only two
// remaining and sets edge0I, edge1I to those two.
bool Foam::removeFaces::getUnremovedEdges
(
    const Map<edge>& edgesToRemove,
    const label vertI,

    label& edge0I,
    label& edge1I
) const
{
    const labelList& pEdges = mesh_.pointEdges()[vertI];

    edge0I = -1;
    edge1I = -1;

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if (!edgesToRemove.found(edgeI))
        {
            if (edge0I == -1)
            {
                edge0I = edgeI;
            }
            else if (edge1I == -1)
            {
                edge1I = edgeI;
            }
            else
            {
                // Already two unremoved edges stored, this is third one.
                return false;
            }
        }
    }

    if (edge1I != -1)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Checks vertI for if it can be removed. (so if only two edges using it
// and these two are in line).
bool Foam::removeFaces::canBeRemoved
(
    const Map<edge>& edgesToRemove,
    const label vertI
) const
{
    label edge0I = -1;
    label edge1I = -1;

    if (getUnremovedEdges(edgesToRemove, vertI, edge0I, edge1I))
    {
        // Have two remaining edges: edge0I, edge1I. Can be merged if
        // in same direction.

        const edge& e0 = mesh_.edges()[edge0I];

        vector e0Vec = e0.vec(mesh_.points());
        e0Vec /= mag(e0Vec) + VSMALL;

        const edge& e1 = mesh_.edges()[edge1I];

        vector e1Vec = e1.vec(mesh_.points());
        e1Vec /= mag(e1Vec) + VSMALL;

        scalar cosAngle = e0Vec & e1Vec;

        if ((e0.start() == e1.start()) || (e0.end() == e1.end()))
        {
            cosAngle = -cosAngle;
        }

        if (cosAngle > minCos_)
        {
            return true;
        }
    }
    return false;
}


// Return face with all pointsToRemove removed.
Foam::face Foam::removeFaces::filterFace
(
    const labelHashSet& pointsToRemove,
    const label faceI
) const
{
    const face& f = mesh_.faces()[faceI];

    labelList newFace(f.size(), -1);

    label newFp = 0;

    forAll(f, fp)
    {
        label vertI = f[fp];

        if (!pointsToRemove.found(vertI))
        {
            newFace[newFp++] = vertI;
        }
    }

    newFace.setSize(newFp);

    if (debug && (newFp != f.size()))
    {
        Pout<< "Filtered face " << faceI << " old verts:" << f
            << " new verts:" << newFace << endl;
    }

    return face(newFace);
}


// Combines two faces, leaving out points on the edge they share and all points
// in pointsToRemove. Ordering of points is as face0.
Foam::face Foam::removeFaces::mergeFace
(
    const labelHashSet& pointsToRemove,
    const label face0I,
    const label face1I,
    const label sharedEdgeI
) const
{
    const edge& e = mesh_.edges()[sharedEdgeI];

    const face& f0 = mesh_.faces()[face0I];

    const face& f1 = mesh_.faces()[face1I];

    labelList newFace(f0.size() + f1.size() - 2);

    label newFp = 0;

    label face0Start = findIndex(f0, e.start());
    label next0 = nexti(f0.size(), face0Start);

    bool edgeOrderOK;

    if (f0[next0] == e.end())
    {
        edgeOrderOK = true;

        // Walk from next0 to face0Start
        label fp = next0;

        do
        {
            label vertI = f0[fp];

            if (!pointsToRemove.found(vertI))
            {
                newFace[newFp++] = vertI;
            }
            fp = nexti(f0.size(), fp);
        }
        while (fp != next0);
    }
    else
    {
        edgeOrderOK = false;

        // Walk from face0Start to next0
        label fp = face0Start;

        do
        {
            label vertI = f0[fp];

            if (!pointsToRemove.found(vertI))
            {
                newFace[newFp++] = vertI;
            }
            fp = nexti(f0.size(), fp);
        }
        while (fp != face0Start);
    }

    // Step to face1. Determine whether same topological direction.
    label face1Start = findIndex(f1, e.start());
    label next1 = nexti(f1.size(), face1Start);


    // No need to swap face1 order if face1
    // uses point of edgeI in different order compared to face0 so
    // if face0 uses edge from start-end and face1 from end-start (or vice
    // versa)

    bool face1OrderOk = (edgeOrderOK ^ (f1[next1] == e.end()));

    if (face1OrderOk)
    {
        label fp = face1Start;

        do
        {
            label vertI = f1[fp];

            if
            (
                (vertI != e.end())         // Lazy, could work this out above.
             && (vertI != e.start())
             && !pointsToRemove.found(vertI)
            )
            {
                newFace[newFp++] = vertI;
            }

            fp = nexti(f1.size(), fp);
        }
        while (fp != face1Start);
    }
    else
    {
        // Reverse face1

        label fp = face1Start;

        do
        {
            label vertI = f1[fp];

            if
            (
                (vertI != e.end())         // Lazy, could work this out above.
             && (vertI != e.start())
             && !pointsToRemove.found(vertI)
            )
            {
                newFace[newFp++] = vertI;
            }
            fp = previ(f1.size(), fp);
        }
        while (fp != face1Start);
    }
    newFace.setSize(newFp);

    return face(newFace);
}
        

// Wrapper for meshMod.modifyFace. Reverses face if own>nei.
void Foam::removeFaces::modFace
(
    const face& f,
    const label masterFaceID,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label newPatchID,
    const bool removeFromZone,
    const label zoneID,
    const bool zoneFlip,

    polyTopoChange& meshMod
) const
{
    if ((nei == -1) || (own < nei))
    {
        if (debug)
        {
            Pout<< "ModifyFace (unreversed) :"
                << "  faceI:" << masterFaceID
                << "  f:" << f
                << "  own:" << own
                << "  nei:" << nei
                << "  flipFaceFlux:" << flipFaceFlux
                << "  newPatchID:" << newPatchID
                << "  removeFromZone:" << removeFromZone
                << "  zoneID:" << zoneID
                << "  zoneFlip:" << zoneFlip
                << endl;
        }

        meshMod.setAction
        (
            polyModifyFace
            (
                f,              // modified face
                masterFaceID,   // label of face being modified
                own,            // owner
                nei,            // neighbour
                flipFaceFlux,   // face flip
                newPatchID,     // patch for face
                removeFromZone, // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }
    else
    {
        if (debug)
        {
            Pout<< "ModifyFace (!reversed) :"
                << "  faceI:" << masterFaceID
                << "  f:" << f.reverseFace()
                << "  own:" << nei
                << "  nei:" << own
                << "  flipFaceFlux:" << flipFaceFlux
                << "  newPatchID:" << newPatchID
                << "  removeFromZone:" << removeFromZone
                << "  zoneID:" << zoneID
                << "  zoneFlip:" << zoneFlip
                << endl;
        }

        meshMod.setAction
        (
            polyModifyFace
            (
                f.reverseFace(),// modified face
                masterFaceID,   // label of face being modified
                nei,            // owner
                own,            // neighbour
                flipFaceFlux,   // face flip
                newPatchID,     // patch for face
                removeFromZone, // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }
}


// Changes region of connected set of cells
void Foam::removeFaces::changeCellRegion
(
    const label cellI,
    const label oldRegion,
    const label newRegion,
    labelList& cellRegion
) const
{
    if (cellRegion[cellI] == oldRegion)
    {
        cellRegion[cellI] = newRegion;

        // Step to neighbouring cells

        const labelList& cCells = mesh_.cellCells()[cellI];

        forAll(cCells, i)
        {
            changeCellRegion(cCells[i], oldRegion, newRegion, cellRegion);
        }
    }
}


// Mark all faces affected in any way by
// - removal of cells
// - removal of faces
// - removal of edges
// - removal of points
Foam::labelHashSet Foam::removeFaces::getFacesAffected
(
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    const labelHashSet& facesToRemove,
    const Map<edge>& edgesToRemove,
    const labelHashSet& pointsToRemove
) const
{
    labelHashSet facesAffected
    (
        cellRegionMaster.size()
      + facesToRemove.size()
      + edgesToRemove.size()
      + pointsToRemove.size()
    );

    // Mark faces affected by removal of cells
    forAll(cellRegion, cellI)
    {
        label region = cellRegion[cellI];

        if (region != -1 && (cellI != cellRegionMaster[region]))
        {
            const labelList& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                facesAffected.insert(cFaces[cFaceI]);
            }
        }
    }


    // Mark faces affected by removal of face. (note: could have used
    // faceLabels)
    for
    (
        labelHashSet::const_iterator iter = facesToRemove.begin();
        iter != facesToRemove.end();
        ++iter
    )
    {
         facesAffected.insert(iter.key());
    }


    //  Mark faces affected by removal of edges
    for
    (
        Map<edge>::const_iterator iter = edgesToRemove.begin();
        iter != edgesToRemove.end();
        ++iter
    )
    {
        label edgeI = iter.key();

        const labelList& eFaces = mesh_.edgeFaces()[edgeI];

        forAll(eFaces, eFaceI)
        {
            facesAffected.insert(eFaces[eFaceI]);
        }
    }

    // Mark faces affected by removal of points
    for
    (
        labelHashSet::const_iterator iter = pointsToRemove.begin();
        iter != pointsToRemove.end();
        ++iter
    )
    {
        label pointI = iter.key();

        const labelList& pFaces = mesh_.pointFaces()[pointI];

        forAll(pFaces, pFaceI)
        {
            facesAffected.insert(pFaces[pFaceI]);
        }
    }
    return facesAffected;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::removeFaces::removeFaces(const polyMesh& mesh, const scalar minCos)
:
    mesh_(mesh),
    minCos_(minCos)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Removing face connects cells. This function works out a consistent set of
// cell regions. 
// - returns faces to remove. Can be extended with additional faces
//   (if owner would become neighbour)
// - sets cellRegion to -1 or to region number
// - regionMaster contains for every region number a master cell.
Foam::label Foam::removeFaces::compatibleRemoves
(
    const labelList& facesToRemove,
    labelList& cellRegion,
    labelList& regionMaster,
    labelList& newFacesToRemove
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    cellRegion.setSize(mesh_.nCells());
    cellRegion = -1;

    regionMaster.setSize(mesh_.nCells());
    regionMaster = -1;

    label nRegions = 0;

    forAll(facesToRemove, i)
    {
        label faceI = facesToRemove[i];

        if (!mesh_.isInternalFace(faceI))
        {
            FatalErrorIn
            (
                "removeFaces::compatibleRemoves(const labelList&"
                ", labelList&, labelList&, labelList&)"
            )   << "Not internal face:" << faceI << abort(FatalError);
        }


        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        if (own >= nei)
        {
            FatalErrorIn("removeFaces::compatibleRemoves")
                << "Problem2" << abort(FatalError);
        }

        label region0 = cellRegion[own];
        label region1 = cellRegion[nei];

        if (region0 == -1)
        {
            if (region1 == -1)
            {
                // Create new region
                cellRegion[own] = nRegions;
                cellRegion[nei] = nRegions;

                // Make owner (lowest numbered!) the master of the region
                regionMaster[nRegions] = own;
                nRegions++;
            }
            else
            {
                // Add owner to neighbour region
                cellRegion[own] = region1;
                // See if owner becomes the master of the region
                regionMaster[region1] = min(own, regionMaster[region1]);
            }
        }
        else
        {
            if (region1 == -1)
            {
                // Add neighbour to owner region
                cellRegion[nei] = region0;
                // nei is higher numbered than own so guaranteed not lower
                // than master of region0.
            }
            else if (region0 != region1)
            {
                // Both have regions. Keep lowest numbered region and master.
                label freedRegion = -1;
                label keptRegion = -1;

                if (region0 < region1)
                {
                    changeCellRegion
                    (
                        nei,
                        region1,    // old region
                        region0,    // new region
                        cellRegion
                    );

                    keptRegion = region0;
                    freedRegion = region1;
                }
                else if (region1 < region0)
                {
                    changeCellRegion
                    (
                        own,
                        region0,    // old region
                        region1,    // new region
                        cellRegion
                    );

                    keptRegion = region1;
                    freedRegion = region0;
                }

                label master0 = regionMaster[region0];
                label master1 = regionMaster[region1];

                regionMaster[freedRegion] = -1;
                regionMaster[keptRegion] = min(master0, master1);
            }
        }
    }

    regionMaster.setSize(nRegions);


    // Various checks
    // - master is lowest numbered in any region 
    // - regions have more than 1 cell
    {
        labelList nCells(regionMaster.size(), 0);

        forAll(cellRegion, cellI)
        {
            label r = cellRegion[cellI];

            if (r != -1)
            {
                nCells[r]++;

                if (cellI < regionMaster[r])
                {
                    FatalErrorIn
                    (
                        "removeFaces::compatibleRemoves(const labelList&"
                        ", labelList&, labelList&, labelList&)"
                    )   << "Not lowest numbered : cell:" << cellI
                        << " region:" << r
                        << " regionmaster:" << regionMaster[r]
                        << abort(FatalError);
                }
            }
        }

        forAll(nCells, region)
        {
            if (nCells[region] == 1)
            {
                FatalErrorIn
                (
                    "removeFaces::compatibleRemoves(const labelList&"
                    ", labelList&, labelList&, labelList&)"
                )   << "Region " << region
                    << " has only " << nCells[region] << " cells in it"
                    << abort(FatalError);
            }
            else if (nCells[region] != 2 && nCells[region] != 0)
            {
                WarningIn("removeFaces::compatibleRemoves")
                    << " region " << region << " has "
                    << nCells[region] << " cells in it" << endl
                    << "Removing will fail if faces share a common point!"
                    << endl;
            }
        }
    }




    // Count number of used regions
    label nUsedRegions = 0;

    forAll(regionMaster, i)
    {
        if (regionMaster[i] != -1)
        {
            nUsedRegions++;
        }
    }

    // Recreate facesToRemove to be consistent with the cellRegions.
    DynamicList<label> allFacesToRemove(facesToRemove.size());

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        if (cellRegion[own] != -1 && cellRegion[own] == cellRegion[nei])
        {
            // Both will become the same cell so add face to list of faces
            // to be removed.
            allFacesToRemove.append(faceI);
        }
    }
    allFacesToRemove.shrink();
    newFacesToRemove.transfer(allFacesToRemove);
    allFacesToRemove.clear();

    return nUsedRegions;
}


void Foam::removeFaces::setRefinement
(
    const labelList& faceLabels,
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    polyTopoChange& meshMod
) const
{



    // Make map of all faces to be removed
    labelHashSet facesToRemove(faceLabels);


    //
    // From faces to be removed get all edges affected. From these check the
    // ones that will have only 2 non-removed faces using them and see
    // if they can be merged.
    //

    // From edge that can be removed to the two faces that can be merged.
    Map<edge> edgesToRemove(2*faceLabels.size());

    {
        labelHashSet affectedEdges(4*faceLabels.size());

        forAll(faceLabels, labelI)
        {
            label faceI = faceLabels[labelI];

            const labelList& fEdges = mesh_.faceEdges()[faceI];

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                affectedEdges.insert(edgeI);
            }
        }


        for
        (
            labelHashSet::const_iterator iter = affectedEdges.begin();
            iter != affectedEdges.end();
            ++iter
        )
        {
            label edgeI = iter.key();

            label face0I = -1;
            label face1I = -1;

            if
            (
               !edgesToRemove.found(edgeI)
             && getUnremovedFaces
                (
                    facesToRemove,
                    edgeI,
                    face0I,
                    face1I
                )
            )
            {
                // Edge has two remaining faces. Check whether can be merged

                if (mesh_.isInternalFace(face0I))
                {
                    if (mesh_.isInternalFace(face1I))
                    {
                        // Internal faces. Mark edge for removal (regardless
                        // of angle) since the faces will share the same
                        // owner and neighbour so should be merged.
                        edgesToRemove.insert(edgeI, edge(face0I, face1I));
                    }
                }
                else if (!mesh_.isInternalFace(face1I))
                {
                    // Merge external face. Check angle.
                    if
                    (
                        mesh_.boundaryMesh().whichPatch(face0I)
                     == mesh_.boundaryMesh().whichPatch(face1I)
                    )
                    {
                        // Get angle between two faces. Since both external
                        // should have same owner so no check on ordering.
                        scalar cosAngle =
                            faceCos
                            (
                                face0I,
                                face1I,
                                false,
                                edgeI
                            );
                
                        if (cosAngle > minCos_)
                        {
                            // Faces on same patch and aligned. Merge.
                            edgesToRemove.insert(edgeI, edge(face0I, face1I));
                        }
                    }
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "edgesToRemove:" << edgesToRemove << endl;
    }


    //
    // From edges to be removed get all points affected. From these check the
    // ones that will have only 2 non-removed edges using them.
    //

    labelHashSet pointsToRemove(4*faceLabels.size());

    for
    (
        Map<edge>::const_iterator iter = edgesToRemove.begin();
        iter != edgesToRemove.end();
        ++iter
    )
    {
        label edgeI = iter.key();

        const edge& e = mesh_.edges()[edgeI];

        // Note order of checking to bypass expensive canBeRemoved check.
        if
        (
            !pointsToRemove.found(e.start())
         && canBeRemoved(edgesToRemove, e.start())
        )
        {
            pointsToRemove.insert(e.start());
        }

        if
        (
            !pointsToRemove.found(e.end())
         && canBeRemoved(edgesToRemove, e.end())
        )
        {
            pointsToRemove.insert(e.end());
        }
    }

    if (debug)
    {
        Pout<< "pointsToRemove:" << pointsToRemove << endl;
    }

    // Get all faces affected in any way by removal of points/edges/faces/cells
    labelHashSet facesAffected
    (
        getFacesAffected
        (
            cellRegion,
            cellRegionMaster,
            facesToRemove,
            edgesToRemove,
            pointsToRemove
        )
    );


    //
    // Now we know
    // - pointsToRemove : points to remove
    // - edgesToRemove  : edges to remove i.e. merge faces
    // - cellRegion/Master  : cells to remove
    // - facesAffected  : faces with points removed and/or owner/neighbour
    //                    changed
    

    // Start modifying mesh and keep track of faces changed.

    // Remove split faces.
    forAll(faceLabels, labelI)
    {
        label faceI = faceLabels[labelI];

        // Remove face if not yet uptodate (which is never; but want to be
        // consistent with rest of face removals/modifications)
        if (facesAffected.erase(faceI))
        {
            if (debug)
            {
                Pout<< "Actually removing split face " << faceI
                    << " verts:" << mesh_.faces()[faceI] << endl;
            }

            meshMod.setAction(polyRemoveFace(faceI));
        }
    }


    // Remove points.
    for
    (
        labelHashSet::const_iterator iter = pointsToRemove.begin();
        iter != pointsToRemove.end();
        ++iter
    )
    {
        label pointI = iter.key();

        if (debug)
        {
            Pout<< "Actually removing point " << pointI << " coord:"
                << mesh_.points()[pointI] << endl;
        }

        meshMod.setAction(polyRemovePoint(pointI));
    }


    // Remove cells.
    forAll(cellRegion, cellI)
    {
        label region = cellRegion[cellI];

        if (region != -1 && (cellI != cellRegionMaster[region]))
        {
            if (debug)
            {
                Pout<< "Actually removing cell " << cellI
                    << " master:" << cellRegionMaster[region]
                    << endl;
            }

            meshMod.setAction(polyRemoveCell(cellI));
        }
    }


    // Get the merged faces
    // - get all points in order of face0
    // - delete face using slave cell of face0 owner side
    // - modify face   ,,  master cell of face0 owner side.

    for
    (
        Map<edge>::const_iterator iter = edgesToRemove.begin();
        iter != edgesToRemove.end();
        ++iter
    )
    {
        label edgeI = iter.key();

        const edge& twoFaces = iter();

        label face0I = twoFaces.start();
        label face1I = twoFaces.end();


        // Remove face1 if not yet updated
        if (facesAffected.erase(face1I))
        {
            if (debug)
            {
                Pout<< "Actually removing merged face " << face1I
                    << " verts:" << mesh_.faces()[face1I]
                    << endl;
            }

            meshMod.setAction(polyRemoveFace(face1I));
        }


        // Modify face0 if not yet updated
        if (facesAffected.erase(face0I))
        {
            // Get merged face, oriented acc. to owner of face0
            face mergedFace
            (
                mergeFace
                (
                    pointsToRemove,
                    face0I,
                    face1I,
                    edgeI
                )
            );

            label own = mesh_.faceOwner()[face0I];

            if (cellRegion[own] != -1)
            {
                own = cellRegionMaster[cellRegion[own]];
            }


            label patchID, zoneID, zoneFlip;

            getFaceInfo(face0I, patchID, zoneID, zoneFlip);


            label nei = -1;

            if (mesh_.isInternalFace(face0I))
            {
                nei = mesh_.faceNeighbour()[face0I];

                if (cellRegion[nei] != -1)
                {
                    nei = cellRegionMaster[cellRegion[nei]];
                }
            }

            // So now we have the merged face and the remaining cells on both
            // sides. Merged face points are according to old owner.

            if (debug)
            {
                Pout<< "Modifying mergedface " << face0I << " for new verts:"
                    << mergedFace
                    << " possibly new owner " << own << " or new nei " << nei
                    << endl;
            }

            modFace
            (
                mergedFace,         // modified face
                face0I,             // label of face being modified
                own,                // owner
                nei,                // neighbour
                false,              // face flip
                patchID,            // patch for face
                false,              // remove from zone
                zoneID,             // zone for face
                zoneFlip,           // face flip in zone

                meshMod
            );
        }
    }


    // Check if any remaining faces have not been updated for new slave/master
    // or points removed.
    for
    (
        labelHashSet::const_iterator iter = facesAffected.begin();
        iter != facesAffected.end();
        ++iter
    )
    {
        label faceI = iter.key();

        if (debug)
        {
            Pout<< "Remaining affected face:" << faceI << endl;
        }

        face f = filterFace(pointsToRemove, faceI);

        label own = mesh_.faceOwner()[faceI];

        if (cellRegion[own] != -1)
        {
            own = cellRegionMaster[cellRegion[own]];
        }

        label patchID, zoneID, zoneFlip;

        getFaceInfo(faceI, patchID, zoneID, zoneFlip);

        label nei = -1;

        if (mesh_.isInternalFace(faceI))
        {
            nei = mesh_.faceNeighbour()[faceI];

            if (cellRegion[nei] != -1)
            {
                nei = cellRegionMaster[cellRegion[nei]];
            }
        }

        if (debug)
        {
            Pout<< "Modifying " << faceI << " for new verts:" << f
                << " or for new owner " << own << " or for new nei " << nei
                << endl;
        }

        modFace
        (
            f,                  // modified face
            faceI,              // label of face being modified
            own,                // owner
            nei,                // neighbour
            false,              // face flip
            patchID,            // patch for face
            false,              // remove from zone
            zoneID,             // zone for face
            zoneFlip,           // face flip in zone

            meshMod
        );
    }
}


// ************************************************************************* //
