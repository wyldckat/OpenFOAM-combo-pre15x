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

\*---------------------------------------------------------------------------*/

#include "directPolyTopoChange.H"
#include "SortableList.H"
#include "polyMesh.H"
#include "polyAddPoint.H"
#include "polyModifyPoint.H"
#include "polyRemovePoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyAddCell.H"
#include "polyModifyCell.H"
#include "polyRemoveCell.H"
#include "objectMap.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directPolyTopoChange, 0);
}



const Foam::point Foam::directPolyTopoChange::greatPoint
(
    GREAT,
    GREAT,
    GREAT
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Renumber
void Foam::directPolyTopoChange::renumber
(
    const labelList& map,
    DynamicList<label>& elems
)
{
    forAll(elems, elemI)
    {
        if (elems[elemI] >= 0)
        {
            elems[elemI] = map[elems[elemI]];
        }
    }
}

void Foam::directPolyTopoChange::renumber
(
    const labelList& map,
    labelHashSet& elems
)
{
    labelHashSet newElems(elems.size());

    forAllConstIter(labelHashSet, elems, iter)
    {
        label newElem = map[iter.key()];

        if (newElem >= 0)
        {
            newElems.insert(newElem);
        }
    }

    elems.transfer(newElems);
}


// Renumber and remove -1 elements.
void Foam::directPolyTopoChange::renumberCompact
(
    const labelList& map,
    labelList& elems
)
{
    label newElemI = 0;

    forAll(elems, elemI)
    {
        label newVal = map[elems[elemI]];

        if (newVal != -1)
        {
            elems[newElemI++] = newVal;
        }
    }
    elems.setSize(newElemI);
}


void Foam::directPolyTopoChange::checkFace
(
    const face& f,
    const label faceI,
    const label own,
    const label nei,
    const label patchI,
    const label zoneI
) const
{
    if (nei == -1)
    {
        if (own == -1 && zoneI != -1)
        {
            // retired face
        }
        else if (patchI == -1 || patchI >= nPatches_)
        {
            FatalErrorIn
            (
                "directPolyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Face has no neighbour (so external) but does not have"
                << " a valid patch" << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << abort(FatalError);
        }
    }
    else
    {
        if (patchI != -1)
        {
            FatalErrorIn
            (
                "directPolyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Cannot both have valid patchI and neighbour" << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << abort(FatalError);
        }

        if (nei <= own)
        {
            FatalErrorIn
            (
                "directPolyTopoChange::checkFace(const face&, const label"
                ", const label, const label, const label)"
            )   << "Owner cell label should be less than neighbour cell label"
                << nl
                << "f:" << f
                << " faceI(-1 if added face):" << faceI
                << " own:" << own << " nei:" << nei
                << " patchI:" << patchI << abort(FatalError);
        }
    }

    if (f.size() < 3 || findIndex(f, -1) != -1)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::checkFace(const face&, const label"
            ", const label, const label, const label)"
        )   << "Illegal vertices in face"
            << nl
            << "f:" << f
            << " faceI(-1 if added face):" << faceI
            << " own:" << own << " nei:" << nei
            << " patchI:" << patchI << abort(FatalError);
    }
}


void Foam::directPolyTopoChange::makeCells
(
    const label nActiveFaces,
    labelList& cellFaces,
    labelList& cellFaceOffsets
) const
{
    cellFaces.setSize(2*nActiveFaces);
    cellFaceOffsets.setSize(cellMap_.size() + 1);

    // Faces per cell
    labelList nNbrs(cellMap_.size(), 0);

    // 1. Count faces per cell

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        nNbrs[faceOwner_[faceI]]++;
    }
    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        if (faceNeighbour_[faceI] >= 0)
        {
            nNbrs[faceNeighbour_[faceI]]++;
        }
    }

    // 2. Calculate offsets

    cellFaceOffsets[0] = 0;
    forAll (nNbrs, cellI)
    {
        cellFaceOffsets[cellI+1] = cellFaceOffsets[cellI] + nNbrs[cellI];
    }

    // 3. Fill faces per cell

    // reset the whole list to use as counter
    nNbrs = 0;

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        label cellI = faceOwner_[faceI];

        cellFaces[cellFaceOffsets[cellI] + nNbrs[cellI]++] = faceI;
    }

    for (label faceI = 0; faceI < nActiveFaces; faceI++)
    {
        label cellI = faceNeighbour_[faceI];

        if (cellI >= 0)
        {
            cellFaces[cellFaceOffsets[cellI] + nNbrs[cellI]++] = faceI;
        }
    }

    // Last offset points to beyond end of cellFaces.
    cellFaces.setSize(cellFaceOffsets[cellMap_.size()]);
}


// Determine order for faces:
// - upper-triangular order for internal faces
// - external faces after internal faces and in patch order.
Foam::labelList Foam::directPolyTopoChange::getFaceOrder
(
    const label nActiveFaces,
    const labelList& cellFaces,
    const labelList& cellFaceOffsets
) const
{
    //const SubList<label> activeOwner(faceOwner_, nActiveFaces, 0);
    //const SubList<label> activeNeighbour(faceNeighbour_, nActiveFaces, 0);
    //const SubList<label> activeRegion(region_, nActiveFaces, 0);
    //
    //labelList newToOld(identity(nActiveFaces));
    //
    //// Sort according to owner so all faces of same cell are clustered.
    //Foam::sort
    //(
    //    newToOld,
    //    faceLess
    //    (
    //        activeOwner,
    //        activeNeighbour,
    //        activeRegion
    //    )
    //);
    //
    //labelList oldToNew(invert(faceOwner_.size(), newToOld));
    //
    //// Retired faces.
    //for (label faceI = nActiveFaces; faceI < oldToNew.size(); faceI++)
    //{
    //    oldToNew[faceI] = faceI;
    //}

    labelList oldToNew(faceOwner_.size(), -1);

    // First unassigned face
    label newFaceI = 0;

    forAll(cellMap_, cellI)
    {
        label startOfCell = cellFaceOffsets[cellI];
        label nFaces = cellFaceOffsets[cellI+1] - startOfCell;

        // Neighbouring cells
        SortableList<label> nbr(nFaces);

        for (label i = 0; i < nFaces; i++)
        {
            label faceI = cellFaces[startOfCell + i];

            label nbrCellI = faceNeighbour_[faceI];

            if (faceI >= nActiveFaces)
            {
                // Retired face.
                nbr[i] = -1;
            }
            else if (nbrCellI != -1)
            {
                // Internal face. Get cell on other side.
                if (nbrCellI == cellI)
                {
                    nbrCellI = faceOwner_[faceI];
                }

                if (cellI < nbrCellI)
                {
                    // CellI is master
                    nbr[i] = nbrCellI;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNew[cellFaces[startOfCell + nbr.indices()[i]]] =
                    newFaceI++;
            }
        }
    }


    // Pick up all patch faces in patch face order. (note: loops over all
    // faces for all patches. Not very efficient)
    for (label patchI = 0; patchI < nPatches_; patchI++)
    {
        for (label faceI = 0; faceI < nActiveFaces; faceI++)
        {
            if (region_[faceI] == patchI)
            {
                oldToNew[faceI] = newFaceI++;
            }
        }
    }

    // Retired faces.
    for (label faceI = nActiveFaces; faceI < oldToNew.size(); faceI++)
    {
        oldToNew[faceI] = faceI;
    }


    // Check done all faces.
    forAll(oldToNew, faceI)
    {
        if (oldToNew[faceI] == -1)
        {
            FatalErrorIn
            (
                "directPolyTopoChange::getFaceOrder"
                "(const label, const labelList&, const labelList&)"
                " const"
            )   << "Did not determine new position"
                << " for face " << faceI
                << abort(FatalError);
        }
    }

    return oldToNew;
}


// Compact and reorder faces according to map.
void Foam::directPolyTopoChange::reorderCompactFaces
(
    const label newSize,
    const labelList& oldToNew
)
{
    reorder(oldToNew, faces_);
    faces_.setSize(newSize);
    faces_.shrink();

    reorder(oldToNew, region_);
    region_.setSize(newSize);
    region_.shrink();

    reorder(oldToNew, faceOwner_);
    faceOwner_.setSize(newSize);
    faceOwner_.shrink();

    reorder(oldToNew, faceNeighbour_);
    faceNeighbour_.setSize(newSize);
    faceNeighbour_.shrink();

    // Update faceMaps.
    reorder(oldToNew, faceMap_);
    faceMap_.setSize(newSize);
    faceMap_.shrink();
    renumber(oldToNew, reverseFaceMap_);

    renumberKey(oldToNew, faceFromPoint_);
    renumberKey(oldToNew, faceFromEdge_);
    renumber(oldToNew, flipFaceFlux_);
    renumberKey(oldToNew, faceZone_);
    renumberKey(oldToNew, faceZoneFlip_);
}


// Compact all and orders faces:
// - internalfaces upper-triangular
// - externalfaces after internal ones.
void Foam::directPolyTopoChange::compact()
{
    points_.shrink();
    pointMap_.shrink();
    reversePointMap_.shrink();

    faces_.shrink();
    region_.shrink();
    faceOwner_.shrink();
    faceNeighbour_.shrink();
    faceMap_.shrink();
    reverseFaceMap_.shrink();

    cellMap_.shrink();
    reverseCellMap_.shrink();
    cellZone_.shrink();


    if (debug)
    {
        Pout<< "Before compacting:" << nl
            << "    points:" << points_.size() << nl
            << "    faces :" << faces_.size() << nl
            << "    cells :" << cellMap_.size() << nl
            << endl;
    }

    // Compact points
    label nActivePoints = -1;
    {
        labelList localPointMap(points_.size(), -1);
        label newPointI = 0;

        forAll(points_, pointI)
        {
            if (!pointRemoved(pointI) && !retiredPoints_.found(pointI))
            {
                localPointMap[pointI] = newPointI++;
            }
        }
        nActivePoints = newPointI;

        forAllConstIter(labelHashSet, retiredPoints_, iter)
        {
            localPointMap[iter.key()] = newPointI++;
        }

        if (debug)
        {
            Pout<< "Points : active:" << nActivePoints
                << "  retired:" << newPointI - nActivePoints
                << "  total:" << newPointI << endl;
        }

        reorder(localPointMap, points_);
        points_.setSize(newPointI);
        points_.shrink();

        // Update pointMaps
        reorder(localPointMap, pointMap_);
        pointMap_.setSize(newPointI);
        pointMap_.shrink();
        renumber(localPointMap, reversePointMap_);

        renumberKey(localPointMap, pointZone_);
        renumber(localPointMap, retiredPoints_);

        // Use map to relabel face vertices
        forAll(faces_, faceI)
        {
            face& f = faces_[faceI];

            labelList oldF(f);

            renumberCompact(localPointMap, f);

            if (!faceRemoved(faceI) && f.size() < 3)
            {
                FatalErrorIn("directPolyTopoChange::compact()")
                    << "Created illegal face " << f
                    << " from face " << oldF
                    << " at position:" << faceI
                    << " when filtering removed points"
                    << abort(FatalError);
            }
        }
    }


    // Compact faces.
    {
        labelList localFaceMap(faces_.size(), -1);
        label newFaceI = 0;

        forAll(faces_, faceI)
        {
            if (!faceRemoved(faceI) && faceOwner_[faceI] >= 0)
            {
                localFaceMap[faceI] = newFaceI++;
            }
        }
        nActiveFaces_ = newFaceI;

        forAll(faces_, faceI)
        {
            if (!faceRemoved(faceI) && faceOwner_[faceI] < 0)
            {
                // Retired face
                localFaceMap[faceI] = newFaceI++;
            }
        }

        if (debug)
        {
            Pout<< "Faces : active:" << nActiveFaces_
                << "  retired:" << newFaceI - nActiveFaces_
                << "  total:" << newFaceI << endl;
        }

        // Reorder faces.
        reorderCompactFaces(newFaceI, localFaceMap);
    }

    // Compact cells.
    {
        labelList localCellMap(cellMap_.size(), -1);
        label newCellI = 0;

        forAll(cellMap_, cellI)
        {
            if (!cellRemoved(cellI))      // removed cells get marked with -2
            {
                localCellMap[cellI] = newCellI++;
            }
        }

        if (newCellI != cellMap_.size())
        {
            reorder(localCellMap, cellMap_);
            cellMap_.setSize(newCellI);
            cellMap_.shrink();
            renumber(localCellMap, reverseCellMap_);

            reorder(localCellMap, cellZone_);
            cellZone_.setSize(newCellI);
            cellZone_.shrink();

            renumberKey(localCellMap, cellFromPoint_);
            renumberKey(localCellMap, cellFromEdge_);
            renumberKey(localCellMap, cellFromFace_);

            renumber(localCellMap, faceOwner_);
            renumber(localCellMap, faceNeighbour_);
        }

        if (debug)
        {
            Pout<< "Cells : (never retired):" << cellMap_.size() << endl;
        }
    }

    // Reorder faces into upper-triangular and patch ordering
    {
        // Create cells (packed storage)
        labelList cellFaces;
        labelList cellFaceOffsets;
        makeCells(nActiveFaces_, cellFaces, cellFaceOffsets);

        // Do upper triangular order.
        labelList localFaceMap
        (
            getFaceOrder
            (
                nActiveFaces_,
                cellFaces,
                cellFaceOffsets
            )
        );

        // Reorder faces.
        reorderCompactFaces(localFaceMap.size(), localFaceMap);
    }
}


void Foam::directPolyTopoChange::calcPatchSizes
(
    labelList& patchSizes,
    labelList& patchStarts
) const
{
    patchStarts.setSize(nPatches_);
    patchStarts = 0;
    patchSizes.setSize(nPatches_);
    patchSizes = 0;


    label& nInternalFaces = patchStarts[0];

    for (label faceI = 0; faceI < nActiveFaces_; faceI++)
    {
        if (region_[faceI] >= 0)
        {
            patchSizes[region_[faceI]]++;
        }
        else
        {
            nInternalFaces++;
        } 
    }

    label faceI = nInternalFaces;

    forAll(patchStarts, patchI)
    {
        patchStarts[patchI] = faceI;
        faceI += patchSizes[patchI];
    }

    if (debug)
    {
        Pout<< "patchSizes:" << patchSizes << nl
            << "patchStarts:" << patchStarts << endl;
    }
}


Foam::labelList Foam::directPolyTopoChange::selectFaces
(
    const primitiveMesh& mesh,
    const labelList& faceLabels,
    const bool internalFacesOnly
)
{
    label nFaces = 0;

    forAll(faceLabels, i)
    {
        label faceI = faceLabels[i];

        if (internalFacesOnly == mesh.isInternalFace(faceI))
        {
            nFaces++;
        }
    }

    labelList collectedFaces(nFaces);

    nFaces = 0;

    forAll(faceLabels, i)
    {
        label faceI = faceLabels[i];

        if (internalFacesOnly == mesh.isInternalFace(faceI))
        {
            collectedFaces[nFaces++] = faceI;
        }
    }

    return collectedFaces;
}


void Foam::directPolyTopoChange::calcPatchPointMap
(
    const List<Map<label> >& oldPatchMeshPointMaps,
    const labelList& oldPatchNMeshPoints,
    const polyBoundaryMesh& boundary,

    labelListList& patchPointRenumber
) const
{
    patchPointRenumber.setSize(boundary.size());

    forAll(boundary, patchI)
    {
        const labelList& newPatchMeshPoints = boundary[patchI].meshPoints();

        const Map<label>& oldZoneMeshPointMap = oldPatchMeshPointMaps[patchI];
        const label oldSize = oldPatchNMeshPoints[patchI];

        labelList& curPatchPointRnb = patchPointRenumber[patchI];

        curPatchPointRnb.setSize(newPatchMeshPoints.size());

        forAll(newPatchMeshPoints, pointI)
        {
            if (newPatchMeshPoints[pointI] < oldSize)
            {
                Map<label>::const_iterator ozmpmIter =
                    oldZoneMeshPointMap.find
                    (
                        pointMap_[newPatchMeshPoints[pointI]]
                    );

                if (ozmpmIter != oldZoneMeshPointMap.end())
                {
                    curPatchPointRnb[pointI] = ozmpmIter();
                }
                else
                {
                    curPatchPointRnb[pointI] = -1;
                }
            }
            else
            {
                curPatchPointRnb[pointI] = -1;
            }
        }
    }
}


void Foam::directPolyTopoChange::calcFaceInflationMaps
(
    const polyMesh& mesh,
    List<objectMap>& facesFromPoints,
    List<objectMap>& facesFromEdges
) const
{
    // Faces inflated from points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    facesFromPoints.setSize(faceFromPoint_.size());

    if (faceFromPoint_.size() > 0)
    {
        label nFacesFromPoints = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, faceFromPoint_, iter)
        {
            label newFaceI = iter.key();

            if (region_[newFaceI] == -1)
            {
                // Get internal faces using point on old mesh
                facesFromPoints[nFacesFromPoints++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.pointFaces()[iter()],
                        true
                    )
                );
            }
            else
            {
                // Get patch faces using point on old mesh
                facesFromPoints[nFacesFromPoints++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.pointFaces()[iter()],
                        false
                    )
                );
            }
        }
    }


    // Faces inflated from edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    facesFromEdges.setSize(faceFromEdge_.size());

    if (faceFromEdge_.size() > 0)
    {
        label nFacesFromEdges = 0;

        // Collect all still existing faces connected to this edge.
        forAllConstIter(Map<label>, faceFromEdge_, iter)
        {
            label newFaceI = iter.key();

            if (region_[newFaceI] == -1)
            {
                // Get internal faces using edge on old mesh
                facesFromEdges[nFacesFromEdges++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.edgeFaces()[iter()],
                        true
                    )
                );
            }
            else
            {
                // Get patch faces using edge on old mesh
                facesFromEdges[nFacesFromEdges++] = objectMap
                (
                    newFaceI,
                    selectFaces
                    (
                        mesh,
                        mesh.edgeFaces()[iter()],
                        false
                    )
                );
            }
        }
    }
}


void Foam::directPolyTopoChange::calcCellInflationMaps
(
    const polyMesh& mesh,
    List<objectMap>& cellsFromPoints,
    List<objectMap>& cellsFromEdges,
    List<objectMap>& cellsFromFaces
) const
{
    cellsFromPoints.setSize(cellFromPoint_.size());

    if (cellFromPoint_.size() > 0)
    {
        label nCellsFromPoints = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, cellFromPoint_, iter)
        {
            cellsFromPoints[nCellsFromPoints++] = objectMap
            (
                iter.key(),
                mesh.pointCells()[iter()]
            );
        }
    }


    cellsFromEdges.setSize(cellFromEdge_.size());

    if (cellFromEdge_.size() > 0)
    {
        label nCellsFromEdges = 0;

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, cellFromEdge_, iter)
        {
            cellsFromEdges[nCellsFromEdges++] = objectMap
            (
                iter.key(),
                mesh.edgeCells()[iter()]
            );
        }
    }


    cellsFromFaces.setSize(cellFromFace_.size());

    if (cellFromFace_.size() > 0)
    {
        label nCellsFromFaces = 0;

        labelList cellsAroundFace(2, -1);

        // Collect all still existing faces connected to this point.
        forAllConstIter(Map<label>, cellFromFace_, iter)
        {
            label oldFaceI = iter();

            if (mesh.isInternalFace(oldFaceI))
            {
                cellsAroundFace[0] = mesh.faceOwner()[oldFaceI];
                cellsAroundFace[1] = mesh.faceNeighbour()[oldFaceI];
            }
            else
            {
                cellsAroundFace[0] = mesh.faceOwner()[oldFaceI];
                cellsAroundFace[1] = -1;
            }

            cellsFromFaces[nCellsFromFaces++] = objectMap
            (
                iter.key(),
                cellsAroundFace
            );
        }
    }
}


void Foam::directPolyTopoChange::resetZones
(
    polyMesh& mesh,
    labelListList& pointZoneMap,
    labelListList& faceZoneFaceMap,
    labelListList& cellZoneMap
) const
{
    // pointZones
    // ~~~~~~~~~~

    pointZoneMap.setSize(mesh.pointZones().size());
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        // Count points per zone

        labelList nPoints(pointZones.size(), 0);

        forAllConstIter(Map<label>, pointZone_, iter)
        {
            label zoneI = iter();

            if (zoneI < 0 || zoneI >= pointZones.size())
            {
                FatalErrorIn
                (
                    "resetZones(polyMesh& mesh, labelListList&"
                    "labelListList&, labelListList&)"
                )   << "Illegal zoneID " << zoneI << " for point "
                    << iter.key() << " coord " << mesh.allPoints()[iter.key()]
                    << abort(FatalError);
            }
            nPoints[zoneI]++;
        }

        // Distribute points per zone

        labelListList addressing(pointZones.size());
        forAll(addressing, zoneI)
        {
            addressing[zoneI].setSize(nPoints[zoneI]);
        }
        nPoints = 0;

        forAllConstIter(Map<label>, pointZone_, iter)
        {
            label zoneI = iter();

            addressing[zoneI][nPoints[zoneI]++] = iter.key();
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get pointZoneMap.
        forAll(addressing, zoneI)
        {
            const pointZone& oldZone = pointZones[zoneI];
            const labelList& newZoneAddr = addressing[zoneI];

            labelList& curPzRnb = pointZoneMap[zoneI];
            curPzRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < pointMap_.size())
                {
                    curPzRnb[i] =
                        oldZone.whichPoint(pointMap_[newZoneAddr[i]]);
                }
                else
                {
                    curPzRnb[i] = -1;
                }
            }
        }

        // Reset the addresing on the zone
        mesh.pointZones().clearAddressing();
        forAll(mesh.pointZones(), zoneI)
        {
            if (debug)
            {
                Pout<< "pointZone:" << zoneI
                    << "  name:" << mesh.pointZones()[zoneI].name()
                    << "  size:" << addressing[zoneI].size()
                    << endl;
            }

            mesh.pointZones()[zoneI] = addressing[zoneI];
        }
    }


    // faceZones
    // ~~~~~~~~~

    faceZoneFaceMap.setSize(mesh.faceZones().size());
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList nFaces(faceZones.size(), 0);

        forAllConstIter(Map<label>, faceZone_, iter)
        {
            label zoneI = iter();

            if (zoneI < 0 || zoneI >= faceZones.size())
            {
                FatalErrorIn
                (
                    "resetZones(polyMesh& mesh, labelListList&"
                    "labelListList&, labelListList&)"
                )   << "Illegal zoneID " << zoneI << " for face "
                    << iter.key()
                    << abort(FatalError);
            }
            nFaces[zoneI]++;
        }

        labelListList addressing(faceZones.size());
        boolListList flipMode(faceZones.size());

        forAll(addressing, zoneI)
        {
            addressing[zoneI].setSize(nFaces[zoneI]);
            flipMode[zoneI].setSize(nFaces[zoneI]);
        }
        nFaces = 0;

        forAllConstIter(Map<label>, faceZone_, iter)
        {
            label zoneI = iter();
            label faceI = iter.key();

            label index = nFaces[zoneI]++;

            addressing[zoneI][index] = faceI;
            flipMode[zoneI][index] = faceZoneFlip_[faceI];
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get faceZoneFaceMap.
        forAll(addressing, zoneI)
        {
            const faceZone& oldZone = faceZones[zoneI];
            const labelList& newZoneAddr = addressing[zoneI];

            labelList& curFzFaceRnb = faceZoneFaceMap[zoneI];

            curFzFaceRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < faceMap_.size())
                {
                    curFzFaceRnb[i] =
                        oldZone.whichFace(faceMap_[newZoneAddr[i]]);
                }
                else
                {
                    curFzFaceRnb[i] = -1;
                }
            }
        }


        // Reset the addresing on the zone
        mesh.faceZones().clearAddressing();
        forAll(mesh.faceZones(), zoneI)
        {
            if (debug)
            {
                Pout<< "faceZone:" << zoneI
                    << "  name:" << mesh.faceZones()[zoneI].name()
                    << "  size:" << addressing[zoneI].size()
                    << endl;
            }

            mesh.faceZones()[zoneI].resetAddressing
            (
                addressing[zoneI],
                flipMode[zoneI]
            );
        }
    }


    // cellZones
    // ~~~~~~~~~

    cellZoneMap.setSize(mesh.cellZones().size());
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList nCells(cellZones.size(), 0);

        forAll(cellZone_, cellI)
        {
            label zoneI = cellZone_[cellI];

            if (zoneI >= cellZones.size())
            {
                FatalErrorIn
                (
                    "resetZones(polyMesh& mesh, labelListList&"
                    "labelListList&, labelListList&)"
                )   << "Illegal zoneID " << zoneI << " for cell "
                    << cellI << abort(FatalError);
            }

            if (zoneI >= 0)
            {
                nCells[zoneI]++;
            }
        }

        labelListList addressing(cellZones.size());
        forAll(addressing, zoneI)
        {
            addressing[zoneI].setSize(nCells[zoneI]);
        }
        nCells = 0;

        forAll(cellZone_, cellI)
        {
            label zoneI = cellZone_[cellI];

            if (zoneI >= 0)
            {
                addressing[zoneI][nCells[zoneI]++] = cellI;
            }
        }

        // So now we both have old zones and the new addressing.
        // Invert the addressing to get cellZoneMap.
        forAll(addressing, zoneI)
        {
            const cellZone& oldZone = cellZones[zoneI];
            const labelList& newZoneAddr = addressing[zoneI];

            labelList& curCellRnb = cellZoneMap[zoneI];

            curCellRnb.setSize(newZoneAddr.size());

            forAll(newZoneAddr, i)
            {
                if (newZoneAddr[i] < cellMap_.size())
                {
                    curCellRnb[i] =
                        oldZone.whichCell(cellMap_[newZoneAddr[i]]);
                }
                else
                {
                    curCellRnb[i] = -1;
                }
            }
        }

        // Reset the addresing on the zone
        mesh.cellZones().clearAddressing();
        forAll(mesh.cellZones(), zoneI)
        {
            if (debug)
            {
                Pout<< "cellZone:" << zoneI
                    << "  name:" << mesh.cellZones()[zoneI].name()
                    << "  size:" << addressing[zoneI].size()
                    << endl;
            }

            mesh.cellZones()[zoneI] = addressing[zoneI];
        }
    }
}


void Foam::directPolyTopoChange::calcFaceZonePointMap
(
    polyMesh& mesh,
    const List<Map<label> >& oldFaceZoneMeshPointMaps,
    labelListList& faceZonePointMap
) const
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    faceZonePointMap.setSize(faceZones.size());

    forAll(faceZones, zoneI)
    {
        const faceZone& newZone = faceZones[zoneI];

        const labelList& newZoneMeshPoints = newZone().meshPoints();

        const Map<label>& oldZoneMeshPointMap = oldFaceZoneMeshPointMaps[zoneI];

        labelList& curFzPointRnb = faceZonePointMap[zoneI];

        curFzPointRnb.setSize(newZoneMeshPoints.size());

        forAll(newZoneMeshPoints, pointI)
        {
            if (newZoneMeshPoints[pointI] < pointMap_.size())
            {
                Map<label>::const_iterator ozmpmIter =
                    oldZoneMeshPointMap.find
                    (
                        pointMap_[newZoneMeshPoints[pointI]]
                    );

                if (ozmpmIter != oldZoneMeshPointMap.end())
                {
                    curFzPointRnb[pointI] = ozmpmIter();
                }
                else
                {
                    curFzPointRnb[pointI] = -1;
                }
            }
            else
            {
                curFzPointRnb[pointI] = -1;
            }
        }
    }
}


Foam::face Foam::directPolyTopoChange::rotateFace
(
    const face& f,
    const label nPos
)
{
    face newF(f.size());

    forAll(f, fp)
    {
        label fp1 = (fp + nPos) % f.size();

        if (fp1 < 0)
        {
            fp1 += f.size();
        }

        newF[fp1] = f[fp];
    }

    return newF;
}


void Foam::directPolyTopoChange::reorderCoupledFaces
(
    const polyBoundaryMesh& boundary,
    const labelList& patchStarts,
    const labelList& patchSizes,
    const pointField& points
)
{
    // Mapping for faces (old to new). Extends over all mesh faces for
    // convenience (could be just the external faces)
    labelList oldToNew(faces_.size());
    forAll(oldToNew, faceI)
    {
        oldToNew[faceI] = faceI;
    }

    // Rotation on new faces.
    labelList rotation(faces_.size(), 0);

    // Send ordering
    forAll(boundary, patchI)
    {
        boundary[patchI].initOrder
        (
            primitivePatch
            (
                SubList<face>
                (
                    faces_,
                    patchSizes[patchI],
                    patchStarts[patchI]
                ),
                points
            )
        );
    }

    // Receive and calculate ordering

    bool anyChanged = false;

    forAll(boundary, patchI)
    {
        labelList patchFaceMap(patchSizes[patchI], -1);
        labelList patchFaceRotation(patchSizes[patchI], 0);

        bool changed = boundary[patchI].order
        (
            primitivePatch
            (
                SubList<face>
                (
                    faces_,
                    patchSizes[patchI],
                    patchStarts[patchI]
                ),
                points
            ),
            patchFaceMap,
            patchFaceRotation
        );

        if (changed)
        {
            // Merge patch face reordering into mesh face reordering table
            label start = patchStarts[patchI];

            forAll(patchFaceMap, patchFaceI)
            {
                oldToNew[patchFaceI + start] = start + patchFaceMap[patchFaceI];
            }

            forAll(patchFaceRotation, patchFaceI)
            {
                rotation[patchFaceI + start] = patchFaceRotation[patchFaceI];
            }

            anyChanged = true;
        }
    }

    reduce(anyChanged, orOp<bool>());

    if (anyChanged)
    {
        // Reorder faces according to oldToNew.
        reorderCompactFaces(oldToNew.size(), oldToNew);

        // Rotate faces (rotation is already in new face indices).
        forAll(rotation, faceI)
        {
            if (rotation[faceI] != 0)
            {
                faces_[faceI] = rotateFace(faces_[faceI], rotation[faceI]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::directPolyTopoChange::directPolyTopoChange
(
    const polyMesh& mesh,
    const bool strict
)
:
    strict_(strict),
    nPatches_(0),
    points_(0),
    pointMap_(0),
    pointZone_(0),
    retiredPoints_(0),
    faces_(0),
    region_(0),
    faceOwner_(0),
    faceNeighbour_(0),
    faceMap_(0),
    faceFromPoint_(0),
    faceFromEdge_(0),
    flipFaceFlux_(0),
    faceZone_(0),
    faceZoneFlip_(0),
    nActiveFaces_(0),
    cellMap_(0),
    cellFromPoint_(0),
    cellFromEdge_(0),
    cellFromFace_(0),
    cellZone_(0)
{
    addMesh
    (
        mesh,
        identity(mesh.boundaryMesh().size()),
        identity(mesh.pointZones().size()),
        identity(mesh.faceZones().size()),
        identity(mesh.cellZones().size())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directPolyTopoChange::clear()
{
    points_.clear();
    points_.setSize(0);
    pointMap_.clear();
    pointMap_.setSize(0);
    reversePointMap_.clear();
    reversePointMap_.setSize(0);
    pointZone_.clear();
    pointZone_.resize(0);
    retiredPoints_.clear();
    retiredPoints_.resize(0);

    faces_.clear();
    faces_.setSize(0);
    region_.clear();
    region_.setSize(0);
    faceOwner_.clear();
    faceOwner_.setSize(0);
    faceNeighbour_.clear();
    faceNeighbour_.setSize(0);
    faceMap_.clear();
    faceMap_.setSize(0);
    reverseFaceMap_.clear();
    reverseFaceMap_.setSize(0);
    faceFromPoint_.clear();
    faceFromPoint_.resize(0);
    faceFromEdge_.clear();
    faceFromEdge_.resize(0);
    flipFaceFlux_.clear();
    flipFaceFlux_.resize(0);
    faceZone_.clear();
    faceZone_.resize(0);
    faceZoneFlip_.clear();
    faceZoneFlip_.resize(0);
    nActiveFaces_ = 0;

    cellMap_.clear();
    cellMap_.setSize(0);    
    reverseCellMap_.clear();
    reverseCellMap_.setSize(0);    
    cellZone_.clear();
    cellZone_.setSize(0);
    cellFromPoint_.clear();
    cellFromPoint_.resize(0);
    cellFromEdge_.clear();
    cellFromEdge_.resize(0);
    cellFromFace_.clear();
    cellFromFace_.resize(0);
}


void Foam::directPolyTopoChange::addMesh
(
    const polyMesh& mesh,
    const labelList& patchMap,
    const labelList& pointZoneMap,
    const labelList& faceZoneMap,
    const labelList& cellZoneMap
)
{
    label maxRegion = nPatches_ - 1;
    forAll(patchMap, i)
    {
        maxRegion = max(maxRegion, patchMap[i]);
    }
    nPatches_ = maxRegion + 1;


    // Add points
    {
        const pointField& points = mesh.allPoints();
        const pointZoneMesh& pointZones = mesh.pointZones();

        // Resize
        points_.setSize(points_.size() + points.size());
        pointMap_.setSize(pointMap_.size() + points.size());
        pointZone_.resize(pointZone_.size() + points.size()/100);

        // Precalc offset zones
        labelList newZoneID(points.size(), -1);

        forAll(pointZones, zoneI)
        {
            const labelList& pointLabels = pointZones[zoneI];

            forAll(pointLabels, j)
            {
                newZoneID[pointLabels[j]] = pointZoneMap[zoneI];
            }
        }

        // Add points in mesh order
        forAll(points, pointI)
        {
            addPoint
            (
                points[pointI],
                pointI,
                newZoneID[pointI],
                true
            );
        }
    }

    // Add cells
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        // Resize

        // Note: polyMesh does not allow retired cells anymore. So allCells
        // always equals nCells
        label nAllCells = mesh.nCells();

        cellMap_.setSize(cellMap_.size() + nAllCells);
        cellFromPoint_.resize(cellFromPoint_.size() + nAllCells/100);
        cellFromEdge_.resize(cellFromEdge_.size() + nAllCells/100);
        cellFromFace_.resize(cellFromFace_.size() + nAllCells/100);
        cellZone_.setSize(cellZone_.size() + nAllCells);


        // Precalc offset zones
        labelList newZoneID(nAllCells, -1);
        
        forAll(cellZones, zoneI)
        {
            const labelList& cellLabels = cellZones[zoneI];

            forAll(cellLabels, j)
            {
                newZoneID[cellLabels[j]] = cellZoneMap[zoneI];
            }
        }

        // Add cells in mesh order
        for (label cellI = 0; cellI < nAllCells; cellI++)
        {
            // Add cell from cell
            addCell(-1, -1, -1, cellI, newZoneID[cellI]);
        }
    }

    // Add faces
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const faceList& faces = mesh.allFaces();
        const labelList& faceOwner = mesh.allOwner();
        const labelList& faceNeighbour = mesh.allNeighbour();
        const faceZoneMesh& faceZones = mesh.faceZones();

        // Resize
        label nAllFaces = mesh.allFaces().size();

        faces_.setSize(faces_.size() + nAllFaces);
        region_.setSize(region_.size() + nAllFaces);
        faceOwner_.setSize(faceOwner_.size() + nAllFaces);
        faceNeighbour_.setSize(faceNeighbour_.size() + nAllFaces);
        faceMap_.setSize(faceMap_.size() + nAllFaces);
        faceFromPoint_.resize(faceFromPoint_.size() + nAllFaces/100);
        faceFromEdge_.resize(faceFromEdge_.size() + nAllFaces/100);
        flipFaceFlux_.resize(flipFaceFlux_.size() + nAllFaces/100);
        faceZone_.resize(faceZone_.size() + nAllFaces/100);
        faceZoneFlip_.resize(faceZoneFlip_.size() + nAllFaces/100);


        // Precalc offset zones
        labelList newZoneID(nAllFaces, -1);
        boolList zoneFlip(nAllFaces, false);

        forAll(faceZones, zoneI)
        {
            const labelList& faceLabels = faceZones[zoneI];
            const boolList& flipMap = faceZones[zoneI].flipMap();

            forAll(faceLabels, j)
            {
                newZoneID[faceLabels[j]] = faceZoneMap[zoneI];
                zoneFlip[faceLabels[j]] = flipMap[j];
            }
        }

        // Add faces in mesh order

        // 1. Internal faces
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            addFace
            (
                faces[faceI],
                faceOwner[faceI],
                faceNeighbour[faceI],
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceI,                      // masterFaceID
                false,                      // flipFaceFlux
                -1,                         // patchID
                newZoneID[faceI],           // zoneID
                zoneFlip[faceI]             // zoneFlip
            );
        }

        // 2. Patch faces
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.start() != faces_.size())
            {
                FatalErrorIn
                (
                    "directPolyTopoChange::directPolyTopoChange"
                    "(const polyMesh& mesh, const bool strict)"
                )   << "Problem : "
                    << "Patch " << pp.name() << " starts at " << pp.start()
                    << endl
                    << "Current face counter at " << faces_.size() << endl
                    << "Are patches in incremental order?"
                    << abort(FatalError);
            }
            forAll(pp, patchFaceI)
            {
                label faceI = pp.start() + patchFaceI;

                addFace
                (
                    faces[faceI],
                    faceOwner[faceI],
                    -1,                         // neighbour
                    -1,                         // masterPointID
                    -1,                         // masterEdgeID
                    faceI,                      // masterFaceID
                    false,                      // flipFaceFlux
                    patchMap[patchI],           // patchID
                    newZoneID[faceI],           // zoneID
                    zoneFlip[faceI]             // zoneFlip
                );
            }
        }

        // 3. Retired faces
        for (label faceI = mesh.nFaces(); faceI < nAllFaces; faceI++)
        {
            addFace
            (
                faces[faceI],
                -1,
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                faceI,                      // masterFaceID
                false,                      // flipFaceFlux
                -1,                         // patchID
                newZoneID[faceI],           // zoneID
                zoneFlip[faceI]             // zoneFlip
            );
        }
    }
}


Foam::label Foam::directPolyTopoChange::setAction(const topoAction& action)
{
    if (isType<polyAddPoint>(action))
    {
        const polyAddPoint& pap = refCast<const polyAddPoint>(action);

        return addPoint
        (
            pap.newPoint(),
            pap.masterPointID(),
            pap.zoneID(),
            pap.inCell()
        );
    }
    else if (isType<polyModifyPoint>(action))
    {
        const polyModifyPoint& pmp = refCast<const polyModifyPoint>(action);

        modifyPoint
        (
            pmp.pointID(),
            pmp.newPoint(),
            pmp.zoneID(),
            pmp.inCell()
        );

        return -1;
    }
    else if (isType<polyRemovePoint>(action))
    {
        const polyRemovePoint& prp = refCast<const polyRemovePoint>(action);

        removePoint(prp.pointID());

        return -1;
    }
    else if (isType<polyAddFace>(action))
    {
        const polyAddFace& paf = refCast<const polyAddFace>(action);

        return addFace
        (
            paf.newFace(),
            paf.owner(),
            paf.neighbour(),
            paf.masterPointID(),
            paf.masterEdgeID(),
            paf.masterFaceID(),
            paf.flipFaceFlux(),
            paf.patchID(),
            paf.zoneID(),
            paf.zoneFlip()
        );
    }
    else if (isType<polyModifyFace>(action))
    {
        const polyModifyFace& pmf = refCast<const polyModifyFace>(action);

        modifyFace
        (
            pmf.newFace(),
            pmf.faceID(),
            pmf.owner(),
            pmf.neighbour(),
            pmf.flipFaceFlux(),
            pmf.patchID(),
            pmf.zoneID(),
            pmf.zoneFlip()
        );

        return -1;
    }
    else if (isType<polyRemoveFace>(action))
    {
        const polyRemoveFace& prf = refCast<const polyRemoveFace>(action);

        removeFace(prf.faceID());

        return -1;
    }
    else if (isType<polyAddCell>(action))
    {
        const polyAddCell& pac = refCast<const polyAddCell>(action);

        return addCell
        (
            pac.masterPointID(),
            pac.masterEdgeID(),
            pac.masterFaceID(),
            pac.masterCellID(),
            pac.zoneID()
        );
    }
    else if (isType<polyModifyCell>(action))
    {
        const polyModifyCell& pmc = refCast<const polyModifyCell>(action);

        if (pmc.removeFromZone())
        {
            modifyCell(pmc.cellID(), -1);
        }
        else
        {
            modifyCell(pmc.cellID(), pmc.zoneID());
        }

        return -1;
    }
    else if (isType<polyRemoveCell>(action))
    {
        const polyRemoveCell& prc = refCast<const polyRemoveCell>(action);

        removeCell(prc.cellID());

        return -1;
    }
    else
    {
        FatalErrorIn
        (
            "label directPolyTopoChange::setAction(const topoAction& action)"
        )   << "Unknown type of topoChange: " << action.type()
            << abort(FatalError);

        // Dummy return to keep compiler happy
        return -1;
    }
}


Foam::label Foam::directPolyTopoChange::addPoint
(
    const point& pt,
    const label masterPointID,
    const label zoneID,
    const bool inCell
)
{
    label pointI = points_.size();

    points_.append(pt);
    pointMap_.append(masterPointID);
    reversePointMap_.append(pointI);

    if (zoneID >= 0)
    {
        pointZone_.insert(pointI, zoneID);
    }

    if (!inCell)
    {
        retiredPoints_.insert(pointI);
    }

    return pointI;
}


void Foam::directPolyTopoChange::modifyPoint
(
    const label pointI,
    const point& pt,
    const label newZoneID,
    const bool inCell
)
{
    if (pointI < 0 || pointI >= points_.size())
    {
        FatalErrorIn
        (
            "directPolyTopoChange::modifyPoint(const label, const point&)"
        )   << "illegal point label " << pointI << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }
    if (pointRemoved(pointI) || pointMap_[pointI] == -1)
    {
        FatalErrorIn
        (
            "directPolyTopoChange::modifyPoint(const label, const point&)"
        )   << "point " << pointI << " already marked for removal"
            << abort(FatalError);
    }
    points_[pointI] = pt;

    Map<label>::iterator pointFnd = pointZone_.find(pointI);

    if (pointFnd != pointZone_.end())
    {
        if (newZoneID >= 0)
        {
            pointFnd() = newZoneID;
        }
        else
        {
            pointZone_.erase(pointFnd);
        }
    }
    else if (newZoneID >= 0)
    {
        pointZone_.insert(pointI, newZoneID);
    }

    if (inCell)
    {
        retiredPoints_.erase(pointI);
    }
    else
    {
        retiredPoints_.insert(pointI);
    }
}


//void Foam::directPolyTopoChange::movePoints(const pointField& newPoints)
//{
//    if (newPoints.size() != points_.size())
//    {
//        FatalErrorIn("directPolyTopoChange::movePoints(const pointField&)")
//            << "illegal pointField size." << endl
//            << "Size:" << newPoints.size() << endl
//            << "Points in mesh:" << points_.size()
//            << abort(FatalError);
//    }
//
//    forAll(points_, pointI)
//    {
//        points_[pointI] = newPoints[pointI];
//    }
//}


void Foam::directPolyTopoChange::removePoint(const label pointI)
{
    if (pointI < 0 || pointI >= points_.size())
    {
        FatalErrorIn("directPolyTopoChange::removePoint(const label)")
            << "illegal point label " << pointI << endl
            << "Valid point labels are 0 .. " << points_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (pointRemoved(pointI) || pointMap_[pointI] == -1)
    )
    {
        FatalErrorIn("directPolyTopoChange::removePoint(const label)")
            << "point " << pointI << " already marked for removal" << nl
            << "Point:" << points_[pointI] << " pointMap:" << pointMap_[pointI]
            << abort(FatalError);
    }

    points_[pointI] = greatPoint;
    pointMap_[pointI] = -1;
    reversePointMap_[pointI] = -1;
    pointZone_.erase(pointI);
    retiredPoints_.erase(pointI);
}


Foam::label Foam::directPolyTopoChange::addFace
(
    const face& f,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const bool flipFaceFlux,
    const label patchID,
    const label zoneID,
    const bool zoneFlip
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, -1, own, nei, patchID, zoneID);
    }

    label faceI = faces_.size();

    faces_.append(f);
    region_.append(patchID);
    faceOwner_.append(own);
    faceNeighbour_.append(nei);

    if (masterPointID >= 0)
    {
        faceMap_.append(-1);
        faceFromPoint_.insert(faceI, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        faceMap_.append(-1);
        faceFromEdge_.insert(faceI, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        faceMap_.append(masterFaceID);
    }
    else
    {
        // Allow inflate-from-nothing?
        //FatalErrorIn("directPolyTopoChange::addFace")
        //    << "Need to specify a master point, edge or face"
        //    << "face:" << f << " own:" << own << " nei:" << nei
        //    << abort(FatalError);
        faceMap_.append(-1);
    }
    reverseFaceMap_.append(faceI);

    if (flipFaceFlux)
    {
        flipFaceFlux_.insert(faceI);
    }

    if (zoneID >= 0)
    {
        faceZone_.insert(faceI, zoneID);
        faceZoneFlip_.insert(faceI, zoneFlip);
    }

    return faceI;
}    


void Foam::directPolyTopoChange::modifyFace
(
    const face& f,
    const label faceI,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label patchID,
    const label zoneID,
    const bool zoneFlip
)
{
    // Check validity
    if (debug)
    {
        checkFace(f, faceI, own, nei, patchID, zoneID);
    }

    faces_[faceI] = f;
    faceOwner_[faceI] = own;
    faceNeighbour_[faceI] = nei;
    region_[faceI] = patchID;

    if (flipFaceFlux)
    {
        flipFaceFlux_.insert(faceI);
    }
    else
    {
        flipFaceFlux_.erase(faceI);
    }

    Map<label>::iterator faceFnd = faceZone_.find(faceI);

    if (faceFnd != faceZone_.end())
    {
        if (zoneID >= 0)
        {
            faceFnd() = zoneID;
            faceZoneFlip_.find(faceI)() = zoneFlip;
        }
        else
        {
            faceZone_.erase(faceFnd);
            faceZoneFlip_.erase(faceI);
        }
    }
    else if (zoneID >= 0)
    {
        faceZone_.insert(faceI, zoneID);
        faceZoneFlip_.insert(faceI, zoneFlip);
    }
}    


void Foam::directPolyTopoChange::removeFace(const label faceI)
{
    if (faceI < 0 || faceI >= faces_.size())
    {
        FatalErrorIn("directPolyTopoChange::removeFace(const label)")
            << "illegal face label " << faceI << endl
            << "Valid face labels are 0 .. " << faces_.size()-1
            << abort(FatalError);
    }

    if
    (
        strict_
     && (faceRemoved(faceI) || faceMap_[faceI] == -1)
    )
    {
        FatalErrorIn("directPolyTopoChange::removeFace(const label)")
            << "face " << faceI
            << " already marked for removal"
            << abort(FatalError);
    }

    faces_[faceI].setSize(0);
    region_[faceI] = -1;
    faceOwner_[faceI] = -1;
    faceNeighbour_[faceI] = -1;
    faceMap_[faceI] = -1;
    reverseFaceMap_[faceI] = -1;
    faceFromEdge_.erase(faceI);
    faceFromPoint_.erase(faceI);
    flipFaceFlux_.erase(faceI);
    faceZone_.erase(faceI);
    faceZoneFlip_.erase(faceI);
}


Foam::label Foam::directPolyTopoChange::addCell
(
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label masterCellID,
    const label zoneID
)
{
    label cellI = cellMap_.size();

    if (masterPointID >= 0)
    {
        cellMap_.append(-1);
        cellFromPoint_.insert(cellI, masterPointID);
    }
    else if (masterEdgeID >= 0)
    {
        cellMap_.append(-1);
        cellFromEdge_.insert(cellI, masterEdgeID);
    }
    else if (masterFaceID >= 0)
    {
        cellMap_.append(-1);
        cellFromFace_.insert(cellI, masterFaceID);
    }
    else
    {
        cellMap_.append(masterCellID);
    }
    reverseCellMap_.append(cellI);
    cellZone_.append(zoneID);

    return cellI;
}


void Foam::directPolyTopoChange::modifyCell
(
    const label cellI,
    const label zoneID
)
{
    cellZone_[cellI] = zoneID;
}


void Foam::directPolyTopoChange::removeCell(const label cellI)
{
    if (cellI < 0 || cellI >= cellMap_.size())
    {
        FatalErrorIn("directPolyTopoChange::removeCell(const label)")
            << "illegal cell label " << cellI << endl
            << "Valid cell labels are 0 .. " << cellMap_.size()-1
            << abort(FatalError);
    }

    if (strict_ && cellMap_[cellI] == -2)
    {
        FatalErrorIn("directPolyTopoChange::removeCell(const label)")
            << "cell " << cellI
            << " already marked for removal"
            << abort(FatalError);
    }

    cellMap_[cellI] = -2;
    cellFromPoint_.erase(cellI);
    cellFromEdge_.erase(cellI);
    cellFromFace_.erase(cellI);
    cellZone_[cellI] = -1;
}


Foam::autoPtr<Foam::mapPolyMesh>
 Foam::directPolyTopoChange::changeMesh(polyMesh& mesh)
{
    if (debug)
    {
        Pout<< "directPolyTopoChange::changeMesh" << endl;
    }


    if (mesh.boundaryMesh().size() != nPatches_)
    {
        FatalErrorIn("directPolyTopoChange::changeMesh(polyMesh&)")
            << "directPolyTopoChange was constructed with a mesh with "
            << nPatches_ << " patches." << endl
            << "The mesh now provided has a different number of patches "
            << mesh.boundaryMesh().size()
            << " which is illegal" << endl
            << abort(FatalError);
    }

    // Remove any holes from points/faces/cells and sort faces.
    // Sets nActiveFaces_.
    compact();


    // Transfer points to pointField. points_ are now cleared!
    // Only done since e.g. reorderCoupledFaces requires pointField.
    pointField newPoints;
    newPoints.transfer(points_);


    // Get patch sizes
    labelList patchSizes;
    labelList patchStarts;
    calcPatchSizes(patchSizes, patchStarts);

    // Reorder any coupled faces
    reorderCoupledFaces
    (
        mesh.boundaryMesh(),
        patchStarts,
        patchSizes,
        newPoints
    );



    // Calculate face inflation maps
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are for the new face the old faces whose value needs to be
    // averaged to get the new value. These old faces come either from
    // the old edgeFaces (inflate from edge) or the old pointFaces (inflate
    // from point). As an additional complexity will use only internal faces
    // to create new value for internal face and vice versa only patch
    // faces to to create patch face value.

    List<objectMap> facesFromPoints(faceFromPoint_.size());
    List<objectMap> facesFromEdges(faceFromEdge_.size());

    calcFaceInflationMaps(mesh, facesFromPoints, facesFromEdges);


    // Calculate cell inflation maps
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<objectMap> cellsFromPoints(cellFromPoint_.size());
    List<objectMap> cellsFromEdges(cellFromEdge_.size());
    List<objectMap> cellsFromFaces(cellFromFace_.size());

    calcCellInflationMaps
    (
        mesh,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces
    );

    // Clear inflation info
    {
        faceFromPoint_.clear();
        faceFromPoint_.resize(0);
        faceFromEdge_.clear();
        faceFromEdge_.resize(0);

        cellFromPoint_.clear();
        cellFromPoint_.resize(0);
        cellFromEdge_.clear();
        cellFromEdge_.resize(0);
        cellFromFace_.clear();
        cellFromFace_.resize(0);
    }

    // Remove demand driven storage
    mesh.clearOut();



    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    // Grab patch mesh point maps
    List<Map<label> > oldPatchMeshPointMaps(boundary.size());
    labelList oldPatchNMeshPoints(boundary.size());
    labelList oldPatchStarts(boundary.size());

    forAll(boundary, patchI)
    {
        // Copy old face zone mesh point maps
        oldPatchMeshPointMaps[patchI] = boundary[patchI].meshPointMap();
        oldPatchNMeshPoints[patchI] = boundary[patchI].meshPoints().size();
        oldPatchStarts[patchI] = boundary[patchI].start();
    }

    // Grab old face zone mesh point maps.
    // These need to be saved before resetting the mesh and are used
    // later on to calculate the faceZone pointMaps.
    List<Map<label> > oldFaceZoneMeshPointMaps(mesh.faceZones().size());

    forAll(mesh.faceZones(), zoneI)
    {
        const faceZone& oldZone = mesh.faceZones()[zoneI];

        oldFaceZoneMeshPointMaps[zoneI] = oldZone().meshPointMap();
    }


    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());


    // Change the mesh
    // ~~~~~~~~~~~~~~~
    // This will invalidate any addressing so better make sure you have
    // all the information you need!!!

    mesh.resetPrimitives
    (
        nActiveFaces_,
        newPoints,
        faces_,
        faceOwner_,
        faceNeighbour_,
        patchSizes,
        patchStarts
    );

    // Clear out primitives
    {
        newPoints.clear();
        retiredPoints_.clear();
        retiredPoints_.resize(0);

        faces_.clear();
        faces_.setSize(0);
        region_.clear();
        region_.setSize(0);
        faceOwner_.clear();
        faceOwner_.setSize(0);
        faceNeighbour_.clear();
        faceNeighbour_.setSize(0);        
    }


    // Zones
    // ~~~~~

    // Inverse of point/face/cell zone addressing. 
    // For every preserved point/face/cells in zone give the old position.
    // For added points, the index is set to -1
    labelListList pointZoneMap(mesh.pointZones().size());
    labelListList faceZoneFaceMap(mesh.faceZones().size());
    labelListList cellZoneMap(mesh.cellZones().size());

    resetZones(mesh, pointZoneMap, faceZoneFaceMap, cellZoneMap);

    // Clear zone info
    {
        pointZone_.clear();
        pointZone_.resize(0);

        faceZone_.clear();
        faceZone_.resize(0);

        faceZoneFlip_.clear();
        faceZoneFlip_.resize(0);

        cellZone_.clear();
        cellZone_.setSize(0);
    }

    // Patch point renumbering
    // For every preserved point on a patch give the old position.
    // For added points, the index is set to -1
    labelListList patchPointMap(boundary.size());
    calcPatchPointMap
    (
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        boundary,
        patchPointMap
    );

    // Create the face zone mesh point renumbering
    labelListList faceZonePointMap(mesh.faceZones().size());
    calcFaceZonePointMap(mesh, oldFaceZoneMeshPointMaps, faceZonePointMap);


    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            mesh,
            nOldPoints,
            nOldFaces,
            nOldCells,

            pointMap_,
            faceMap_,
            facesFromPoints,
            facesFromEdges,

            cellMap_,
            cellsFromPoints,
            cellsFromEdges,
            cellsFromFaces,

            reversePointMap_,
            reverseFaceMap_,
            reverseCellMap_,

            flipFaceFlux_,

            patchPointMap,

            pointZoneMap,

            faceZonePointMap,
            faceZoneFaceMap,
            cellZoneMap,

            pointField(0),          // so no inflation
            oldPatchStarts,
            oldPatchNMeshPoints
        )
    );
}


// ************************************************************************* //
