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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "mapPolyMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::mapPolyMesh::mapPolyMesh
(
    const polyMesh& mesh,
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const labelList& pointMap,
    const labelList& faceMap,
    const List<objectMap>& facesFromPoints,
    const List<objectMap>& facesFromEdges,
    const labelList& cellMap,
    const List<objectMap>& cellsFromPoints,
    const List<objectMap>& cellsFromEdges,
    const List<objectMap>& cellsFromFaces,
    const labelList& reversePointMap,
    const labelList& reverseFaceMap,
    const labelList& reverseCellMap,
    const labelHashSet& flipFaceFlux,
    const labelListList& patchPointMap,
    const labelListList& pointZoneMap,
    const labelListList& faceZonePointMap,
    const labelListList& faceZoneFaceMap,
    const labelListList& cellZoneMap,
    const pointField& preMotionPoints,
    const labelList& oldPatchStarts,
    const labelList& oldPatchNMeshPoints
)
:
    mesh_(mesh),
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    pointMap_(pointMap),
    faceMap_(faceMap),
    facesFromPointsMap_(facesFromPoints),
    facesFromEdgesMap_(facesFromEdges),
    cellMap_(cellMap),
    cellsFromPointsMap_(cellsFromPoints),
    cellsFromEdgesMap_(cellsFromEdges),
    cellsFromFacesMap_(cellsFromFaces),
    reversePointMap_(reversePointMap),
    reverseFaceMap_(reverseFaceMap),
    reverseCellMap_(reverseCellMap),
    flipFaceFlux_(flipFaceFlux),
    patchPointMap_(patchPointMap),
    pointZoneMap_(pointZoneMap),
    faceZonePointMap_(faceZonePointMap),
    faceZoneFaceMap_(faceZoneFaceMap),
    cellZoneMap_(cellZoneMap),
    preMotionPoints_(preMotionPoints),
    oldPatchSizes_(oldPatchStarts.size()),
    oldPatchStarts_(oldPatchStarts),
    oldPatchNMeshPoints_(oldPatchNMeshPoints)
{
    // Calculate old patch sizes
    for (label patchI = 0; patchI < oldPatchStarts_.size() - 1; patchI++)
    {
        oldPatchSizes_[patchI] =
            oldPatchStarts_[patchI + 1] - oldPatchStarts_[patchI];
    }

    // Set the last one by hand
    const label lastPatchID = oldPatchStarts_.size() - 1;

    oldPatchSizes_[lastPatchID] = nOldFaces_ - oldPatchStarts_[lastPatchID];

    if (polyMesh::debug)
    {
        if (min(oldPatchSizes_) < 0)
        {
            FatalErrorIn("mapPolyMesh::mapPolyMesh(...)")
                << "Calculated negative old patch size.  Error in mapping data"
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapPolyMesh::reorderPatchFaces(const labelList& oldToNew)
{
    // faceMap_
    labelList newFaceMap(faceMap_);

    forAll (faceMap_, i)
    {
        newFaceMap[oldToNew[i]] = faceMap_[i];
    }

    faceMap_.transfer(newFaceMap);


    // facesFromPointsMap_
    forAll (facesFromPointsMap_, i)
    {
        const objectMap& obj = facesFromPointsMap_[i];

        facesFromPointsMap_[i] =
            objectMap
            (
                oldToNew[obj.index()],  // Renumbered face label
                obj.masterObjects()     // Original faces on old mesh
            );
    }

    // facesFromEdgesMap_
    forAll (facesFromEdgesMap_, i)
    {
        const objectMap& obj = facesFromEdgesMap_[i];

        facesFromEdgesMap_[i] =
            objectMap
            (
                oldToNew[obj.index()],  // Renumbered face label
                obj.masterObjects()     // Original faces on old mesh
            );
    }

    //- reverseFaceMap_
    //  (new face for every old face)
    forAll (reverseFaceMap_, faceI)
    {
        reverseFaceMap_[faceI] = oldToNew[reverseFaceMap_[faceI]];
    }
    
    //- flipFaceFlux_
    labelHashSet oldFFFlux(flipFaceFlux_);
    flipFaceFlux_.clear();

    for
    (
        labelHashSet::const_iterator iter = oldFFFlux.begin();
        iter != oldFFFlux.end();
        ++iter
    )
    {
        flipFaceFlux_.insert(oldToNew[iter.key()]);
    }

    // faceZoneFaceMap_
    labelListList newFaceZoneFaceMap(faceZoneFaceMap_);

    forAll (faceZoneFaceMap_, faceI)
    {
        newFaceZoneFaceMap[oldToNew[faceI]].transfer(faceZoneFaceMap_[faceI]);
    }

    faceZoneFaceMap_.transfer(newFaceZoneFaceMap);
}


// ************************************************************************* //
