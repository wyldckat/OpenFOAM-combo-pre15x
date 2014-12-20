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

Class
    mapAddedPolyMesh

\*----------------------------------------------------------------------------*/

#include "mapAddedPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::mapAddedPolyMesh::calcReverseMap
(
    const labelList& map,
    const bool getOldMesh,
    const label nOldElems
)
{
    labelList reverseMap(nOldElems, -1);

    if (getOldMesh)
    {
        // Get all negative entries
        forAll(map, i)
        {
            label index = map[i];

            if (index < 0)
            {
                label oldI = -index-1;

                reverseMap[oldI] = i;
            }
        }
    }
    else
    {
        // Get all positive entries
        forAll(map, i)
        {
            label index = map[i];

            if (index > 0)
            {
                label oldI = index-1;

                reverseMap[oldI] = i;
            }
        }
    }
    return reverseMap;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::mapAddedPolyMesh::mapAddedPolyMesh()
:
    nOldPoints_(0),
    nOldFaces_(0),
    nOldCells_(0),
    nAddedPoints_(0),
    nAddedFaces_(0),
    nAddedCells_(0),
    pointMap_(0),
    faceMap_(0),
    cellMap_(0),
    oldPatchMap_(0),
    addedPatchMap_(0)
{}


// Construct from components
Foam::mapAddedPolyMesh::mapAddedPolyMesh
(
    const label nOldPoints,
    const label nOldFaces,
    const label nOldCells,
    const label nAddedPoints,
    const label nAddedFaces,
    const label nAddedCells,
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap,
    const labelList& oldPatchMap,
    const labelList& addedPatchMap
)
:
    nOldPoints_(nOldPoints),
    nOldFaces_(nOldFaces),
    nOldCells_(nOldCells),
    nAddedPoints_(nAddedPoints),
    nAddedFaces_(nAddedFaces),
    nAddedCells_(nAddedCells),
    pointMap_(pointMap),
    faceMap_(faceMap),
    cellMap_(cellMap),
    oldPatchMap_(oldPatchMap),
    addedPatchMap_(addedPatchMap)
{}


// Construct as copy
Foam::mapAddedPolyMesh::mapAddedPolyMesh(const mapAddedPolyMesh& map)
:
    nOldPoints_(map.nOldPoints()),
    nOldFaces_(map.nOldFaces()),
    nOldCells_(map.nOldCells()),
    nAddedPoints_(map.nAddedPoints()),
    nAddedFaces_(map.nAddedFaces()),
    nAddedCells_(map.nAddedCells()),
    pointMap_(map.pointMap()),
    faceMap_(map.faceMap()),
    cellMap_(map.cellMap()),
    oldPatchMap_(map.oldPatchMap()),
    addedPatchMap_(map.addedPatchMap())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::mapAddedPolyMesh::combine
(
    const label nNew,
    const labelList& oldToNew,
    const labelList& addedToNew
)
{
    labelList map(nNew, -1);

    forAll(oldToNew, i)
    {
        if (oldToNew[i] != -1)
        {
            map[oldToNew[i]] = -i-1;
        }
    }

    forAll(addedToNew, i)
    {
        if (addedToNew[i] != -1)
        {
            map[addedToNew[i]] = i+1;
        }
    }

    return map;
}


Foam::labelList Foam::mapAddedPolyMesh::combine
(
    const label nNew,
    const label nOld,
    const labelList& addedToNew
)
{
    labelList map(nNew, -1);

    for (label i = 0; i < nOld; i++)
    {
        map[i] = -i-1;
    }

    forAll(addedToNew, i)
    {
        if (addedToNew[i] != -1)
        {
            map[addedToNew[i]] = i+1;
        }
    }

    return map;
}


Foam::labelList Foam::mapAddedPolyMesh::oldCellMap() const
{
    return calcReverseMap(cellMap_, true, nOldCells_);
}

//- From added mesh cells to new mesh cells. -1 if cell not added.
Foam::labelList Foam::mapAddedPolyMesh::addedCellMap() const
{
    return calcReverseMap(cellMap_, false, nAddedCells_);
}


//- From old mesh faces to new mesh faces.
Foam::labelList Foam::mapAddedPolyMesh::oldFaceMap() const
{
    return calcReverseMap(faceMap_, true, nOldFaces_);
}

//- From added mesh faces to new mesh facess. -1 if face not added.
Foam::labelList Foam::mapAddedPolyMesh::addedFaceMap() const
{
    return calcReverseMap(faceMap_, false, nAddedFaces_);
}


//- From old mesh points to new mesh points.
Foam::labelList Foam::mapAddedPolyMesh::oldPointMap() const
{
    return calcReverseMap(pointMap_, true, nOldPoints_);
}

//- From added mesh points to new mesh points. -1 if pt not added.
Foam::labelList Foam::mapAddedPolyMesh::addedPointMap() const
{
    return calcReverseMap(pointMap_, false, nAddedPoints_);
}


// ************************************************************************* //
