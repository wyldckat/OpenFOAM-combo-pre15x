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
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

\*---------------------------------------------------------------------------*/

#include "meshSubset.H"
#include "boolList.H"
#include "emptyPolyPatch.H"
#include "demandDrivenData.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshSubset::checkCellSubset() const
{
    if (!meshSubsetPtr_)
    {
        FatalErrorIn("bool meshSubset::checkCellSubset() const")
            << "Mesh subset not set.  Please set the cell map using "
            << "void setCellSubset(const labelHashSet& cellsToSubset)" << endl
            << "before attempting to access subset data"
            << abort(FatalError);

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::meshSubset::markPoints
(
    const labelList& curPoints,
    Map<label>& pointMap
)
{
    forAll (curPoints, pointI)
    {
        if (!pointMap.found(curPoints[pointI]))
        {
            pointMap.insert(curPoints[pointI], 0);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshSubset::meshSubset(const IOobject& io)
:
    fvMesh(io),
    meshSubsetPtr_(NULL),
    pointMap_(0),
    faceMap_(0),
    cellMap_(0),
    patchMap_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshSubset::~meshSubset()
{
    deleteDemandDrivenData(meshSubsetPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshSubset::setCellSubset
(
    const labelHashSet& globalCellMap,
    const label patchID
)
{
    // Initial check on patches before doing anything time consuming.
    label wantedPatchID = patchID;

    if (wantedPatchID == -1)
    {
        // No explicit patch specified. Put in oldInternalFaces patch.
        // Check if patch with this name already exists.
        wantedPatchID = boundaryMesh().findPatchID("oldInternalFaces");
    }
    else if (wantedPatchID < 0 || wantedPatchID >= boundary().size())
    {
        FatalErrorIn
        (
            "meshSubset::setCellSubset(const labelHashSet&"
            ", const label patchID)"
        )   << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << boundary().size()-1
            << abort(FatalError);
    }


    cellMap_ = (const labelList&)globalCellMap.toc();

    // Sort the cell map in the ascending order
    sort(cellMap_);

    // Approximate sizing parameters for face and point lists
    const label avgNFacesPerCell = 6;
    const label avgNPointsPerFace = 4;

    const cellList& oldCells = allCells();
    const faceList& oldFaces = allFaces();
    const pointField& oldPoints = allPoints();
    const labelList& oldAllOwner = allOwner();
    const labelList& oldNeighbour = neighbour();
    const polyBoundaryMesh& oldPatches = boundaryMesh();

    label nCellsInSet = cellMap_.size();

    // Mark all used faces

    Map<label> facesToSubset(avgNFacesPerCell*nCellsInSet);

    forAll (cellMap_, cellI)
    {
        // Mark all faces from the cell
        const labelList& curFaces = oldCells[cellMap_[cellI]];

        forAll (curFaces, faceI)
        {
            if (!facesToSubset.found(curFaces[faceI]))
            {
                facesToSubset.insert(curFaces[faceI], 1);
            }
            else
            {
                facesToSubset[curFaces[faceI]]++;
            }
        }
    }

    // Mark all used points and make a global-to-local face map
    Map<label> globalFaceMap(facesToSubset.size());

    // Make a global-to-local point map
    Map<label> globalPointMap(avgNPointsPerFace*facesToSubset.size());

    // This is done in two goes, so that the boundary faces are last
    // in the list.  Because of this, I need to create the face map
    // along the way rather than just grab the table of contents.
    labelList facesToc = (const labelList&)facesToSubset.toc();
    sort(facesToc);
    faceMap_.setSize(facesToc.size());

    // 1. Get all faces that will be internal to the submesh.
    forAll (facesToc, faceI)
    {
        if (facesToSubset[facesToc[faceI]] == 2)
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[faceI];
            globalFaceMap.insert(facesToc[faceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[faceI]], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    label nInternalFaces = globalFaceMap.size();


    // Where to insert old internal faces.
    label oldPatchStart = labelMax;
    if (wantedPatchID != -1)
    {
        oldPatchStart = oldPatches[wantedPatchID].start();
    }


    label faceI = 0;

    // 2. Boundary faces up to where we want to insert old internal faces
    for (; faceI< facesToc.size(); faceI++)
    {
        if (facesToc[faceI] >= oldPatchStart)
        {
            break;
        }
        if
        (
            !isInternalFace(facesToc[faceI])
         && facesToSubset[facesToc[faceI]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[faceI];
            globalFaceMap.insert(facesToc[faceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[faceI]], globalPointMap);
        }
    }

    // 3. old internal faces
    forAll(facesToc, intFaceI)
    {
        if
        (
            isInternalFace(facesToc[intFaceI])
         && facesToSubset[facesToc[intFaceI]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[intFaceI];
            globalFaceMap.insert(facesToc[intFaceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[intFaceI]], globalPointMap);
        }
    }

    // 4. Remaining boundary faces
    for (; faceI< facesToc.size(); faceI++)
    {
        if
        (
            !isInternalFace(facesToc[faceI])
         && facesToSubset[facesToc[faceI]] == 1
        )
        {
            // Mark face and increment number of points in set
            faceMap_[globalFaceMap.size()] = facesToc[faceI];
            globalFaceMap.insert(facesToc[faceI], globalFaceMap.size());

            // Mark all points from the face
            markPoints(oldFaces[facesToc[faceI]], globalPointMap);
        }
    }



    // Grab the points map
    pointMap_ = (const labelList&)globalPointMap.toc();
    sort(pointMap_);

    forAll (pointMap_, pointI)
    {
        globalPointMap[pointMap_[pointI]] = pointI;
    }

    Info << "Number of cells in new mesh: " << nCellsInSet << endl;
    Info << "Number of faces in new mesh: " << globalFaceMap.size() << endl;
    Info << "Number of points in new mesh: " << globalPointMap.size() << endl;

    // Make a new mesh
    pointField newPoints(globalPointMap.size());

    label nNewPoints = 0;

    forAll (pointMap_, pointI)
    {
        newPoints[nNewPoints] = oldPoints[pointMap_[pointI]];
        nNewPoints++;
    }

    faceList newFaces(globalFaceMap.size());

    label nNewFaces = 0;

    // Make internal faces
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const face& oldF = oldFaces[faceMap_[faceI]];

        face newF(oldF.size());

        forAll (newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }

    // Make boundary faces

    label nbSize = boundary().size();
    label oldInternalPatchID  = -1;

    if (wantedPatchID == -1)
    {
        // Create 'oldInternalFaces' patch at the end
        // and put all exposed internal faces in there.
        oldInternalPatchID = nbSize;
        nbSize++;

    }
    else
    {
        oldInternalPatchID = wantedPatchID;
    }


    // Grad size and start of each patch on the fly.  Because of the
    // structure of the underlying mesh, the patches will appear in the
    // ascending order
    labelList boundaryPatchSizes(nbSize, 0);

    // Assign boundary faces. Visited in order of faceMap_.
    for (label faceI = nInternalFaces; faceI < faceMap_.size(); faceI++)
    {
        label oldFaceI = faceMap_[faceI];

        face oldF = oldFaces[oldFaceI];

        // Turn the faces as necessary to point outwards
        if (isInternalFace(oldFaceI))
        {
            // Internal face. Possibly turned the wrong way round
            if
            (
                !globalCellMap.found(oldAllOwner[oldFaceI])
             && globalCellMap.found(oldNeighbour[oldFaceI])
            )
            {
                oldF = oldFaces[oldFaceI].reverseFace();
            }

            // Update count for patch
            boundaryPatchSizes[oldInternalPatchID]++;
        }
        else
        {
            // Boundary face. Increment the appropriate patch
            label patchOfFace = boundaryMesh().whichPatch(oldFaceI);

            // Update count for patch
            boundaryPatchSizes[patchOfFace]++;
        }

        face newF(oldF.size());

        forAll (newF, i)
        {
            newF[i] = globalPointMap[oldF[i]];
        }

        newFaces[nNewFaces] = newF;
        nNewFaces++;
    }



    // Create cells
    cellList newCells(nCellsInSet);

    label nNewCells = 0;

    forAll (cellMap_, cellI)
    {
        const labelList& oldC = oldCells[cellMap_[cellI]];

        labelList newC(oldC.size());

        forAll (newC, i)
        {
            newC[i] = globalFaceMap[oldC[i]];
        }

        newCells[nNewCells] = cell(newC);
        nNewCells++;
    }


    // Make a new mesh
    meshSubsetPtr_ = new fvMesh
    (
        IOobject
        (
            name() + "SubSet",
            time().timeName(),
            time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        newPoints,
        newFaces,
        newCells
    );


    // Add old patches
    List<polyPatch*> newBoundary(nbSize);
    patchMap_.setSize(nbSize);
    label nNewPatches = 0;
    label patchStart = nInternalFaces;

    forAll(oldPatches, patchI)
    {
        if (boundaryPatchSizes[patchI] > 0)
        {
            // Patch still exists. Add it
            newBoundary[nNewPatches] = oldPatches[patchI].clone
            (
                meshSubsetPtr_->boundaryMesh(),
                nNewPatches,
                boundaryPatchSizes[patchI],
                patchStart
            ).ptr();

            patchStart += boundaryPatchSizes[patchI];
            patchMap_[nNewPatches] = patchI;
            nNewPatches++;
        }
    }

    if (wantedPatchID == -1)
    {
        // Newly created patch so is at end. Check if any faces in it.
        if (boundaryPatchSizes[oldInternalPatchID] > 0)
        {
            newBoundary[nNewPatches] = new emptyPolyPatch
            (
                "oldInternalFaces",
                boundaryPatchSizes[oldInternalPatchID],
                patchStart,
                nNewPatches,
                meshSubsetPtr_->boundaryMesh()
            );

            // The index for the first patch is -1 as it originates from
            // the internal faces
            patchMap_[nNewPatches] = -1;
            nNewPatches++;
        }
    }

    // Reset the patch lists
    newBoundary.setSize(nNewPatches);
    patchMap_.setSize(nNewPatches);

    // Add the fvPatches
    meshSubsetPtr_->addFvPatches(newBoundary);
}


const fvMesh& Foam::meshSubset::subMesh() const
{
    checkCellSubset();

    return *meshSubsetPtr_;
}


const labelList& Foam::meshSubset::pointMap() const
{
    checkCellSubset();

    return pointMap_;
}


const labelList& Foam::meshSubset::faceMap() const
{
    checkCellSubset();

    return faceMap_;
}


const labelList& Foam::meshSubset::cellMap() const
{
    checkCellSubset();

    return cellMap_;
}


const labelList& Foam::meshSubset::patchMap() const
{
    checkCellSubset();

    return patchMap_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
