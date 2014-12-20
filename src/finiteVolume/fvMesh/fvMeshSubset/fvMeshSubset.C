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
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

    MJ 23/03/05 on coupled faces change the patch of the face to the
    oldInternalFaces patch.

\*---------------------------------------------------------------------------*/

#include "fvMeshSubset.H"
#include "boolList.H"
#include "emptyPolyPatch.H"
#include "demandDrivenData.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fvMeshSubset::checkCellSubset() const
{
    if (!fvMeshSubsetPtr_)
    {
        FatalErrorIn("bool fvMeshSubset::checkCellSubset() const")
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


void Foam::fvMeshSubset::markPoints
(
    const labelList& curPoints,
    Map<label>& pointMap
)
{
    forAll (curPoints, pointI)
    {
        // Note: insert will only insert if not yet there.
        pointMap.insert(curPoints[pointI], 0);
    }
}


void Foam::fvMeshSubset::markPoints
(
    const labelList& curPoints,
    labelList& pointMap
)
{
    forAll (curPoints, pointI)
    {
        pointMap[curPoints[pointI]] = 0;
    }
}


// Synchronize nCellsUsingFace on both sides of coupled patches. Marks
// faces that become 'uncoupled' with 3.
void Foam::fvMeshSubset::doCoupledPatches(labelList& nCellsUsingFace) const
{
    const polyBoundaryMesh& oldPatches = boundaryMesh();

    label nUncoupled = 0;

    if (Pstream::parRun())
    {
        // Send face usage across processor patches
        forAll (oldPatches, oldPatchI)
        {
            const polyPatch& pp = oldPatches[oldPatchI];

            if (typeid(pp) == typeid(processorPolyPatch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                OPstream toNeighbour(procPatch.neighbProcNo());

                toNeighbour
                    << SubList<label>(nCellsUsingFace, pp.size(), pp.start());
            }
        }


        // Receive face usage count and check for faces that become uncoupled.
        forAll (oldPatches, oldPatchI)
        {
            const polyPatch& pp = oldPatches[oldPatchI];

            if (typeid(pp) == typeid(processorPolyPatch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                IPstream fromNeighbour(procPatch.neighbProcNo());

                labelList nbrCellsUsingFace(fromNeighbour);

                // Combine with this side.

                forAll (pp, i)
                {
                    if
                    (
                        nCellsUsingFace[pp.start()+i] == 1
                     && nbrCellsUsingFace[i] == 0
                    )
                    {
                        // Face's neighbour is no longer there. Mark face off
                        // as coupled
                        nCellsUsingFace[pp.start()+i] = 3;
                        nUncoupled++;
                    }
                }       
            }
        }
    }
 
    // Do same for cyclics.
    forAll (oldPatches, oldPatchI)
    {
        const polyPatch& pp = oldPatches[oldPatchI];

        if (typeid(pp) == typeid(cyclicPolyPatch))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            forAll (cycPatch, i)
            {
                label thisFaceI = cycPatch.start() + i;
                label otherFaceI = cycPatch.transformGlobalFace(thisFaceI);

                if
                (
                    nCellsUsingFace[thisFaceI] == 1
                 && nCellsUsingFace[otherFaceI] == 0
                )
                {
                    nCellsUsingFace[thisFaceI] = 3;
                    nUncoupled++;
                }
            }
        }
    }

    reduce(nUncoupled, sumOp<label>());

    if (nUncoupled > 0)
    {
        Info<< "Uncoupled " << nUncoupled << " faces on coupled patches. "
            << "(processorPolyPatch, cyclicPolyPatch)" << nl
            << "You might need to run couplePatches to restore the patch face"
            << " ordering." << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fvMeshSubset::fvMeshSubset(const IOobject& io)
:
    fvMesh(io),
    fvMeshSubsetPtr_(NULL),
    pointMap_(0),
    faceMap_(0),
    cellMap_(0),
    patchMap_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshSubset::~fvMeshSubset()
{
    deleteDemandDrivenData(fvMeshSubsetPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshSubset::setCellSubset
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
            "fvMeshSubset::setCellSubset(const labelHashSet&"
            ", const label patchID)"
        )   << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << boundary().size()-1
            << abort(FatalError);
    }


    cellMap_ = globalCellMap.toc();

    // Sort the cell map in the ascending order
    sort(cellMap_);

    // Approximate sizing parameters for face and point lists
    const label avgNFacesPerCell = 6;
    const label avgNPointsPerFace = 4;

    const cellList& oldCells = cells();
    const faceList& oldFaces = faces();
    const pointField& oldPoints = points();
    const labelList& oldOwner = owner();
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
    labelList facesToc = facesToSubset.toc();
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
    pointMap_ = globalPointMap.toc();
    sort(pointMap_);

    forAll (pointMap_, pointI)
    {
        globalPointMap[pointMap_[pointI]] = pointI;
    }

    Pout << "Number of cells in new mesh: " << nCellsInSet << endl;
    Pout << "Number of faces in new mesh: " << globalFaceMap.size() << endl;
    Pout << "Number of points in new mesh: " << globalPointMap.size() << endl;

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
                !globalCellMap.found(oldOwner[oldFaceI])
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
    fvMeshSubsetPtr_ = new fvMesh
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
                fvMeshSubsetPtr_->boundaryMesh(),
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
                fvMeshSubsetPtr_->boundaryMesh()
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
    fvMeshSubsetPtr_->addFvPatches(newBoundary);
}


void Foam::fvMeshSubset::setLargeCellSubset
(
    const labelList& region,
    const label currentRegion,
    const label patchID
)
{
    const cellList& oldCells = cells();
    const faceList& oldFaces = faces();
    const pointField& oldPoints = points();
    const labelList& oldOwner = owner();
    const labelList& oldNeighbour = neighbour();
    const polyBoundaryMesh& oldPatches = boundaryMesh();
    const label oldNInternalFaces = nInternalFaces();


    // Initial checks

    if (region.size() != oldCells.size())
    {
        FatalErrorIn
        (
            "fvMeshSubset::setCellSubset(const labelList&"
            ", const label, const label)"
        )   << "Size of region " << region.size()
            << " is not equal to number of cells in mesh " << oldCells.size()
            << abort(FatalError);
    }


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
            "fvMeshSubset::setCellSubset(const labelList&"
            ", const label, const label)"
        )   << "Non-existing patch index " << wantedPatchID << endl
            << "Should be between 0 and " << boundary().size()-1
            << abort(FatalError);
    }


    // Get the cells for the current region.
    cellMap_.setSize(oldCells.size());
    label nCellsInSet = 0;

    forAll (region, oldCellI)
    {
        if (region[oldCellI] == currentRegion)
        {
            cellMap_[nCellsInSet++] = oldCellI;
        }
    }
    cellMap_.setSize(nCellsInSet);


    // Mark all used faces. Count number of cells using them
    // 0: face not used anymore
    // 1: face used by one cell, face becomes/stays boundary face
    // 2: face still used and remains internal face
    // 3: face coupled and used by one cell only (so should become normal,
    //    non-coupled patch face)
    //
    // Note that this is not really nessecary - but means we can size things
    // correctly. Also makes handling coupled faces much easier.

    labelList nCellsUsingFace(oldFaces.size(), 0);

    label nFacesInSet = 0;
    forAll (oldFaces, oldFaceI)
    {
        bool faceUsed = false;

        if (region[oldOwner[oldFaceI]] == currentRegion)
        {
            nCellsUsingFace[oldFaceI]++;
            faceUsed = true;
        }

        if
        (
            isInternalFace(oldFaceI)
         && (region[oldNeighbour[oldFaceI]] == currentRegion)
        )
        {
            nCellsUsingFace[oldFaceI]++;
            faceUsed = true;
        }

        if (faceUsed)
        {
            nFacesInSet++;
        }
    }
    faceMap_.setSize(nFacesInSet);

    // Handle coupled faces. Modifies patch faces to be uncoupled to 3.
    doCoupledPatches(nCellsUsingFace);


    // See which patch to use for exposed internal faces.
    label nbSize = boundary().size();
    label oldInternalPatchID = -1;

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
    labelList boundaryPatchSizes(nbSize, 0);


    // Make a global-to-local point map
    labelList globalPointMap(oldPoints.size(), -1);

    labelList globalFaceMap(oldFaces.size(), -1);
    label faceI = 0;

    // 1. Pick up all preserved internal faces.
    for (label oldFaceI = 0; oldFaceI < oldNInternalFaces; oldFaceI++)
    {
        if (nCellsUsingFace[oldFaceI] == 2)
        {
            globalFaceMap[oldFaceI] = faceI;
            faceMap_[faceI++] = oldFaceI;

            // Mark all points from the face
            markPoints(oldFaces[oldFaceI], globalPointMap);
        }
    }

    // These are all the internal faces in the mesh.
    label nInternalFaces = faceI;

    // 2. Boundary faces up to where we want to insert old internal faces
    for
    (
        label oldPatchI = 0;
        oldPatchI < oldPatches.size()
     && oldPatchI <= oldInternalPatchID;
        oldPatchI++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchI];

        label oldFaceI = oldPatch.start();

        forAll (oldPatch, i)
        {
            if (nCellsUsingFace[oldFaceI] == 1)
            {
                // Boundary face is kept.
                
                // Mark face and increment number of points in set
                globalFaceMap[oldFaceI] = faceI;
                faceMap_[faceI++] = oldFaceI;

                // Mark all points from the face
                markPoints(oldFaces[oldFaceI], globalPointMap);

                // Increment number of patch faces
                boundaryPatchSizes[oldPatchI]++;
            }
            oldFaceI++;
        }
    }

    // 3a. old internal faces that have become exposed.
    for (label oldFaceI = 0; oldFaceI < oldNInternalFaces; oldFaceI++)
    {
        if (nCellsUsingFace[oldFaceI] == 1)
        {
            globalFaceMap[oldFaceI] = faceI;
            faceMap_[faceI++] = oldFaceI;

            // Mark all points from the face
            markPoints(oldFaces[oldFaceI], globalPointMap);

            // Increment number of patch faces
            boundaryPatchSizes[oldInternalPatchID]++;
        }
    }

    // 3b. coupled patch faces that have become uncoupled.
    for
    (
        label oldFaceI = oldNInternalFaces;
        oldFaceI < oldFaces.size();
        oldFaceI++
    )
    {
        if (nCellsUsingFace[oldFaceI] == 3)
        {
            globalFaceMap[oldFaceI] = faceI;
            faceMap_[faceI++] = oldFaceI;

            // Mark all points from the face
            markPoints(oldFaces[oldFaceI], globalPointMap);

            // Increment number of patch faces
            boundaryPatchSizes[oldInternalPatchID]++;
        }
    }

    // 4. Remaining boundary faces
    for
    (
        label oldPatchI = oldInternalPatchID+1;
        oldPatchI < oldPatches.size();
        oldPatchI++
    )
    {
        const polyPatch& oldPatch = oldPatches[oldPatchI];

        label oldFaceI = oldPatch.start();

        forAll (oldPatch, i)
        {
            if (nCellsUsingFace[oldFaceI] == 1)
            {
                // Boundary face is kept.
                
                // Mark face and increment number of points in set
                globalFaceMap[oldFaceI] = faceI;
                faceMap_[faceI++] = oldFaceI;

                // Mark all points from the face
                markPoints(oldFaces[oldFaceI], globalPointMap);

                // Increment number of patch faces
                boundaryPatchSizes[oldPatchI]++;
            }
            oldFaceI++;
        }
    }

    if (faceI != nFacesInSet)
    {
        FatalErrorIn
        (
            "fvMeshSubset::setCellSubset(const labelList&"
            ", const label, const label)"
        )   << "Problem" << abort(FatalError);
    }


    // Grab the points map
    label nPointsInSet = 0;

    forAll(globalPointMap, pointI)
    {
        if (globalPointMap[pointI] != -1)
        {
            nPointsInSet++;
        }
    }
    pointMap_.setSize(nPointsInSet);

    nPointsInSet = 0;

    forAll(globalPointMap, pointI)
    {
        if (globalPointMap[pointI] != -1)
        {
            pointMap_[nPointsInSet] = pointI;
            globalPointMap[pointI] = nPointsInSet;
            nPointsInSet++;
        }
    }

    Pout<< "Number of cells in new mesh : " << cellMap_.size() << endl;
    Pout<< "Number of faces in new mesh : " << faceMap_.size() << endl;
    Pout<< "Number of points in new mesh: " << pointMap_.size() << endl;

    // Make a new mesh
    pointField newPoints(pointMap_.size());

    label nNewPoints = 0;

    forAll (pointMap_, pointI)
    {
        newPoints[nNewPoints] = oldPoints[pointMap_[pointI]];
        nNewPoints++;
    }

    faceList newFaces(faceMap_.size());

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


    // Make boundary faces. (different from internal since might need to be
    // flipped)
    for (label faceI = nInternalFaces; faceI < faceMap_.size(); faceI++)
    {
        label oldFaceI = faceMap_[faceI];

        face oldF = oldFaces[oldFaceI];

        // Turn the faces as necessary to point outwards
        if (isInternalFace(oldFaceI))
        {
            // Was internal face. Possibly turned the wrong way round
            if
            (
                region[oldOwner[oldFaceI]] != currentRegion
             && region[oldNeighbour[oldFaceI]] == currentRegion
            )
            {
                oldF = oldFaces[oldFaceI].reverseFace();
            }
        }

        // Relabel vertices of the (possibly turned) face.
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
    fvMeshSubsetPtr_ = new fvMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
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

    Pout<< "Number of faces in new patches:" << nl;

    forAll (oldPatches, patchI)
    {
        if (boundaryPatchSizes[patchI] > 0)
        {
            // Patch still exists. Add it
            newBoundary[nNewPatches] = oldPatches[patchI].clone
            (
                fvMeshSubsetPtr_->boundaryMesh(),
                nNewPatches,
                boundaryPatchSizes[patchI],
                patchStart
            ).ptr();

            Pout<< "    " << oldPatches[patchI].name() << " : "
                << boundaryPatchSizes[patchI] << endl;

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
                fvMeshSubsetPtr_->boundaryMesh()
            );

            Pout<< "    oldInternalFaces : "
                << boundaryPatchSizes[oldInternalPatchID] << endl;

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
    fvMeshSubsetPtr_->addFvPatches(newBoundary);
}


void Foam::fvMeshSubset::setLargeCellSubset
(
    const labelHashSet& globalCellMap,
    const label patchID
)
{
    labelList region(nCells(), 0);

    forAllConstIter (labelHashSet, globalCellMap, iter)
    {
        region[iter.key()] = 1;
    }
    setLargeCellSubset(region, 1, patchID);
}


const fvMesh& Foam::fvMeshSubset::subMesh() const
{
    checkCellSubset();

    return *fvMeshSubsetPtr_;
}


const labelList& Foam::fvMeshSubset::pointMap() const
{
    checkCellSubset();

    return pointMap_;
}


const labelList& Foam::fvMeshSubset::faceMap() const
{
    checkCellSubset();

    return faceMap_;
}


const labelList& Foam::fvMeshSubset::cellMap() const
{
    checkCellSubset();

    return cellMap_;
}


const labelList& Foam::fvMeshSubset::patchMap() const
{
    checkCellSubset();

    return patchMap_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
