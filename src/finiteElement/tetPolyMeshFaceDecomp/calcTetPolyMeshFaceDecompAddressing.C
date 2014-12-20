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
    The insertion of edges in the edge list is dictated by the upper
    triangular ordering. The current decomposition of a polyhedral cell into
    tetrahedra requires insertion of face centres and cell centres. The points
    are ordered in the following way:
        1) points of the polyMesh (supporting polyhedral cells)
        2) face centres
        3) cell centres

    The algorithm for owner-neighbour insertion first adds all the points the
    owner point shares the edge with (only the ones with the higher label than
    the owner point), followed by the face centres of the pointFaces, followed
    by the pointCells. This is because the face and the the cell centres are
    guaranteed to have a higher index than the internal vertices.

    Note:
    It is assumed that the edges are constructed such that the start label
    is lower than the end label and that pointFaces and pointCells lists are
    ordered.


\*---------------------------------------------------------------------------*/

#include "error.H"

#include "tetPolyMeshFaceDecomp.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void tetPolyMeshFaceDecomp::calcAddressing() const
{
    if (debug)
    {
        Info<< "void tetPolyMesh::calcAddressing() const : "
            << "Calculating tetPolyMesh addressing" << endl;
    }

    if (lduPtr_)
    {
        FatalErrorIn
        (
            "void tetPolyMeshFaceDecomp::calcAddressing() const"
        )   << "addressing already calculated"
            << abort(FatalError);
    }

    // Create owner and neighbour addressing list.
    // At the same time fill in the owner start lookup list
    if (debug)
    {       
        Info<< "face offset: " << faceOffset_ << " cell offset: " << cellOffset_
            << " n edges: " << nEdges() << endl;
    }

    labelList owner(nEdges(), -1);

    labelList neighbour(nEdges(), -1);

    // Get reference to edges and pointEdges
    const edgeList& meshEdges = mesh_.edges();

    // Get references to pointFaces and pointCells
    const labelListList& pointFaces = mesh_.pointFaces();
    const labelListList& pointCells = mesh_.pointCells();

    // Loop through all points
    label nCreatedEdges = 0;
    label curOwner = 0;
    label edgeI = 0;

    // Loop through all points
    forAll (pointFaces, pointI)
    {
        // Add the point neighbours

        // The edge construction is such that the end label (neighbour)
        // is higher than the start label (owner)
        while
        (
            edgeI < meshEdges.size()
         && meshEdges[edgeI].start() == pointI
        )
        {
            owner[nCreatedEdges] = curOwner;
            neighbour[nCreatedEdges] = meshEdges[edgeI].end();
            nCreatedEdges++;

            edgeI++;
        }

        // Add the face neighbours

        // Get the sorted list of pointFaces
        const labelList& curPointFaces = pointFaces[pointI];

        forAll (curPointFaces, faceI)
        {
            // add as neighbour
            owner[nCreatedEdges] = curOwner;
            neighbour[nCreatedEdges] = faceOffset_ + curPointFaces[faceI];
            nCreatedEdges++;
        }

        // Add the cell neighbours

        // Get the list of sorted pointCells
        const labelList& curPointCells = pointCells[pointI];

        forAll (curPointCells, cellI)
        {
            // Add as neighbour
            owner[nCreatedEdges] = curOwner;
            neighbour[nCreatedEdges] = cellOffset_ + curPointCells[cellI];
            nCreatedEdges++;
        }

        // Increment the current owner node
        curOwner++;
    }

    // Loop through all internal faces and add owner and neighbour of the face
    const unallocLabelList& meshOwner = mesh_.faceOwner();

    const unallocLabelList& meshNeighbour = mesh_.faceNeighbour();

    forAll (meshOwner, faceI)
    {
        // Add owner cell centre
        owner[nCreatedEdges] = curOwner;

        neighbour[nCreatedEdges] = cellOffset_ + meshOwner[faceI];

        nCreatedEdges++;

        // Inelegant. Change.
        if (faceI < meshNeighbour.size())
        {
            // Add neighbour cell centre
            owner[nCreatedEdges] = curOwner;
            neighbour[nCreatedEdges] = cellOffset_ + meshNeighbour[faceI];
            nCreatedEdges++;
        }

        curOwner++;
    }

    if (nCreatedEdges != nEdges())
    {
        FatalErrorIn("void tetPolyMesh::calcAddressing() const")
            << "Problem with edge counting in lduAddressing."
            << abort(FatalError);
    }

    // Add dummy boundary addressing
    labelListList dummyPatchAddr(boundary().size());

    forAll (dummyPatchAddr, patchI)
    {
        dummyPatchAddr[patchI].setSize(0);
    }

    lduPtr_ =
        new lduAddressingStore
        (
            nPoints(),
            owner,
            neighbour,
            dummyPatchAddr
        );

    if (debug)
    {
        Info<< "void tetPolyMesh::calcAddressing() const : "
            << "Finished calculating tetPolyMesh addressing" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const lduAddressing& tetPolyMeshFaceDecomp::ldu() const
{
    if (!lduPtr_)
    {
        calcAddressing();
    }

    return *lduPtr_;
}


label tetPolyMeshFaceDecomp::maxNPointsForCell() const
{
    if (maxNPointsForCell_ < 0)
    {
        const faceList& meshFaces = mesh_.faces();
        const cellList& polyCells = mesh_.cells();

        forAll (polyCells, cellI)
        {
            maxNPointsForCell_ =
                max
                (
                    maxNPointsForCell_,
                    polyCells[cellI].labels(meshFaces).size()
                  + polyCells[cellI].size()
                  + 1
                );
        }
    }

    return maxNPointsForCell_;
}


// Fill buffer with addressing for the cell
label tetPolyMeshFaceDecomp::addressing
(
    const label cellID,
    labelList& localToGlobalBuffer,
    labelList& globalToLocalBuffer
) const
{
    const unallocFaceList& meshFaces = mesh_.faces();

    const labelList& cellFaces = mesh_.cells()[cellID];

    label nextLocal = 0;

    // First mark up the vertices
    forAll (cellFaces, faceI)
    {
        const face& curFace = meshFaces[cellFaces[faceI]];

        forAll (curFace, pointI)
        {
            // If the point has not been already inserted into the local
            // buffer, add it
            if (globalToLocalBuffer[curFace[pointI]] == -1)
            {
                localToGlobalBuffer[nextLocal] = curFace[pointI];
                globalToLocalBuffer[curFace[pointI]] = nextLocal;
                nextLocal++;
            }
        }
    }

    // Mark up face centres
    forAll (cellFaces, faceI)
    {
        const label curFaceIndex = cellFaces[faceI] + faceOffset();

        // Mark up the face
        if (globalToLocalBuffer[curFaceIndex] == -1)
        {
            localToGlobalBuffer[nextLocal] = curFaceIndex;
            globalToLocalBuffer[curFaceIndex] = nextLocal;
            nextLocal++;
        }
    }

    // Mark up the cell centre
    localToGlobalBuffer[nextLocal] = cellOffset() + cellID;
    globalToLocalBuffer[cellOffset() + cellID] = nextLocal;
    nextLocal++;

    // Return size of addressing
    return nextLocal;
}


// Clear global to local addressing
void tetPolyMeshFaceDecomp::clearAddressing
(
    const label cellID,
    const label nCellPoints,
    labelList& localToGlobalBuffer,
    labelList& globalToLocalBuffer
) const
{
    // Only clear the places that have been used.  The rest of the buffer
    // is already initiated to -1
    for (label localI = 0; localI < nCellPoints; localI++)
    {
        globalToLocalBuffer[localToGlobalBuffer[localI]] = -1;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
