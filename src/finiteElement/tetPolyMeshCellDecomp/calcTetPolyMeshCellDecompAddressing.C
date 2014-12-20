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
    tetrahedra requires insertion of cell centres. The points
    are ordered in the following way:
        1) points of the polyMesh (supporting polyhedral cells)
        2) cell centres

    The algorithm for owner-neighbour insertion first adds all the
    points the owner point shares the edge with (only the ones with
    the higher label than the owner point), followed by the
    pointCells. This is because the the cell centres are
    guaranteed to have a higher index than the internal vertices.

    Note:
    It is assumed that the edges are constructed such that the start
    label is lower than the end label and that the pointCells list are
    ordered.

\*---------------------------------------------------------------------------*/

#include "tetPolyMeshCellDecomp.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void tetPolyMeshCellDecomp::calcAddressing() const
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
            "void tetPolyMeshCellDecomp::calcAddressing() const"
        )   << "lduPtr_ already allocated"
            << abort(FatalError);
    }

    // Create owner and neighbour addressing list.
    // At the same time fill in the owner start lookup list

    const faceList& faces = mesh_.faces();

    const labelListList& pointFaces = mesh_.pointFaces();
    const labelListList& pointCells = mesh_.pointCells();

    labelList owner(nEdges(), -1);
    labelList neighbour(nEdges(), -1);

    // Count the added owners and neighbours
    label nCreatedEdges = 0;

    label pointInFace, prev, next, f0;

    // Loop through all points
    forAll (pointFaces, pointI)
    {
        const labelList& curFaces = pointFaces[pointI];

        // Create a list of labels to keep the neighbours that
        // have already been added.  Size is estimated
        labelHashSet addedNeighbours
        (
            2*curFaces.size()*primitiveMesh::pointsPerFace_
        );

        forAll (curFaces, faceI)
        {
            const face& f = faces[curFaces[faceI]];

            // Grab zeroth label of face
            f0 = f[0];

            labelList neighbourPointsFromFace(f.size() - 1, -1);

            if (f0 == pointI)
            {
                // If the current point is the zero point of the face,
                // it is connected to all other points
                for (label nbrI = 1; nbrI < f.size(); nbrI++)
                {
                    if (f[nbrI] > pointI)
                    {
                        addedNeighbours.insert(f[nbrI]);
                    }
                }
            }
            else if (f[1] == pointI)
            {
                // If the current point is the first point of the face,
                // it is connected to all other points
                // if it is the last point, it is connected to point zero
                // and the penultimate point
                if (f0 > pointI)
                {
                    addedNeighbours.insert(f0);
                }

                if (f[2] > pointI)
                {
                    addedNeighbours.insert(f[2]);
                }
            }
            else if (f[f.size() - 1] == pointI)
            {
                // If it is the last point, it is connected to point zero
                // and the penultimate point
                if (f[f.size() - 2] > pointI)
                {
                    addedNeighbours.insert(f[f.size() - 2]);
                }

                if (f0 > pointI)
                {
                    addedNeighbours.insert(f0);
                }
            }
            else
            {
                // Otherwise, it is connected to the previous and the next
                // point and additionally to point zero
                pointInFace = f.which(pointI);
                prev = f.prevLabel(pointInFace);
                next = f.nextLabel(pointInFace);

                if (prev > pointI)
                {
                    addedNeighbours.insert(prev);
                }

                if (next > pointI)
                {
                    addedNeighbours.insert(next);
                }

                if (f0 > pointI)
                {
                    addedNeighbours.insert(f0);
                }
            }
        }

        // All neighbours for the current point found. Before adding
        // them to the list, it is necessary to sort them in the
        // increasing order of the neighbouring point.

        // Make real list out of SLList to simplify the manipulation.
        labelList an(addedNeighbours.toc());

        // Use a simple sort to sort the an list.
        sort(an);

        // Adding the neighbours
        forAll (an, edgeI)
        {
            owner[nCreatedEdges] = pointI;
            neighbour[nCreatedEdges] = an[edgeI];
            nCreatedEdges++;
        }

        // Now add cell neighbours
        const labelList& curPointCells = pointCells[pointI];

        forAll (curPointCells, cellI)
        {
            // add as neighbour
            owner[nCreatedEdges] = pointI;

            neighbour[nCreatedEdges] = cellOffset_ + curPointCells[cellI];

            nCreatedEdges++;
        }
    }

    // Add dummy boundary addressing
    labelListList dummyPatchAddr(boundary().size());
    forAll (dummyPatchAddr, patchI)
    {
        dummyPatchAddr[patchI].setSize(0);
    }

    if (nCreatedEdges != nEdges())
    {
        FatalErrorIn("void tetPolyMesh::calcAddressing() const")
            << "Problem with edge counting in lduAddressing: "
            << "the cell decomposition is multiply connected or otherwise "
            << "invalid.  Please Use face decomposition instead."
            << abort(FatalError);
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

const lduAddressing& tetPolyMeshCellDecomp::ldu() const
{
    if (!lduPtr_)
    {
        calcAddressing();
    }

    return *lduPtr_;
}


label tetPolyMeshCellDecomp::maxNPointsForCell() const
{
    if (maxNPointsForCell_ < 0)
    {
        const cellList& polyCells = mesh_.cells();
        const faceList& faces = mesh_.faces();

        forAll (polyCells, cellI)
        {
            maxNPointsForCell_ =
                max
                (
                    maxNPointsForCell_,
                    polyCells[cellI].labels(faces).size() + 1
                );
        }
    }

    return maxNPointsForCell_;
}


// Fill buffer with addressing for the cell
label tetPolyMeshCellDecomp::addressing
(
    const label cellID,
    labelList& localToGlobalBuffer,
    labelList& globalToLocalBuffer
) const
{
    const unallocFaceList& meshFaces = mesh_.faces();

    const labelList& cellFaces = mesh_.cells()[cellID];

    label nextLocal = 0;

    forAll (cellFaces, faceI)
    {
        const face& curFace = meshFaces[cellFaces[faceI]];

        // First mark up the vertices
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

    // Mark up the cell centre
    localToGlobalBuffer[nextLocal] = cellOffset() + cellID;
    globalToLocalBuffer[cellOffset() + cellID] = nextLocal;
    nextLocal++;

    // Return size of addressing
    return nextLocal;
}


// Clear global to local addressing
void tetPolyMeshCellDecomp::clearAddressing
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
