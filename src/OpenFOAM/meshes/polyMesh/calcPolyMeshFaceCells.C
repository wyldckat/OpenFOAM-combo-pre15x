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

#include "polyMesh.H"
#include "primitiveMesh.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMesh::calcFaceCells() const
{
    // Loop through all faces and mark up cells. The first label in the faceCell
    // list is the lower.

    if (debug)
    {
        Info<< "void polyMesh::calcFaceCells() : "
            << "calculating faceCells" << endl;
    }

    // It is an error to attempt to recalculate faceCells
    // if the pointer is already set
    if (allOwnerPtr_ || allNeighbourPtr_)
    {
        FatalErrorIn("polyMesh::calcFaceCells() const")
            << "owner and neighbour already calculated"
            << abort(FatalError);
    }
    else
    {
        const pointField& p = allPoints();
        const faceList& f = allFaces();
        const cellList& c = allCells();

        allOwnerPtr_ = new labelList(f.size(), -1);
        allNeighbourPtr_ = new labelList(f.size(), -1);
        labelList& own = *allOwnerPtr_;
        labelList& nei = *allNeighbourPtr_;

        boolList markedFaces(f.size(), false);

        label nInternalFaces = 0;

        forAll (c, cellI)
        {
            // get reference to face labels for current cell
            const labelList& cellfaces = c[cellI];

            forAll (cellfaces, faceI)
            {
                if (!markedFaces[cellfaces[faceI]])
                {
                    // First visit: owner
                    own[cellfaces[faceI]] = cellI;

                    markedFaces[cellfaces[faceI]] = true;
                }
                else
                {
                    // Second visit: neighbour
                    nei[cellfaces[faceI]] = cellI;
                    nInternalFaces++;
                }
            }
        }

        // Count the number of real faces.
        // Note: if there are unused faces in the mesh, they should be
        // clustered at the end of the list

        label nUsedPoints = p.size();
        label nUsedFaces = own.size();

        forAll (own, faceI)
        {
            if (own[faceI] < 0)
            {
                nUsedFaces = faceI;
                break;
            }
        }

        // If not all faces are being used, check that all unused
        // faces are at the back of the list and check the number of
        // used points.
        if (nUsedFaces < own.size())
        {
            if (debug || morphDebug)
            {
                Info<< "void polyMesh::calcFaceCells() : "
                    << "unused faces detected.  "
                    << "Number of used faces: " << nUsedFaces 
                    << ".  Total number of faces: " << own.size() << endl;
            }

            for (label i = nUsedFaces; i < own.size(); i++)
            {
                if (own[i] >= 0)
                {
                    FatalErrorIn
                    (
                        "void polyMesh::calcFaceCells() const"
                    )   << "Error in face ordering: mixed used and unused "
                        << "faces at the end of face list." << nl
                        << "Number of used faces: " << nUsedFaces
                        << "  and face " << i << " is owned by cell " << own[i]
                        << abort(FatalError);
                }
            }

            // Mark the points used by live faces.
            boolList usedPoints(p.size(), false);

            for (label faceI = 0; faceI < nUsedFaces; faceI++)
            {
                const face& curFace = f[faceI];

                forAll (curFace, pointI)
                {
                    usedPoints[curFace[pointI]] = true;
                }
            }

            forAll (usedPoints, pointI)
            {
                if (!usedPoints[pointI])
                {
                    nUsedPoints = pointI;
                    break;
                }
            }

            if (nUsedPoints < p.size())
            {
                if (debug || morphDebug)
                {
                    Info<< "void polyMesh::calcFaceCells() : unused points "
                        << "detected.  Number of used points: "
                        << nUsedPoints << ". Total number of points: "
                        << p.size() << endl;
                }

                for (label i = nUsedPoints; i < p.size(); i++)
                {
                    if (usedPoints[i])
                    {
                        FatalErrorIn
                        (
                            "void polyMesh::calcFaceCells() const"
                        )   << "Error in point ordering: mixed used and unused "
                            << "points at the end of point list." << nl
                            << "Number of used points: " << nUsedPoints
                            << "  and point " << i << " is used by a live face."
                            << abort(FatalError);
                    }
                }
            }
        }

        // Reset the primitiveMesh
        const_cast<polyMesh&>(*this).primitiveMesh::reset
        (
            nUsedPoints,
            nInternalFaces,
            nUsedFaces,
            c.size(),
            p,
            f,
            own,
            nei
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& polyMesh::allOwner() const
{
    if (!allOwnerPtr_)
    {
        calcFaceCells();
    }

    return *allOwnerPtr_;
}

const labelList& polyMesh::allNeighbour() const
{
    if (!allNeighbourPtr_)
    {
        calcFaceCells();
    }

    return *allNeighbourPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
