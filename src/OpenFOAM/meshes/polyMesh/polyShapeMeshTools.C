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
    Helper functions for construction of polyMesh from cell and patch shapes

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "primitiveMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

labelListList polyMesh::cellShapePointCells(const cellShapeList& c) const
{
    List<DynamicList<label, primitiveMesh::cellsPerPoint_> > 
        pc(allPoints().size());

    // For each cell
    forAll(c, i)
    {
        // For each vertex
        const labelList& labels = c[i];

        forAll(labels, j)
        {
            // Set working point label
            label curPoint = labels[j];
            DynamicList<label, primitiveMesh::cellsPerPoint_>& curPointCells =
                pc[curPoint];

            // Enter the cell label in the point's cell list
            curPointCells.append(i);
        }
    }

    labelListList pointCellAddr(pc.size());

    forAll (pc, pointI)
    {
        pointCellAddr[pointI].transfer(pc[pointI].shrink());
    }

    return pointCellAddr;
}


labelList polyMesh::facePatchFaceCells
(
    const faceList& patchFaces,
    const labelListList& pointCells,
    const faceListList& cellsFaceShapes,
    const label patchID
) const
{
    register bool found;

    labelList FaceCells(patchFaces.size());

    forAll(patchFaces, fI)
    {
        found = false;

        const face& curFace = patchFaces[fI];
        const labelList& facePoints = patchFaces[fI];

        forAll(facePoints, pointI)
        {
            const labelList& facePointCells = pointCells[facePoints[pointI]];

            forAll(facePointCells, cellI)
            {
                faceList cellFaces = cellsFaceShapes[facePointCells[cellI]];

                forAll(cellFaces, cellFace)
                {
                    if (cellFaces[cellFace] == curFace)
                    {
                        // Found the cell corresponding to this face
                        FaceCells[fI] = facePointCells[cellI];

                        found = true;
                    }
                    if (found) break;
                }
                if (found) break;
            }
            if (found) break;
        }

        if (!found)
        {
            FatalErrorIn
            (
                "polyMesh::facePatchFaceCells(const faceList& patchFaces,"
                "const labelListList& pointCells,"
                "const faceListList& cellsFaceShapes,"
                "const label patchID)"
            )   << "face " << fI << " in patch " << patchID
                << " does not have neighbour cell"
                << " face: " << patchFaces[fI]
                << abort(FatalError);
        }
    }

    return FaceCells;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
