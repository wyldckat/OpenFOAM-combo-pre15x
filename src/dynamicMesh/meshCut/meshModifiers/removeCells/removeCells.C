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

#include "removeCells.H"
#include "polyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(removeCells, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::removeCells::removeCells(const polyMesh& mesh, const word& patchName)
:
    mesh_(mesh),
    patchName_(patchName)
{
    label patchI = mesh.boundaryMesh().findPatchID(patchName_);

    if (patchI == -1)
    {
        FatalErrorIn("removeCells::removeCells(const polyMesh&, const word&)")
            << "Cannot find patch " << patchName << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::removeCells::setRefinement
(
    const labelList& cellLabels,
    polyTopoChange& meshMod
) const
{
    label patchI = mesh_.boundaryMesh().findPatchID(patchName_);

    boolList removedCell(mesh_.nCells(), false);

    // Go from labelList of cells-to-remove to a boolList and remove all
    // cells mentioned.
    forAll(cellLabels, i)
    {
        label cellI = cellLabels[i];

        removedCell[cellI] = true;

        meshMod.setAction(polyRemoveCell(cellI));
    }


    // Remove faces that are no longer used. Modify faces that
    // are used by one cell only.

    const faceList& faces = mesh_.faces();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        const face& f = faces[faceI];
        label own = faceOwner[faceI];
        label nei = faceNeighbour[faceI];

        if (removedCell[own])
        {
            if (removedCell[nei])
            {
                // Face no longer used
                meshMod.setAction(polyRemoveFace(faceI));
            }
            else
            {
                // nei is remaining cell. FaceI becomes external cell
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        f.reverseFace(),        // modified face
                        faceI,                  // label of face being modified
                        nei,                    // owner
                        -1,                     // neighbour
                        false,                  // face flip
                        patchI,                 // patch for face
                        false,                  // remove from zone
                        -1,                     // zone for face
                        false                   // face flip in zone
                    )
                );
            }
        }
        else if (removedCell[nei])
        {
            // own is remaining cell. FaceI becomes external cell.
            meshMod.setAction
            (
                polyModifyFace
                (
                    f,                      // modified face
                    faceI,                  // label of face being modified
                    own,                    // owner
                    -1,                     // neighbour
                    false,                  // face flip
                    patchI,                 // patch for face
                    false,                  // remove from zone
                    -1,                     // zone for face
                    false                   // face flip in zone
                )
            );
        }
    }

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        if (removedCell[faceOwner[faceI]])
        {
            meshMod.setAction(polyRemoveFace(faceI));
        }
    }        


    // Remove points that are no longer used.

    const labelListList& pointFaces = mesh_.pointFaces();

    forAll(pointFaces, pointI)
    {
        // Check all faces using point to see whether they have not been
        // deleted. Only if all faces have been deleted can this point be
        // deleted.
        const labelList& pFaces = pointFaces[pointI];

        bool pointUsed = false;

        forAll(pFaces, i)
        {
            label faceI = pFaces[i];

            if
            (
                !removedCell[faceOwner[faceI]]
             || (
                    mesh_.isInternalFace(faceI)
                 && !removedCell[faceNeighbour[faceI]]
                )
            )
            {
                pointUsed = true;
                break;
            }
        }


        if (!pointUsed)
        {
            meshMod.setAction(polyRemovePoint(pointI));
        }
    }
}


// ************************************************************************* //
