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

#include "addPatchCellLayer.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "ListOps.H"
#include "polyAddFace.H"
#include "IndirectList.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addPatchCellLayer::addPatchCellLayer(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::addPatchCellLayer::setRefinement
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp,
    polyTopoChange& meshMod
) const
{
    if (pp.nPoints() != patchDisp.size())
    {
        FatalErrorIn
        (
            "addPatchCellLayer::setRefinement(const indirectPrimitivePatch&"
            ", const pointField&, polyTopoChange&)"
        )   << "Size of new points is not same as number of points used by"
            << " the face subset" << endl
            << "  patch.nPoints:" << pp.nPoints()
            << "  points:" << patchDisp.size()
            << abort(FatalError);
    }

    // From master point (in patch point label) to added point (in mesh point
    // label)
    labelList addedPoints(pp.nPoints(), -1);


    //
    // Create new points
    //

    forAll(patchDisp, patchPointI)
    {
        if (mag(patchDisp[patchPointI]) > SMALL)
        {
            label meshPointI = pp.meshPoints()[patchPointI];

            label addedVertI = 
                meshMod.setAction
                (
                    polyAddPoint
                    (
                        pp.localPoints()[patchPointI] + patchDisp[patchPointI], // point
                        meshPointI,             // master point
                        -1,                     // zone for point
                        true                    // supports a cell
                    )
                );

            addedPoints[patchPointI] = addedVertI;
        }
    }


    //
    // Add cells to all boundaryFaces
    //

    labelList addedCells(pp.size(), -1);

    forAll(pp, patchFaceI)
    {
        // Add a cell (layer) to faceI

        const face& f = pp.localFaces()[patchFaceI];

        // Do not add cell if none of the vertices of the face gets duplicated
        bool needCell = false;

        forAll(f, fp)
        {
            if (addedPoints[f[fp]] != -1)
            {
                needCell = true;

                break;       
            }
        }

        if (needCell)
        {
            label meshFaceI = pp.addressing()[patchFaceI];

            addedCells[patchFaceI] =
                meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,             // master point
                        -1,             // master edge
                        meshFaceI,      // master face
                        -1,             // master cell id
                        -1              // zone for cell
                    )
                );
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Precalculated patchID for each patch face
    labelList patchID(pp.size());

    forAll(pp, patchFaceI)
    {
        label meshFaceI = pp.addressing()[patchFaceI];

        patchID[patchFaceI] = patches.whichPatch(meshFaceI);
    }


    //
    // Create new faces on the outside. These should have outwards orientation
    // so just create them as (renumbered) copies of old face since these
    // already have outwards pointing normals.
    //

    forAll(pp.localFaces(), patchFaceI)
    {
        if (addedCells[patchFaceI] != -1)
        {
            // Get duplicated vertices on the patch face.
            // Since at least one unique (see test before) there will be
            // a valid unique face.
            const face& f = pp.localFaces()[patchFaceI];

            face newFace(f.size());

            forAll(f, fp)
            {
                if (addedPoints[f[fp]] == -1)
                {
                    newFace[fp] = pp.meshPoints()[f[fp]];
                }
                else
                {
                    newFace[fp] = addedPoints[f[fp]];
                }
            }

            label meshFaceI = pp.addressing()[patchFaceI];

            meshMod.setAction
            (
                polyAddFace
                (
                    newFace,                // face
                    addedCells[patchFaceI], // owner
                    -1,                     // neighbour
                    -1,                     // master point
                    -1,                     // master edge
                    meshFaceI,              // master face for addition
                    false,                  // flux flip
                    patchID[patchFaceI],    // patch for face
                    -1,                     // zone for face
                    false                   // face zone flip
                )
            );
        }
    }

    //
    // Modify old patch faces to be on the inside
    //
    forAll(pp, patchFaceI)
    {
        if (addedCells[patchFaceI] != -1)
        {
            label meshFaceI = pp.addressing()[patchFaceI];

            meshMod.setAction
            (
                polyModifyFace
                (
                    pp[patchFaceI],                 // modified face
                    meshFaceI,                      // label of face being modified
                    mesh_.faceOwner()[meshFaceI],   // owner
                    addedCells[patchFaceI],         // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    -1,                             // zone for face
                    false                           // face flip in zone
                )
            );
        }
    }

    //
    // Create 'side' faces, one per edge that is being extended.
    //

    forAll(pp.edges(), patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];

        // Get corresponding mesh edge label
        label v0 = pp.meshPoints()[e.start()];
        label v1 = pp.meshPoints()[e.end()];
        const edge meshE(v0, v1);

        label meshEdgeI = -1;

        const labelList& v0Edges = mesh_.pointEdges()[meshE.start()];

        forAll(v0Edges, v0EdgeI)
        {
            label edgeI = v0Edges[v0EdgeI];

            if (mesh_.edges()[edgeI] == meshE)
            {
                meshEdgeI = edgeI;

                break;
            }
        }

        if (meshEdgeI == -1)
        {
            FatalErrorIn
            (
                "addPatchCellLayer::setRefinement"
                "(const pointField&, polyTopoChange&"
            )   << "Problem: cannot find mesh edge using points "
                << v0 << " and " << v1
                << abort(FatalError);
        }



        // Count number of internal faces used by meshEdgeI
        // If this is 0 we should not really inflate face from it (since mapping
        // might not make sense). Instead inflate face from nothing.

        label masterEdgeI = -1;

        // Mesh faces using edge
        const labelList& meshFaces = mesh_.edgeFaces()[meshEdgeI];

        forAll(meshFaces, i)
        {
            if (mesh_.isInternalFace(meshFaces[i]))
            {
                // meshEdge uses internal faces so ok to inflate from it
                masterEdgeI = meshEdgeI;

                break;
            }
        }


        // Patch faces using edge
        const labelList& patchFaces = pp.edgeFaces()[patchEdgeI];

        // Designate 'owner' of edge
        label patchFaceI = patchFaces[0];

        // Find order of v0, v1 in patchFaces

        const face& f = pp[patchFaceI];

        bool reverseFace = false;

        forAll(f, fp)
        {
            if (f[fp] == v0)
            {
                label fp1 = (fp + 1) % f.size();

                if (f[fp1] != v1)
                {
                    reverseFace = true;
                }
                break;
            }
        }

        // Create face consistent with having patchFaceI as owner
        face newFace(4);
        label newFp = 0;

        if (reverseFace)
        {
            newFace[newFp++] = v0;
            if (addedPoints[e.start()] != -1)
            {
                newFace[newFp++] = addedPoints[e.start()];
            }
            if (addedPoints[e.end()] != -1)
            {
                newFace[newFp++] = addedPoints[e.end()];
            }
            newFace[newFp++] = v1;
        }
        else
        {
            newFace[newFp++] = v0;
            newFace[newFp++] = v1;
            if (addedPoints[e.end()] != -1)
            {
                newFace[newFp++] = addedPoints[e.end()];
            }
            if (addedPoints[e.start()] != -1)
            {
                newFace[newFp++] = addedPoints[e.start()];
            }
        }

        if (newFp >= 3)
        {
            newFace.setSize(newFp);

            // Is edge external edge?
            if (patchFaces.size() == 1)
            {
                // External edge so external face. Patch id is obtained from
                // any other patch connected to edge.
                label meshFaceI = pp.addressing()[patchFaceI];

                // Loop over all faces connected to edge to inflate and
                // see if any boundary face (but not meshFaceI)
                label otherPatchID = patchID[patchFaceI];

                forAll(meshFaces, i)
                {
                    label faceI = meshFaces[i];

                    if (!mesh_.isInternalFace(faceI) && faceI != meshFaceI)
                    {
                        otherPatchID = patches.whichPatch(faceI);

                        break;
                    }
                }


                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                // face
                        addedCells[patchFaceI], // owner
                        -1,                     // neighbour
                        -1,                     // master point
                        masterEdgeI,            // master edge
                        -1,                     // master face for addition
                        false,                  // flux flip
                        otherPatchID,           // patch for face
                        -1,                     // zone for face
                        false                   // face zone flip
                    )
                );
            }
            else
            {
                label nbrFaceI = patchFaces[1];

                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,                // face
                        addedCells[patchFaceI], // owner
                        addedCells[nbrFaceI],   // neighbour
                        -1,                     // master point
                        masterEdgeI,            // master edge
                        -1,                     // master face for addition
                        false,                  // flux flip
                        -1,                     // patch for face
                        -1,                     // zone for face
                        false                   // face zone flip
                    )
                );
            }
        }
    }
}


// ************************************************************************* //
