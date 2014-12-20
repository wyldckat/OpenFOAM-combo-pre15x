/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "directCombineFaces.H"
#include "polyMesh.H"
#include "directPolyTopoChange.H"
#include "polyRemoveFace.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(directCombineFaces, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Walk face-edge-face to regionize faces.
void Foam::directCombineFaces::regionFaces
(
    const indirectPrimitivePatch& fp,
    const scalar featureCos,
    const label regionI,
    const label faceI,
    labelList& region
)
{
    if (region[faceI] == -1)
    {
        region[faceI] = regionI;

        const labelList& fEdges = fp.faceEdges()[faceI];

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            const labelList& eFaces = fp.edgeFaces()[edgeI];

            if (eFaces.size() == 2)
            {
                label otherI = eFaces[0];

                if (otherI == faceI)
                {
                    otherI = eFaces[1];
                }

                if (featureCos < (-1+SMALL))
                {
                    regionFaces(fp, featureCos, regionI, otherI, region);
                }
                else
                {
                    // They are both patch faces so normals pointing outwards.
                    scalar cosAng =
                        fp.faceNormals()[faceI]
                      & fp.faceNormals()[otherI];

                    if (cosAng > featureCos)
                    {
                        regionFaces(fp, featureCos, regionI, otherI, region);
                    }
                }
            }
        }
    }
}


// Test face for (almost) convexeness. Allows certain convexity before
// complaining.
// For every two consecutive edges calculate the normal. If it differs too
// much from the average face normal complain.
bool Foam::directCombineFaces::convexFace
(
    const scalar concaveSin,
    const pointField& points,
    const face& f
)
{
    // Get outwards pointing normal of f.
    vector n = f.normal(points);
    n /= mag(n);

    // Get edge from f[0] to f[size-1];
    vector ePrev(points[f[0]] - points[f[f.size()-1]]);
    scalar magEPrev = mag(ePrev);
    ePrev /= magEPrev + VSMALL;

    forAll(f, fp0)
    {
        // Get vertex after fp
        label fp1 = (fp0 + 1) % f.size();

        // Normalized vector between two consecutive points
        vector e10(points[f[fp1]] - points[f[fp0]]);
        scalar magE10 = mag(e10);
        e10 /= magE10 + VSMALL;

        if (magEPrev > SMALL && magE10 > SMALL)
        {
            vector edgeNormal = ePrev ^ e10;
            scalar magEdgeNormal = mag(edgeNormal);

            if (magEdgeNormal < concaveSin)
            {
                // Edges (almost) aligned -> face is ok.
            }
            else
            {
                // Check normal
                edgeNormal /= magEdgeNormal;

                if ((edgeNormal & n) < SMALL)
                {
                    return false;
                }
            }
        }

        ePrev = e10;
        magEPrev = magE10;
    }

    // Not a single internal angle is concave so face is convex.
    return true;
}


// Gets outside edgeloop as a face
// - in same order as faces
// - in mesh vertex labels
Foam::face Foam::directCombineFaces::getOutsideFace
(
    const indirectPrimitivePatch& fp
)
{
    if (fp.edgeLoops().size() != 1)
    {
        FatalErrorIn
        (
            "directCombineFaces::sameOrder(const indirectPrimitivePatch&)"
        )   << "Multiple outside loops:" << fp.edgeLoops()
            << abort(FatalError);
    }

    // Get first boundary edge. Since guaranteed one edgeLoop when in here this
    // edge must be on it.
    label bEdgeI = fp.nInternalEdges();

    const edge& e = fp.edges()[bEdgeI];

    const labelList& eFaces = fp.edgeFaces()[bEdgeI];

    if (eFaces.size() != 1)
    {
        FatalErrorIn
        (
            "directCombineFaces::sameOrder(const indirectPrimitivePatch&)"
        )   << "boundary edge:" << bEdgeI
            << " points:" << fp.meshPoints()[e[0]]
            << ' ' << fp.meshPoints()[e[1]]
            << " on indirectPrimitivePatch has " << eFaces.size()
            << " faces using it" << abort(FatalError);
    }


    // Outside loop
    const labelList& outsideLoop = fp.edgeLoops()[0];


    // Get order of edge e in outside loop
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // edgeLoopConsistent : true  = edge in same order as outsideloop
    //                      false = opposite order
    bool edgeLoopConsistent = false;

    {
        label index0 = findIndex(outsideLoop, e[0]);
        label index1 = findIndex(outsideLoop, e[1]);

        if (index0 == -1 || index1 == -1)
        {
            FatalErrorIn
            (
                "directCombineFaces::sameOrder(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in edgeLoop:" << outsideLoop << abort(FatalError);
        }
        else if (index1 == outsideLoop.fcIndex(index0))
        {
            edgeLoopConsistent = true;
        }
        else if (index0 == outsideLoop.fcIndex(index1))
        {
            edgeLoopConsistent = false;
        }
        else
        {
            FatalErrorIn
            (
                "directCombineFaces::sameOrder(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " on consecutive points in edgeLoop:"
                << outsideLoop << abort(FatalError);
        }
    }


    // Get face in local vertices
    const face& localF = fp.localFaces()[eFaces[0]];

    // Get order of edge in localF
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // faceEdgeConsistent : true  = edge in same order as localF
    //                      false = opposite order
    bool faceEdgeConsistent = false;

    {
        // Find edge in face.
        label index = findIndex(fp.faceEdges()[eFaces[0]], bEdgeI);

        if (index == -1)
        {
            FatalErrorIn
            (
                "directCombineFaces::sameOrder(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in face:" << eFaces[0]
                << " edges:" << fp.faceEdges()[eFaces[0]]
                << abort(FatalError);
        }
        else if (localF[index] == e[0] && localF.nextLabel(index) == e[1])
        {
            faceEdgeConsistent = true;
        }
        else if (localF[index] == e[1] && localF.nextLabel(index) == e[0])
        {
            faceEdgeConsistent = false;
        }
        else
        {
            FatalErrorIn
            (
                "directCombineFaces::sameOrder(const indirectPrimitivePatch&)"
            )   << "Cannot find boundary edge:" << e
                << " points:" << fp.meshPoints()[e[0]]
                << ' ' << fp.meshPoints()[e[1]]
                << " in face:" << eFaces[0] << " verts:" << localF
                << abort(FatalError);
        }
    }

    // Get face in mesh points.
    face meshFace(renumber(fp.meshPoints(), outsideLoop));

    if (faceEdgeConsistent != edgeLoopConsistent)
    {
        reverse(meshFace);
    }
    return meshFace;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::directCombineFaces::directCombineFaces(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::directCombineFaces::getFaceRegions
(
    const scalar featureCos,
    const scalar concaveSin,
    const labelList& setFaces
) const
{
    // Check
    if (debug)
    {
        label own = -1;

        forAll(setFaces, i)
        {
            label faceI = setFaces[i];

            if (mesh_.isInternalFace(setFaces[i]))
            {
                FatalErrorIn
                (
                    "directCombineFaces::getFaceRegions"
                    "(const scalar, const scalar, const labelListList&)"
                )   << "Can only merge boundary faces. face:" << faceI
                    << abort(FatalError);
            }

            if (own == -1)
            {
                own = mesh_.faceOwner()[faceI];
            }
            else if (mesh_.faceOwner()[faceI] != own)
            {
                FatalErrorIn
                (
                    "directCombineFaces::getFaceRegions"
                    "(const scalar, const scalar, const labelListList&)"
                )   << "Faces not on same cell. Owner1:" << own
                    << " owner2:" << mesh_.faceOwner()[faceI]
                    << " for face:" << faceI
                    << abort(FatalError);
            }
        }
    }


    // Addressing on setFaces
    indirectPrimitivePatch fp
    (
        IndirectList<face>
        (
            mesh_.faces(),
            setFaces
        ),
        mesh_.points()
    );

    // Divide faces into regions.
    labelList region(fp.size(), -1);

    // Current region
    label regionI = 0;

    // Current face (index in fp, not mesh)
    label startFaceI = 0;

    while (true)
    {
        // Find face with unset region
        for
        (
            ;
            startFaceI < region.size() && region[startFaceI] != -1;
            startFaceI++
        );

        if (startFaceI == region.size())
        {
            break;
        }

        // Mark all faces connected to startFaceI with current region number.
        regionFaces(fp, featureCos, regionI, startFaceI, region);

        regionI++;
    }

    // Collect into lists of face sets (in fp numbering)
    labelListList faceSets(invertOneToMany(regionI, region));
    // Renumber to mesh face labels
    forAll(faceSets, i)
    {
        inplaceRenumber(setFaces, faceSets[i]);
    }


    // Filter out non-convex faces.
    label convexI = 0;

    forAll(faceSets, setI)
    {
        labelList& setFaces = faceSets[setI];

        if (setFaces.size() == 1)
        {
            // Only one face in set so ok.
            if (convexI != setI)
            {
                faceSets[convexI++] = setFaces;
            }
        }   
        else
        {
            // Check whether faces combine into convex face

            // Make face out of setFaces
            indirectPrimitivePatch bigFace
            (
                IndirectList<face>
                (
                    mesh_.faces(),
                    setFaces
                ),
                mesh_.points()
            );

            // Get outside vertices (in local vertex numbering)
            const labelListList& edgeLoops = bigFace.edgeLoops();

            // Only store if -only one outside loop -which forms convex face
            if (edgeLoops.size() == 1)
            {
                face f(getOutsideFace(bigFace));

                //if (debug)
                //{
                //    Pout<< "For cell:" << mesh_.faceOwner()[setFaces[0]]
                //        << " trying to merge faces:" << setFaces << nl
                //        << "verts:"
                //        << IndirectList<face>(mesh_.faces(), setFaces)()
                //        << " outside loop:" << f
                //        << " points:"
                //        << IndirectList<point>(mesh_.points(), f)()
                //        << " normal:" << f.normal(mesh_.points())
                //        << endl;
                //}

                if (convexFace(concaveSin, mesh_.points(), f))
                {
                    if (convexI != setI)
                    {
                        faceSets[convexI].transfer(setFaces);
                    }
                    convexI++;
                }
            }
        }
    }
    faceSets.setSize(convexI);

    return faceSets;
}


Foam::labelListList Foam::directCombineFaces::getMergeSets
(
    const scalar featureCos,
    const scalar concaveSin
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Lists of faces that can be merged.
    DynamicList<labelList> allFaceSets;

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!patch.coupled())
        {
            // Create table from owner to mesh faces

            Map<labelList> ownerToFace(2*patch.size());

            forAll(patch, patchFaceI)
            {
                label meshFaceI = patch.start() + patchFaceI;

                label own = mesh_.faceOwner()[meshFaceI];

                Map<labelList>::iterator cellFind = ownerToFace.find(own);

                if (cellFind != ownerToFace.end())
                {
                    labelList& ownerFaces = cellFind();

                    ownerFaces.setSize(ownerFaces.size() + 1);
                    ownerFaces[ownerFaces.size()-1] = meshFaceI;
                }
                else
                {
                    ownerToFace.insert(own, labelList(1, meshFaceI));
                }
            }

            // For all cells with multiple faces on same patch find out
            // whether faces can be merged.

            for
            (
                Map<labelList>::const_iterator iter = ownerToFace.begin();
                iter != ownerToFace.end();
                ++iter
            )
            {
                const labelList& ownerFaces = iter();

                if (ownerFaces.size() > 1)
                {
                    labelListList faceSets
                    (
                        getFaceRegions
                        (
                            featureCos,
                            concaveSin,
                            ownerFaces
                        )
                    );

                    forAll(faceSets, i)
                    {
                        if (faceSets[i].size() > 1)
                        {
                            allFaceSets.append(faceSets[i]);
                        }
                    }
                }
            }
        }
    }
    return allFaceSets.shrink();
}


void Foam::directCombineFaces::setRefinement
(
    const labelListList& faceSets,
    directPolyTopoChange& meshMod
)
{
    boolList pointUsed(mesh_.nPoints(), true);

    forAll(faceSets, setI)
    {
        const labelList& setFaces = faceSets[setI];

        // Check
        forAll(setFaces, i)
        {
            if (mesh_.isInternalFace(setFaces[i]))
            {
                FatalErrorIn
                (
                    "directCombineFaces::setRefinement"
                    "(const labelListList&, directPolyTopoChange&)"
                )   << "Can only merge boundary faces but found internal face:"
                    << setFaces[i] << " in set " << setI
                    << abort(FatalError);
            }
        }


        // Make face out of setFaces
        indirectPrimitivePatch bigFace
        (
            IndirectList<face>
            (
                mesh_.faces(),
                setFaces
            ),
            mesh_.points()
        );

        // Get outside vertices (in local vertex numbering)
        const labelListList& edgeLoops = bigFace.edgeLoops();

        if (edgeLoops.size() != 1)
        {
            FatalErrorIn
            (
                "directCombineFaces::setRefinement"
                "(const labelListList&, directPolyTopoChange&)"
            )   << "Faces to-be-merged " << setFaces
                << " do not form a single big face." << nl
                << abort(FatalError);
        }


        // Delete all faces except master, whose face gets modified.

        // Modify master face
        // ~~~~~~~~~~~~~~~~~~

        label masterFaceI = setFaces[0];

        // Get outside face in mesh vertex labels
        face f(getOutsideFace(bigFace));

        label zoneID = mesh_.faceZones().whichZone(masterFaceI);

        bool zoneFlip = false;

        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(masterFaceI)];
        }

        label patchI = mesh_.boundaryMesh().whichPatch(masterFaceI);

        meshMod.setAction
        (
            polyModifyFace
            (
                f,                              // modified face
                masterFaceI,                    // label of face being modified
                mesh_.faceOwner()[masterFaceI], // owner
                -1,                             // neighbour
                false,                          // face flip
                patchI,                         // patch for face
                false,                          // remove from zone
                zoneID,                         // zone for face
                zoneFlip                        // face flip in zone
            )
        );


        // Delete all non-master faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label i = 1; i < setFaces.size(); i++)
        {
            meshMod.setAction(polyRemoveFace(setFaces[i]));
        }


        // Mark unused points for deletion
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const labelList& meshPoints = bigFace.meshPoints();

        labelHashSet outsidePoints(edgeLoops[0]);

        forAll(meshPoints, i)
        {
            if (!outsidePoints.found(i))
            {
                // Since all faces will be merged into one check if point
                // will still be needed.
                const labelList& pFaces = mesh_.pointFaces()[meshPoints[i]];

                if ((pFaces.size() - setFaces.size()) == 0)
                {
                    pointUsed[meshPoints[i]] = false;
                }
            }
        }
    }

    // If one or more processors use the point it needs to be kept.
    syncTools::syncPointList
    (
        mesh_,
        pointUsed,
        orEqOp<bool>(),
        false,              // null value
        false               // no separation
    );

    forAll(pointUsed, pointI)
    {
        if (!pointUsed[pointI])
        {
            meshMod.setAction(polyRemovePoint(pointI));
        }
    }
}


// ************************************************************************* //
