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
    Checks for neighbouring patch faces on same cell and combines them.

    Faces have to make angle less than the given feature angle to be combined.
    Edges with additional points on them are also recombined.


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "Map.H"
#include "ListOps.H"
#include "primitiveFacePatch.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Sin of angle between two consecutive edges on a face. If sin(angle) larger
// than this the face will be considered concave.
const scalar concaveSin = 30*mathematicalConstant::pi/180.0;


// Test face for (almost) convexeness. Allows certain convexity before
// complaining.
// For every two consecutive edges calculate the normal. If it differs too
// much from the average face normal complain.
bool convexFace(const pointField& points, const labelList& pointLabels)
{
    // Convert into face to be able to use face functionality
    face f(pointLabels);

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


// Walk face-edge-face to regionize faces.
void regionFaces
(
    const primitiveFacePatch& fp,
    const scalar minCos,
    const label currentRegion,
    const label faceI,
    labelList& region
)
{
    if (region[faceI] == -1)
    {
        region[faceI] = currentRegion;

        const labelList& fEdges = fp.faceEdges()[faceI];

        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            const labelList& eFaces = fp.edgeFaces()[edgeI];

            if (eFaces.size() == 2)
            {
                label otherFaceI = eFaces[0];

                if (otherFaceI == faceI)
                {
                    otherFaceI = eFaces[1];
                }

                // They are both patch faces so normals pointing outwards.
                scalar cosAng =
                    fp.faceNormals()[faceI]
                  & fp.faceNormals()[otherFaceI];

                if (cosAng > minCos)
                {
                    regionFaces(fp, minCos, currentRegion, otherFaceI, region);
                }
            }
        }
    }
}


// Cell cellI has multiple faces (meshFaces) on the same patchI. See if they
// can be merged.
label mergeCellFaces
(
    const polyMesh& mesh,
    const scalar minCos,
    const bool allowConcaveFaces,
    const label patchI,
    const label cellI,
    const labelList& meshFaces,
    polyTopoChange& meshMod
)
{
    label nChanged = 0;

    if (meshFaces.size() == 1)
    {
        // Quick rejection. Single face so leave as is.
        return nChanged;
    }

    faceList faces(meshFaces.size());

    forAll(meshFaces, i)
    {
        faces[i] = mesh.faces()[meshFaces[i]];
    }

    // Addressing on faces only in mesh vertices.
    primitiveFacePatch fp(faces, mesh.points());

    // Divide faces into regions.
    labelList region(fp.size(), -1);

    // Current region
    label currentRegion = 0;

    do
    {
        label startFaceI = findIndex(region, -1);

        if (startFaceI == -1)
        {
            break;
        }

        regionFaces(fp, minCos, currentRegion, startFaceI, region);

        currentRegion++;
    } while(true);


    //
    // Collect faces per region and replace by single face
    //

    for (label regionI = 0; regionI < currentRegion; regionI++)
    {
        faceList myFaces(fp.size());
        label myFaceI = 0;
        label masterFaceI = labelMax;

        // Linear search through region to collect all faces belonging to
        // current region. Is linear since usually currentRegion=1 and
        // only one face in fp.
        forAll(region, patchFaceI)
        {
            if (region[patchFaceI] == regionI)
            {
                // Copy face (global vertex numbering)
                myFaces[myFaceI++] = fp[patchFaceI];

                masterFaceI = min(masterFaceI, patchFaceI);
            }
        }

        if (myFaceI == 0)
        {
            FatalErrorIn("mergeCellFaces") << "Problem" << abort(FatalError);
        }
        else if (myFaceI == 1)
        {
            // Region consists of single face. Leave as is.
        }
        else
        {
            // Multiple faces of same region that can possibly be merged.
            myFaces.setSize(myFaceI);

            // Check if cell would get too few faces.
            // (current region (size myFaceI) would be replaced by single face
            //  so might cause cell < tetrahedron)
            if ((mesh.cells()[cellI].size() - (myFaceI-1)) < 4)
            {
                WarningIn
                (
                    "label mergeCellFaces(const polyMesh& mesh, "
                    "const scalar minCos, const bool allowConcaveFaces,"
                    "const label patchI, const label cellI, "
                    "const labelList& meshFaces, polyTopoChange& meshMod)"
                )   << "Not merging since cell would get less than 4 faces"
                    << nl
                    << "Patch:" << patchI
                    << " cell:" << cellI
                    << " region:" << regionI
                    << " masterFaceI:" << meshFaces[masterFaceI]
                    << " nFaces:" << mesh.cells()[cellI].size()
                    << " faces to merge:" << myFaceI
                    << endl;
            }

            // myFaces replacement
            primitiveFacePatch myFace(myFaces, mesh.points());
            
            // Get outside vertices (in local vertex numbering)
            const labelListList& edgeLoops = myFace.edgeLoops();

            if
            (
                edgeLoops.size() != 1
             || (
                    !allowConcaveFaces
                 && !convexFace(myFace.localPoints(), edgeLoops[0])
                )
            )
            {
                // Region is not convex or has hole. Don't merge.
                WarningIn
                (
                    "label mergeCellFaces(const polyMesh& mesh, "
                    "const scalar minCos, const bool allowConcaveFaces,"
                    "const label patchI, const label cellI, "
                    "const labelList& meshFaces, polyTopoChange& meshMod)"
                )   << "Not merging since faces of region do not form a single"
                    << " convex, hole-free face."
                    << nl
                    << "Patch:" << patchI
                    << " cell:" << cellI
                    << " region:" << regionI
                    << " masterFaceI:" << meshFaces[masterFaceI]
                    << " edgeLoops:" << edgeLoops.size()
                    << endl;

                return nChanged;
            }


            // Convert list of local vertex numbers into global face
            face f(edgeLoops[0].size());

            forAll(edgeLoops[0], i)
            {
                f[i] = myFace.meshPoints()[edgeLoops[0][i]];
            }



            // Delete all faces except master, whose face gets modified.

            forAll(region, patchFaceI)
            {
                if (region[patchFaceI] == regionI)
                {
                    // Get mesh face.
                    label faceI = meshFaces[patchFaceI];

                    if (patchFaceI == masterFaceI)
                    {
                        label zoneID = mesh.faceZones().whichZone(faceI);

                        bool zoneFlip = false;

                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = mesh.faceZones()[zoneID];

                            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                        }

                        meshMod.setAction
                        (
                            polyModifyFace
                            (
                                f,              // modified face
                                faceI,          // label of face being modified
                                cellI,          // owner
                                -1,             // neighbour
                                false,          // face flip
                                patchI,         // patch for face
                                false,          // remove from zone
                                zoneID,         // zone for face
                                zoneFlip        // face flip in zone
                            )
                        );

                        nChanged++;
                    }
                    else
                    {
                        meshMod.setAction(polyRemoveFace(faceI));

                        nChanged++;
                    }
                }
            }
        }
    }

    return nChanged;
}


// Remove points not used by any face.
label removeHangingPoints(polyMesh& mesh)
{
    // Topology changes container
    polyTopoChange meshMod(mesh);

    boolList pointUsed(mesh.nPoints(), false);

    forAll(mesh.faces(), faceI)
    {
        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            pointUsed[f[fp]] = true;
        }
    }

    label nRemoved = 0;
    forAll(pointUsed, pointI)
    {
        if (!pointUsed[pointI])
        {
            meshMod.setAction(polyRemovePoint(pointI));
            nRemoved++;
        }
    }

    if (nRemoved > 0)
    {
        Info<< "Morphing to remove unused points ..." << endl;

        autoPtr<mapPolyMesh> morphMap = polyTopoChanger::changeMesh
        (
            mesh,
            meshMod
        );

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }
    }
    else
    {
        Info<< "No unused points removed ..." << endl;
    }

    return nRemoved;
}


// Remove points on straight edges. Works on all edges, not just boundary ones.
label mergeEdges(polyMesh& mesh, const scalar minCos)
{
    const labelListList& pointEdges = mesh.pointEdges();
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();

    // Collect points to remove and faces affected by them
    boolList pointRemoved(mesh.nPoints(), false);
    // Running count
    label nRemoved = 0;
    // Faces affected by points removed
    labelHashSet facesAffected;

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        if (pEdges.size() == 2)
        {
            const edge& e0 = edges[pEdges[0]];
            const edge& e1 = edges[pEdges[1]];

            label v0 = e0.commonVertex(e1);
            label vLeft = e0.otherVertex(v0);
            label vRight = e1.otherVertex(v0);

            vector e0Vec = points[v0] - points[vLeft];
            e0Vec /= mag(e0Vec);

            vector e1Vec = points[vRight] - points[v0];
            e1Vec /= mag(e1Vec);

            if ((e0Vec & e1Vec) > minCos)
            {
                // v0 candidate for removal. Used in straight edge only.
                pointRemoved[v0] = true;
                nRemoved++;

                const labelList& pFaces = mesh.pointFaces()[v0];

                forAll(pFaces, i)
                {
                    facesAffected.insert(pFaces[i]);
                }
            }
        }
    }

    for
    (
        labelHashSet::const_iterator iter = facesAffected.begin();
        iter != facesAffected.end();
        ++iter
    )
    {
        label faceI = iter.key();

        const face& f = mesh.faces()[faceI];

        label nRemaining = f.size();

        // Count number of remaining face points.
        forAll(f, fp)
        {
            if (pointRemoved[f[fp]])
            {
                if (nRemaining == 3)
                {
                    //Info<< "Unremoving point " << f[fp] << " since face "
                    //    << faceI << " verts:" << f << " would get too little"
                    //    << "vertices" << endl;

                    pointRemoved[f[fp]] = false;
                    nRemoved--;
                }
                else
                {
                    nRemaining--;
                }
            }
        }
    }

    // Topology changes container
    polyTopoChange meshMod(mesh);

    // Actually remove points
    forAll(pointRemoved, pointI)
    {
        if (pointRemoved[pointI])
        {
            meshMod.setAction(polyRemovePoint(pointI));
        }
    }

    // Update faces
    for
    (
        labelHashSet::const_iterator iter = facesAffected.begin();
        iter != facesAffected.end();
        ++iter
    )
    {
        label faceI = iter.key();

        const face& f = mesh.faces()[faceI];

        face newFace(f.size());

        label newI = 0;

        forAll(f, fp)
        {
            if (!pointRemoved[f[fp]])
            {
                newFace[newI++] = f[fp];
            }
        }
        newFace.setSize(newI);

        // Get other face data.
        label patchI = -1;
        label owner = mesh.faceOwner()[faceI];
        label neighbour = -1;

        if (mesh.isInternalFace(faceI))
        {
            neighbour = mesh.faceNeighbour()[faceI];
        }
        else
        {
            patchI = mesh.boundaryMesh().whichPatch(faceI);
        }

        label zoneID = mesh.faceZones().whichZone(faceI);

        bool zoneFlip = false;

        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
        }

        meshMod.setAction
        (
            polyModifyFace
            (
                newFace,        // modified face
                faceI,          // label of face being modified
                owner,          // owner
                neighbour,      // neighbour
                false,          // face flip
                patchI,         // patch for face
                false,          // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }

    if (nRemoved > 0)
    {
        Info<< "Morphing to remove straight edge points ..." << endl;

        autoPtr<mapPolyMesh> morphMap = polyTopoChanger::changeMesh
        (
            mesh,
            meshMod
        );

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }
    }
    else
    {
        Info<< "No straight edges simplified ..." << endl;
    }

    return nRemoved;
}


label mergePatchFaces
(
    polyMesh& mesh,
    const scalar minCos,
    const bool allowConcaveFaces
)
{
    // Topology changes container
    polyTopoChange meshMod(mesh);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nChanged = 0;

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!patch.coupled())
        {
            Info<< "Trying to merge faces on patch " << patch.name() << endl;

            // Create table from owner to mesh faces

            Map<labelList> ownerToFace(2*patch.size());

            forAll(patch, patchFaceI)
            {
                label meshFaceI = patch.start() + patchFaceI;

                label own = mesh.faceOwner()[meshFaceI];

                Map<labelList>::iterator cellFind = ownerToFace.find(own);

                if (cellFind != ownerToFace.end())
                {
                    labelList& pFaces = cellFind();

                    pFaces.setSize(pFaces.size() + 1);
                    pFaces[pFaces.size()-1] = meshFaceI;
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
                label cellI = iter.key();

                const labelList& pFaces = iter();

                nChanged +=
                    mergeCellFaces
                    (
                        mesh,
                        minCos,
                        allowConcaveFaces,
                        patchI,
                        cellI,
                        pFaces,
                        meshMod
                    );
            }
        }
    }


    if (nChanged > 0)
    {
        // Do all changes
        Info<< "Morphing to merge faces ..." << endl;

        autoPtr<mapPolyMesh> morphMap = polyTopoChanger::changeMesh
        (
            mesh,
            meshMod
        );

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }
    }
    else
    {
        Info<< "No faces merged ..." << endl;
    }

    return nChanged;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("feature angle [0..180]");
    argList::validOptions.insert("allowConcaveFaces", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    scalar featureAngle(readScalar(IStringStream(args.args()[3])()));

    scalar minCos = Foam::cos(featureAngle * mathematicalConstant::pi/180.0);

    bool allowConcaveFaces = args.options().found("allowConcaveFaces");

    if (allowConcaveFaces)
    {
        Info<< "Merging all faces of a cell" << nl
            << "    - which are on the same patch" << nl
            << "    - which make an angle < " << featureAngle << nl
            << "      (cos:" << minCos << ')' << nl
            << "    - even when resulting face becomes concave" << nl
            << endl;
    }
    else
    {
        Info<< "Merging all faces of a cell" << nl
            << "    - which are on the same patch" << nl
            << "    - which make an angle < " << featureAngle << nl
            << "      (cos:" << minCos << ')' << nl
            << "    - except when resulting face would become concave" << nl
            << endl;
    }


    runTime++;

    // Merge faces on same patch
    label nChanged = mergePatchFaces(mesh, minCos, allowConcaveFaces);

    // Remove unused points
    nChanged += removeHangingPoints(mesh);

    // Merge points on straight edges.
    nChanged += mergeEdges(mesh, minCos);


    if (nChanged > 0)
    {
        Info<< "Writing morphed mesh to time " << runTime.timeName() << endl;

        mesh.write();
    }
    else
    {
        Info<< "Mesh unchanged." << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
