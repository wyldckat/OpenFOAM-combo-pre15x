/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    Implements patch face ordering on cyclic halves.

    initOrder: does nothing
    order:
        - finds unmodified pairs of faces which were opposite each other
        on the old patch. These are the seed faces to start walking from.
        - walk from the seed faces starting from f[0] and in anticlockwise
        (half0) or clockwise (half1) direction and store the face labels
        and the index in face.
        - use the face labels as facemap, use the index-in-face to rotate
        the faces on half1.

    So faces on both halves get reordered and the faces on the second half get
    rotated as well. This reordering on the first half is not really nessecary
    but it gives a nice almost upper-triangular order to the faces.

\*---------------------------------------------------------------------------*/

#include "cyclicPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "patchZones.H"
#include "walkPatch.H"
#include "OFstream.H"
#include "primitiveFacePatch.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar cyclicPolyPatch::featureCos_ = 0.9;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Get (pairs of) faces that are unchanged in the whole morphing process.
//  These are used to determine a relative ordering of connected faces later
//  on.
bool cyclicPolyPatch::getUnmodifiedFacePairs
(
    const polyTopoChange& ref,
    const mapPolyMesh& map,
    List<FixedList<label, 2> >& unmodFaces
) const
{
    const labelList& oldToNew = map.reverseFaceMap();
    const labelList& oldPatchStarts = map.oldPatchStarts();

    // Get modified/unmodified status for every face on old mesh.
    boolList unmodifiedOldFace(getUnmodifiedFaces(ref, map));

    // Find any unmodified old face which is on the current patch
    // so we know what this patch came from. (probably the same as it was but
    // you never know)
    label oldPatchI = -1;
    forAll(unmodifiedOldFace, oldFaceI)
    {
        if (unmodifiedOldFace[oldFaceI] && inPatch(oldToNew, oldFaceI))
        {
            oldPatchI = whichPatch(oldPatchStarts, oldFaceI);

            if (oldPatchI == -1)
            {
                FatalErrorIn("cyclicPolyPatch::getUnmodifiedFacePairs")
                    << "unmodified old face:" << oldFaceI
                    << " new face:" << oldToNew[oldFaceI]
                    << " is in new patch " << name()
                    << " but was not in any patch in old mesh." << endl
                    << "Old patch starting faces:" << oldPatchStarts
                    << abort(FatalError);
            }

            break;
        }
    }


    if (oldPatchI == -1)
    {
        SeriousErrorIn("cyclicPolyPatch::getUnmodifiedFacePairs")
            << " patch:" << name() << " : "
            << "Cannot find unmodified face on patch " << name() << endl
            << "Please use geometric ordering or adapt your topology changes"
            << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }

    const label oldStart = oldPatchStarts[oldPatchI];

    // Get size of old patch

    label oldSize;

    if (oldPatchI == oldPatchStarts.size()-1)
    {
        oldSize = map.nOldFaces() - oldStart;
    }
    else
    {
        oldSize = oldPatchStarts[oldPatchI+1] - oldStart;
    }

    if (polyMesh::morphDebug)
    {
        Info<< "cyclicPolyPatch::getUnmodifiedFacePairs : "
            << "Old patch start:" << oldStart << " size:" << oldSize
            << "  New patch start:" << start() << " size:" << size() << endl;
    }

    // Collect remaining faces which are the unmodified faces.
    // Store corresponding pairs of these (in new mesh labels!)

    label oldHalf = oldSize/2;

    unmodFaces.setSize(oldHalf);
    label unmodI = 0;

    for (label oldFaceA = oldStart; oldFaceA < oldStart+oldHalf; oldFaceA++)
    {
        label oldFaceB = oldFaceA + oldHalf;

        label faceA = oldToNew[oldFaceA];
        label faceB = oldToNew[oldFaceB];

        if
        (
            unmodifiedOldFace[oldFaceA]
         && unmodifiedOldFace[oldFaceB]
         && faceA != -1
         && faceB != -1
        )
        {
            unmodFaces[unmodI][0] = faceA - start();
            unmodFaces[unmodI][1] = faceB - start();
            unmodI++;
        }
    }
    unmodFaces.setSize(unmodI);

    if (unmodFaces.size() == 0)
    {
        SeriousErrorIn("cyclicPolyPatch::getUnmodifiedFacePairs")
            << " patch:" << name() << " : "
            << "Cannot find two corresponding faces (one on each side of"
            << " cyclic half) that are not affected"
            << " by the morph." << endl
            << "This is illegal when morphing cyclic patches; you need at least"
            << " one unchanged pair of faces per cyclic patch area" << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

            return false;
    }

    if (polyMesh::morphDebug)
    {
        Info<< "cyclicPolyPatch::getUnmodifiedFacePairs : Dumping " << unmodI
            << " unmodified face pairs to OBJ file unmodPairs.obj" << endl;

        OFstream unmodStr("unmodPairs.obj");

        // Construct face centres using premotion points.
        const pointField& pts = map.preMotionPoints();
        const primitivePatch& pp = *this;

        label vertI = 0;

        forAll(unmodFaces, pairI)
        {
            const FixedList<label, 2>& pair = unmodFaces[pairI];

            // Write edge between two face centres
            const point c0(pp[pair[0]].centre(pts));
            const point c1(pp[pair[1]].centre(pts));
            writeOBJ(unmodStr, c0, c1, vertI);
        }
    }
    return true;
}


bool cyclicPolyPatch::geometricOrder
(
    const polyTopoChange&,
    const mapPolyMesh& map,
    const patchZones& ppZones,
    const pointField& normals,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (ppZones.nZones() != 2)
    {
        SeriousErrorIn("cyclicPolyPatch::geometricOrder")
            << " patch:" << name() << " : "
            << "Patch gets decomposed in "
            << ppZones.nZones() << " zones" << endl
            << "This means that the patch is either not two separate regions"
            << " or one region where the angle between the different regions"
            << " is not sufficiently sharp." << endl
            << "Please use topological matching or adapt the featureCos()"
            << " setting" << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }

    const primitivePatch& pp = *this;

    // Get the two halves and construct PrimitivePatches on them

    // From half0 to patch index
    labelList half0ToPatch(size(), -1);

    // All faces of half0
    faceList half0Faces(size());
    label n0Faces = 0;
    label face0 = -1;       // Arbitrary face on half0 (used for normal)

    labelList half1ToPatch(size(), -1);
    faceList half1Faces(size());
    label n1Faces = 0;
    label face1 = -1;       // Arbitrary face on half1 (used for normal)

    forAll(ppZones, faceI)
    {
        if (ppZones[faceI] == 0)
        {
            face0 = faceI;
            half0ToPatch[n0Faces] = faceI;
            half0Faces[n0Faces++] = pp[faceI];
        }
        else
        {
            face1 = faceI;
            half1ToPatch[n1Faces] = faceI;
            half1Faces[n1Faces++] = pp[faceI];
        }
    }

    if (n0Faces != n1Faces)
    {
        SeriousErrorIn("cyclicPolyPatch::geometricOrder")
            << " patch:" << name() << " : "
            << "Patch " << name() << " gets decomposed in two zones of"
            << "inequal size: " << n0Faces << " and " << n1Faces << endl
            << "This means that the patch is either not two separate regions"
            << " or one region where the angle between the different regions"
            << " is not sufficiently sharp." << endl
            << "Please use topological matching or adapt the featureCos()"
            << " setting" << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }

    half0Faces.setSize(n0Faces);
    half0ToPatch.setSize(n0Faces);
    half1Faces.setSize(n1Faces);
    half1ToPatch.setSize(n1Faces);


    primitiveFacePatch half0(half0Faces, map.preMotionPoints());
    primitiveFacePatch half1(half1Faces, map.preMotionPoints());


    // Get geometric data on both halves.
    const vector& n0 = normals[face0];
    const vector& n1 = normals[face1];

    pointField half0Ctrs(calcFaceCentres(half0, half0.points()));
    pointField anchors0(getAnchorPoints(half0, half0.points()));
    pointField half1Ctrs(calcFaceCentres(half1, half0.points()));

    if (mag(n0 & n1) < 1-SMALL)
    {
        if (polyMesh::morphDebug)
        {
            Info<< "cyclicPolyPatch::geometricOrder : Rotation :"
                << " n0:" << n0 << " n1:" << n1 << endl;
        }

        // Rotation (around origin)
        const tensor reverseT(transformationTensor(n0, -n1));

        // Rotation
        forAll(half0Ctrs, faceI)
        {
            half0Ctrs[faceI] = Foam::transform(reverseT, half0Ctrs[faceI]);
            anchors0[faceI] = Foam::transform(reverseT, anchors0[faceI]);
        }
    }
    else
    {
        // Parallel translation
        const pointField& half0Pts = half0.localPoints();
        const point ctr0(sum(half0Pts)/half0Pts.size());

        const pointField& half1Pts = half1.localPoints();
        const point ctr1(sum(half1Pts)/half1Pts.size());

        if (polyMesh::morphDebug)
        {
            Info<< "cyclicPolyPatch::geometricOrder : Translation :"
                << " n0:" << n0 << " n1:" << n1
                << " ctr0:" << ctr0 << " ctr1:" << ctr1 << endl;
        }

        half0Ctrs += ctr1 - ctr0;
        anchors0 += ctr1 - ctr0;
    }

    // Calculate typical distance per face
    scalarField tols(calcFaceTol(half1, half1.points(), half1Ctrs));

    // Geometric match of face centre vectors
    labelList from1To0;

    bool matchedAll = matchPoints(half1Ctrs, half0Ctrs, tols, false, from1To0);

    if (polyMesh::morphDebug)
    {
        fileName ccName(name() + "_faceCentres.obj");

        Info<< "cyclicPolyPatch::order : "
            << "Dumping newly found cyclic match as lines between"
            << " corresponding face centres to file " << ccName
            << endl;

        pointField half0Ctrs(calcFaceCentres(half0, half0.points()));
        pointField half1Ctrs(calcFaceCentres(half1, half0.points()));

        OFstream ccStr(ccName);

        label vertI = 0;

        forAll(half1Ctrs, i)
        {
            const point& c1 = half1Ctrs[i];

            if (from1To0[i] != -1)
            {
                // Write edge between c1 and c0
                const point& c0 = half0Ctrs[from1To0[i]];

                writeOBJ(ccStr, c0, c1, vertI);
            }
        }
    }


    if (!matchedAll)
    {
        SeriousErrorIn("cyclicPolyPatch::geometricOrder")
            << " patch:" << name() << " : "
            << "Cannot match vectors to faces on both sides of patch" << endl
            << "half0Ctrs[0]:" << half0Ctrs[0] << endl
            << "half1Ctrs[0]:" << half1Ctrs[0] << endl
            << "Please use topological matching or adapt the featureCos()"
            << " setting" << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

            return false;
    }

    // Set faceMap such that half0 faces get first and corresponding half1
    // faces last.

    forAll(half0ToPatch, half0FaceI)
    {
        // Label in original patch
        label patchFaceI = half0ToPatch[half0FaceI];

        faceMap[patchFaceI] = half0FaceI;

        // No rotation
        rotation[patchFaceI] = 0;
    }

    forAll(from1To0, half1FaceI)
    {
        label patchFaceI = half1ToPatch[half1FaceI];

        // This face has to match the corresponding one on half0.
        label half0FaceI = from1To0[half1FaceI];

        label newFaceI = half0FaceI + size()/2;

        faceMap[patchFaceI] = newFaceI;

        // Rotate patchFaceI such that its f[0] aligns with that of
        // the corresponding face
        // (which after shuffling will be at position half0FaceI)

        const point& wantedAnchor = anchors0[half0FaceI];

        rotation[newFaceI] =
            getRotation
            (
                half1.points(),
                half1[half1FaceI],
                wantedAnchor,
                tols[half1FaceI]
            );

        if (rotation[newFaceI] == -1)
        {
            SeriousErrorIn("cyclicPolyPatch::geometricOrder")
                << " patch:" << name() << " : "
                << "Cannot find point on face " << half1[half1FaceI]
                << " with vertices:"
                << IndirectList<point>(half1.points(), half1[half1FaceI])
                << " that matches point " << wantedAnchor
                << " when matching the halves of cyclic patch " << name()
                << endl
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }
    }
    
    forAll(faceMap, faceI)
    {
        if (faceMap[faceI] != faceI || rotation[faceI] != 0)
        {
            return true;
        }
    }
    return false;
}


bool cyclicPolyPatch::topologicalOrder
(
    const polyTopoChange& ref,
    const mapPolyMesh& map,
    const patchZones& ppZones,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (ppZones.nZones() % 2 == 1)
    {
        SeriousErrorIn("cyclicPolyPatch::topologicalOrder")
            << " patch:" << name() << " : "
            << "Patch gets decomposed in an odd number"
            << " of zones:" << ppZones.nZones() << endl
            << "This means that a geometric test on angle of neigbouring"
            << " faces does not correctly split the patch." << endl
            << "Please correct your geometry or adapt the featureCos()"
            << " setting" << endl
            << "Continuing with incorrect face ordering from now on!" << endl;

        return false;
    }



    // Pairs of corresponding unmodified faces (in new mesh face labels)
    List<FixedList<label, 2> > unmodFaces;

    bool ok = getUnmodifiedFacePairs(ref, map, unmodFaces);

    if (!ok)
    {
        // No unmodified face pairs found.
        return false;
    }

    // Use the unmodFaces as seeds to determine a corresponding ordering on
    // both sides.

    // All visited faces
    boolList visited(size(), false);

    // Running count of visited.
    label nVisited = 0;

    // Has anything changed?
    bool changed = false;

    // New face index
    label newFaceI = 0;

    forAll(unmodFaces, seedI)
    {
        // Walk on patchA starting from seed faces to get relative connectivity

        const label startFaceA = unmodFaces[seedI][0];

        walkPatch halfAWalker
        (
            *this,
            ppZones,
            false,          // Normal face order
            startFaceA,
            localFaces()[startFaceA][0],
            visited
        );

        nVisited += halfAWalker.visitOrder().size();

        if (nVisited == size())
        {
            SeriousErrorIn("cyclicPolyPatch::topologicalOrder")
                << " patch:" << name() << " : "
                << "Starting from patch face " << startFaceA
                << " mesh face:" << startFaceA + start() << " ended up visiting"
                << " the whole patch." << endl
                << "This means that the cyclic patch is not topologically"
                << " separated in two separate halves" << endl
                << "Please adapt the featureCos setting or use geometric"
                << " matching" << endl
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }

        // Walk corresponding path starting from  startFaceB
        const label startFaceB = unmodFaces[seedI][1];

        walkPatch halfBWalker
        (
            *this,
            ppZones,
            true,          // Walk reverse face order
            startFaceB,
            localFaces()[startFaceB][0],
            visited
        );

        nVisited += halfBWalker.visitOrder().size();

        // Now the visiting order becomes the way in which the faces are to be
        // ordered.
        const DynamicList<label>& visitOrderA = halfAWalker.visitOrder();

        forAll(visitOrderA, i)
        {
            label faceA = visitOrderA[i];

            if (faceA != newFaceI)
            {
                changed = true;
            }

            faceMap[faceA] = newFaceI;
            rotation[newFaceI] = 0;     // Do not rotate face
            newFaceI++;
        }
        // Corresponding order on second half
        const DynamicList<label>& visitOrderB = halfBWalker.visitOrder();

        if (polyMesh::morphDebug)
        {
            //Info<< "cyclicPolyPatch::order : "
            //    << "Walked from unmodified face pair."
            //    << " From startFaceA:" << startFaceA
            //    << " visited " << visitOrderA.size() << " faces."
            //    << "  From startFaceB:" << startFaceB
            //    << " visited " << visitOrderB.size() << " faces." << endl;

            pointField ctrs(calcFaceCentres(*this, map.preMotionPoints()));

            fileName ccNameA(name() + "_A" + Foam::name(seedI) + ".obj");

            Info<< "cyclicPolyPatch::order : Dumping visit order for seed "
                << seedI << " on A to OBJ file " << ccNameA << endl;

            OFstream ccStrA(ccNameA);
            writeOBJ(ccStrA, ctrs, visitOrderA);

            fileName ccNameB(name() + "_B" + Foam::name(seedI) + ".obj");

            Info<< "cyclicPolyPatch::order : Dumping visit order for seed "
                << seedI << " on B to OBJ file " << ccNameB << endl;

            OFstream ccStrB(ccNameB);
            writeOBJ(ccStrB, ctrs, visitOrderB);
        }


        if (visitOrderA.size() != visitOrderB.size())
        {
            SeriousErrorIn("cyclicPolyPatch::topologicalOrder")
                << " patch:" << name() << " : "
                << "Walk on both cyclic halves did not produce same walk"
                << " From startFaceA : " << startFaceA
                << " visited " << visitOrderA.size() << " faces" << endl
                << " From startFaceB : " << startFaceB
                << " visited " << visitOrderB.size() << " faces" << endl
                << "Please adapt the featureCos setting or use geometric"
                << " matching" << endl
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }


        const DynamicList<label>& indexInFaceB = halfBWalker.indexInFace();
        const DynamicList<label>& indexInFaceA = halfAWalker.indexInFace();

        forAll(visitOrderB, i)
        {
            label faceB = visitOrderB[i];

            faceMap[faceB] = newFaceI;

            // Reorient face
            const face& f = localFaces()[faceB];

            label nShift =
                ((f.size()-indexInFaceA[i]) - indexInFaceB[i]) % f.size();

            if (faceB != newFaceI || nShift != 0)
            {
                changed = true;
            }

            rotation[newFaceI] = nShift;

            newFaceI++;
        }

        if (nVisited == size())
        {
            // Done all faces; no need to try any other starting pairs
            break;
        }
    }

    // Check if all faces visited.
    forAll(visited, patchFaceI)
    {
        if (!visited[patchFaceI])
        {
            SeriousErrorIn("cyclicPolyPatch::topologicalOrder")
                << " patch:" << name() << " : "
                << "Did not visit mesh face "
                << start() + patchFaceI
                << " on cyclic patch " << name()
                << " from any of the seed faces." << endl
                << "There probably is not an unmodified face on every cyclic"
                << " half." << endl
                << "Please adapt the featureCos setting or use geometric"
                << " matching" << endl
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }
    }

    return changed;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar cyclicPolyPatch::featureCos()
{
    return featureCos_;
}


//- Set the feature angle, returning the previous value
scalar cyclicPolyPatch::setFeatureCos(const scalar t)
{
    if (t < -VSMALL || t > 1 + VSMALL)
    {
        FatalErrorIn
        (
            "scalar cyclicPolyPatch::setFeatureCos(const scalar t)"
        )   << "Illegal value for feature cos:" << t
            << abort(FatalError);
    }

    scalar oldTol = featureCos_;
    featureCos_ = t;

    return oldTol;
}


//- Initialize ordering (on new mesh)
void cyclicPolyPatch::initOrder
(
    const polyTopoChange&,
    const mapPolyMesh&
) const
{}


//- Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool cyclicPolyPatch::order
(
    const polyTopoChange& ref,
    const mapPolyMesh& map,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(size());
    faceMap = -1;

    rotation.setSize(size());
    rotation = 0;

    if (size() == 0)
    {
        // No faces, nothing to change.
        return false;
    }

    // Get geometric zones of patch by looking at normals

    boolList regionEdge(nEdges(), false);

    // Added faces have no updated normal yet (points are vector::zero).
    // So use preMotionPoints
    // and hope for the best. Cyclics are limited to one transformation tensor
    // currently anyway (i.e. straight plane) so should not be too big a
    // problem.
    const pointField& preMotionPoints = map.preMotionPoints();

    vectorField normals(size());

    const primitivePatch& pp = *this;

    forAll(pp, faceI)
    {
        const face& f = pp[faceI];

        normals[faceI] = f.normal(preMotionPoints);
    }
    normals /= mag(normals) + VSMALL;

    const labelListList& eFaces = edgeFaces();

    label nRegionEdges = 0;

    forAll(eFaces, edgeI)
    {
        const labelList& locEFaces = eFaces[edgeI];

        if (locEFaces.size() == 2)
        {
            if ((normals[locEFaces[0]] & normals[locEFaces[1]]) < featureCos())
            {
                regionEdge[edgeI] = true;

                nRegionEdges++;
            }
        }
    }


    // For every face determine zone it is connected to (without crossing
    // any regionEdge)
    patchZones ppZones(*this, regionEdge);

    if (polyMesh::morphDebug)
    {
        Info<< "cyclicPolyPatch::order : "
            << "Found " << nRegionEdges << " on patch " << name()
            << " where the cos of the angle between two connected faces"
            << " was less than " << featureCos() << endl;

        PtrList<OFstream> streams(ppZones.nZones());

        forAll(streams, zoneI)
        {
            fileName zoneName(name() + "_zone_" + Foam::name(zoneI) + ".obj");

            Info<< "cyclicPolyPatch::order : Writing zone " << zoneI
                << " face centres to OBJ file " << zoneName << endl;

            streams.hook(new OFstream(zoneName));
        }

        labelList nZoneFaces(ppZones.nZones(), 0);

        forAll(ppZones, faceI)
        {
            label zoneI = ppZones[faceI];

            const point cc(operator[](faceI).centre(map.preMotionPoints()));

            writeOBJ(streams[zoneI], cc);

            nZoneFaces[zoneI]++;
        }

        Info<< "cyclicPolyPatch::order : Number of faces per zone:"
            << nZoneFaces << endl;
    }


    // Has anything changed?
    bool changed = false;

    if (geometricMatch())
    {
        changed = geometricOrder(ref, map, ppZones, normals, faceMap, rotation);
    }
    else
    {
        changed = topologicalOrder(ref, map, ppZones, faceMap, rotation);
    }

    return changed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
