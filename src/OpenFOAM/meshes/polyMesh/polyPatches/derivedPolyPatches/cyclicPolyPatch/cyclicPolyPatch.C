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

\*---------------------------------------------------------------------------*/

#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "patchZones.H"
#include "matchPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicPolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, word);
addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, Istream);
addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, dictionary);

scalar cyclicPolyPatch::featureCos_ = 0.9;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Force calculation of transformation tensors
void cyclicPolyPatch::calcTransforms()
{
    if (size() > 0)
    {
        const pointField& points = allPoints();

        const face& f0 = static_cast<const faceList&>(*this)[0];
        const face& fn2 = static_cast<const faceList&>(*this)[size()/2];

        vector nf0 = f0.normal(points);
        nf0 /= mag(nf0);

        vector nfn2 = fn2.normal(points);
        nfn2 /= mag(nfn2);

        calcTransformTensors
        (
            f0.centre(points),
            fn2.centre(points),
            nf0,
            nfn2
        );
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

// Construct from components
cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}


// Construct from Istream
cyclicPolyPatch::cyclicPolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(is, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}


// Construct from dictionary
cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}

//- Construct as copy, resetting the boundary mesh
cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}

//- Construct as copy, resetting the face list and boundary mesh data
cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    calcTransforms();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar cyclicPolyPatch::featureCos()
{
    return featureCos_;
}

void cyclicPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}

void cyclicPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}

void cyclicPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}

void cyclicPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    calcTransforms();
}

void cyclicPolyPatch::initUpdateTopology()
{
    polyPatch::initUpdateTopology();
}

void cyclicPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


const edgeList& cyclicPolyPatch::coupledPoints() const
{
    if (!coupledPointsPtr_)
    {
        // Look at cyclic patch as two halves, A and B.
        // Now all we know is that relative face index in halfA is same
        // as coupled face in halfB and also that the 0th vertex
        // corresponds.

        // From halfA point to halfB or -1.
        labelList coupledPoint(nPoints(), -1);

        for (label patchFaceA = 0; patchFaceA < size()/2; patchFaceA++)
        {
            const face& fA = localFaces()[patchFaceA];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (coupledPoint[patchPointA] == -1)
                {
                    const face& fB = localFaces()[patchFaceA + size()/2];

                    label indexB = (fB.size() - indexA) % fB.size();

                    coupledPoint[patchPointA] = fB[indexB];
                }
            }
        }

        coupledPointsPtr_ = new edgeList(nPoints());
        edgeList& connected = *coupledPointsPtr_;

        // Extract coupled points.
        label connectedI = 0;

        forAll(coupledPoint, i)
        {
            if (coupledPoint[i] != -1)
            {
                connected[connectedI++] = edge(i, coupledPoint[i]);
            }
        }

        connected.setSize(connectedI);

        if (debug)
        {
            Info<< "Writing file coupledPoints.obj with coordinates of "
                << "coupled points" << endl;

            OFstream str("coupledPoints.obj");
            label vertI = 0;

            forAll(connected, i)
            {
                const point& a = localPoints()[connected[i][0]];
                const point& b = localPoints()[connected[i][1]];

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledPointsPtr_;
}


const edgeList& cyclicPolyPatch::coupledEdges() const
{
    if (!coupledEdgesPtr_)
    {
        // Build map from points on halfA to points on halfB.
        const edgeList& pointCouples = coupledPoints();

        Map<label> aToB(2*pointCouples.size());

        forAll(pointCouples, i)
        {
            const edge& e = pointCouples[i];

            aToB.insert(e[0], e[1]);
        }

        // Map from edge on half A to points (in halfB indices)
        HashTable<label, edge, Hash<edge> > edgeMap(nEdges());

        for (label patchFaceA = 0; patchFaceA < size()/2; patchFaceA++)
        {
            const labelList& fEdges = faceEdges()[patchFaceA];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Convert edge end points to corresponding points on halfB
                // side.

                edge halfBEdge(aToB[e[0]], aToB[e[1]]);

                edgeMap.insert(halfBEdge, edgeI);
            }
        }

        coupledEdgesPtr_ = new edgeList(nEdges()/2);
        edgeList& coupledEdges = *coupledEdgesPtr_;
        label coupleI = 0;

        for (label patchFaceB = size()/2; patchFaceB < size(); patchFaceB++)
        {
            const labelList& fEdges = faceEdges()[patchFaceB];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Look up edge from HashTable.
                label halfAEdgeI = edgeMap[e];

                // Store correspondence
                coupledEdges[coupleI++] = edge(halfAEdgeI, edgeI);

                // Remove (since should be used only once)
                edgeMap.erase(e);
            }
        }


        // Some checks

        if (coupleI != coupledEdges.size())
        {
            FatalErrorIn("cyclicPolyPatch::coupledEdges() const")
                << "Problem : coupleI:" << coupleI
                << " number of edges per halfpatch:" << coupledEdges.size()
                << abort(FatalError);
        }

        forAll(coupledEdges, i)
        {
            const edge& e = coupledEdges[i];

            if (e[0] == e[1] || e[0] < 0 || e[1] < 0)
            {
                FatalErrorIn("cyclicPolyPatch::coupledEdges() const")
                    << "Problem : at position " << i
                    << " illegal couple:" << e
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            Info<< "Writing file coupledEdges.obj with centres of "
                << "coupled edges" << endl;

            OFstream str("coupledEdges.obj");
            label vertI = 0;

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                const point& a = edges()[e[0]].centre(localPoints());
                const point& b = edges()[e[1]].centre(localPoints());

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledEdgesPtr_;
}


//- Initialize ordering (on new mesh)
void cyclicPolyPatch::initOrder(const primitivePatch& pp) const
{}


//- Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool cyclicPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (pp.size() == 0)
    {
        // No faces, nothing to change.
        return false;
    }

    // Get geometric zones of patch by looking at normals

    boolList regionEdge(pp.nEdges(), false);

    // Added faces have no updated normal yet (points are vector::zero).
    // So use preMotionPoints
    // and hope for the best. Cyclics are limited to one transformation tensor
    // currently anyway (i.e. straight plane) so should not be too big a
    // problem.

    vectorField normals(pp.size());

    forAll(pp, faceI)
    {
        normals[faceI] = pp[faceI].normal(pp.points());
    }
    normals /= mag(normals) + VSMALL;

    const labelListList& eFaces = pp.edgeFaces();

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
    patchZones ppZones(pp, regionEdge);

    if (debug)
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

            const point cc(operator[](faceI).centre(pp.points()));

            writeOBJ(streams[zoneI], cc);

            nZoneFaces[zoneI]++;
        }

        Info<< "cyclicPolyPatch::order : Number of faces per zone:"
            << nZoneFaces << endl;
    }


    // Get the two halves and construct PrimitivePatches on them

    // From half0 to patch index
    labelList half0ToPatch(pp.size(), -1);

    // All faces of half0
    faceList half0Faces(pp.size());
    label n0Faces = 0;
    label face0 = -1;       // Arbitrary face on half0 (used for normal)

    labelList half1ToPatch(pp.size(), -1);
    faceList half1Faces(pp.size());
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


    primitiveFacePatch half0(half0Faces, pp.points());
    primitiveFacePatch half1(half1Faces, pp.points());


    // Get geometric data on both halves.
    const vector& n0 = normals[face0];
    const vector& n1 = normals[face1];

    pointField half0Ctrs(calcFaceCentres(half0, half0.points()));
    pointField anchors0(getAnchorPoints(half0, half0.points()));
    pointField half1Ctrs(calcFaceCentres(half1, half0.points()));

    if (mag(n0 & n1) < 1-SMALL)
    {
        if (debug)
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

        if (debug)
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

    if (debug)
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

        label newFaceI = half0FaceI + pp.size()/2;

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
