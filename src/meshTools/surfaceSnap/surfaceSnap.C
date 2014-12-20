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

#include "surfaceSnap.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "dictionary.H"
#include "features.H"
#include "physicalConstants.H"
#include "meshTools.H"
#include "octree.H"
#include "octreeDataEdges.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfaceSnap, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Calculate minimum length of edges on patch
Foam::scalarField Foam::surfaceSnap::minEdgeLength
(
    const pointField& points
) const
{
    const edgeList& edges = pp_.edges();
    const labelListList& pointEdges = pp_.pointEdges();

    scalarField minEdgeLen(points.size(), GREAT);

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        forAll(pEdges, pEdgeI)
        {
            const edge& e = edges[pEdges[pEdgeI]];

            scalar len = e.mag(points);

            minEdgeLen[pointI] = min(minEdgeLen[pointI], len);
        }
    }
    return minEdgeLen;
}


// Get average position of point in patch. Do not take points on features 
// into account so points move away from them.
Foam::point Foam::surfaceSnap::getAvgPos
(
    const Map<label>& patchToFeatureEdge,
    const Map<label>& patchToFeaturePoint,
    const bool allPoints,
    const label pointI
) const
{
    point avgPos(vector::zero);

    label nPos = 0;

    const labelList& pEdges = pp_.pointEdges()[pointI];

    forAll(pEdges, i)
    {
        const edge& e = pp_.edges()[pEdges[i]];

        label otherVertI = e.otherVertex(pointI);

        if
        (
            allPoints
         || (
                !patchToFeatureEdge.found(otherVertI)
             && !patchToFeaturePoint.found(otherVertI)
            )
        )
        {
            avgPos += points_[otherVertI];
            nPos++;
        }
    }

    if (nPos == 0)
    {
        return points_[pointI];
    }
    else
    {
        return avgPos/nPos;
    }
}


//// Given face with some vertices snapped find out if and which snapped vertex
//// should be unsnapped.
//void Foam::surfaceSnap::unsnapFace
//(
//    const label faceI,
//    Map<label>& patchToFeatureEdge
//) const
//{
//    const triSurface& surf = querySurf_.surface();
//
//    const pointField& surfPoints = surf.localPoints();
//
//
//    const face& f = pp_.localFaces()[faceI];
//
//    scalar oldArea = f.mag(points_);
//
//    // Get face with all points snapped
//    pointField facePts(f.size());
//
//    forAll(f, fp)
//    {
//        Map<label>::const_iterator iter = patchToFeatureEdge.find(f[fp]);
//
//        if (iter != patchToFeatureEdge.end())
//        {
//            label surfEdgeI = iter();
//
//            const edge& e = surf.edges()[surfEdgeI];
//
//            // Get nearest point on surfEdge
//            facePts[fp] =
//                line<point, const point&>
//                (
//                    surfPoints[e.start()],
//                    surfPoints[e.end()]
//                ).nearestDist(points_[f[fp]]).rawPoint();
//        }
//        else
//        {
//            facePts[fp] = points_[f[fp]];
//        }
//    }
//
//    // Create I-face to index facePts
//    face newFace(f.size());
//
//    forAll(f, fp)
//    {
//        newFace[fp] = fp;
//    }
//
//    scalar newArea = newFace.mag(facePts);
//
//    Info<< "Face:" << faceI << " oldArea:" << oldArea << " snappedArea:"
//        << newArea << endl;
//
//    if (newArea > SMALL && (newArea > 0.5 * oldArea))
//    {
//        // newArea isn't too bad.
//        return;
//    }
//
//    // Face with points snapped to feature would cause too small a face.
//    // See if it gets better if we leave out one of the features.
//
//    scalar maxArea = -GREAT;
//    label maxFp = -1;
//
//    forAll(f, fp)
//    {
//        if (patchToFeatureEdge.found(f[fp]))
//        {
//            label fp1 = (fp + 1) % f.size();
//            label fpMin1 = (fp == 0 ? f.size()-1 : fp-1);
//
//            vector a(facePts[fp1] - facePts[fp]);
//            vector b(facePts[fpMin1] - facePts[fp]);
//
//            if ((a & b) <= 0)
//            {
//                // Angle > 90 degrees.
//
//                // Temporarily change facePts[fp] to unsnapped value
//                point oldPt = facePts[fp];
//                facePts[fp] = points_[f[fp]];
//
//                scalar magArea = newFace.mag(facePts);
//
//                Info<< "faceI:" << faceI << " facePts:" << facePts << " area:"
//                    << magArea << endl;
//
//                if (magArea > maxArea)
//                {
//                    maxArea = magArea;
//                    maxFp = fp;
//                }
//
//                // Restore facePts to all snapped
//                facePts[fp] = oldPt;
//            }
//        }
//    }
//
//
//    if (maxFp != -1)
//    {
//        Info<< "Removing patch point " << f[maxFp] << " with area " << maxArea
//            << endl;
//
//        patchToFeatureEdge.erase(f[maxFp]);
//    }
//    else
//    {
//        WarningIn("Foam::surfaceSnap::unsnapFace")
//            << "Did not find angle on face > 90 degrees" << nl
//            << "Face:" << faceI << " unsnapped points:"
//            << IndirectList<point>(points_, f)
//            << " snapped points:" << facePts
//            << endl;
//    }
//}


Foam::point Foam::surfaceSnap::snapPosition
(
    const Map<label>& patchToFeatureEdge,
    const Map<label>& patchToFeaturePoint,
    const label patchPointI,
    const scalar projectDist
) const
{
    const triSurface& surf = querySurf_.surface();
    const pointField& surfPoints = surf.localPoints();

    if (patchToFeaturePoint.found(patchPointI))
    {
        return surfPoints[patchToFeaturePoint[patchPointI]];
    }
    else
    {
        point pos;

        Map<label>::const_iterator iter =
            patchToFeatureEdge.find(patchPointI);

        if (iter != patchToFeatureEdge.end())
        {
            // On feature edge or feature edge endpoint.

            if (relax_ < SMALL)
            {
                pos = points_[patchPointI];
            }
            else
            {
                // Get average position
                pos =
                    (1-relax_)
                  * points_[patchPointI]
                  + relax_
                  * getAvgPos
                    (
                        patchToFeatureEdge,
                        patchToFeaturePoint,
                        true,   // use all points for avg calculation
                        patchPointI
                    );
            }

            label surfEdgeI = iter();

            const edge& e = surf.edges()[surfEdgeI];

            // Get nearest point on surfEdge
            return
                line<point, const point&>
                (
                    surfPoints[e.start()],
                    surfPoints[e.end()]
                ).nearestDist(pos).rawPoint();
        }
        else
        {
            if (relax_ < SMALL)
            {
                pos = points_[patchPointI];
            }
            else
            {
                // Get average position
                pos =
                    (1-relax_)
                  * points_[patchPointI]
                  + relax_
                  * getAvgPos
                    (
                        patchToFeatureEdge,
                        patchToFeaturePoint,
                        false,              // Do not use points on
                                            // feature edges to get avg.
                        patchPointI
                    );
            }

            if (projectDist > SMALL)
            {
                // Snap pos to surface.
                tmp<pointField> tnearestSurfPts =
                    querySurf_.calcNearest
                    (
                        pointField(1, pos),
                        triSurfaceSearch::greatPoint
                    );
                const pointField& nearestSurfPts = tnearestSurfPts();

                vector displacement = nearestSurfPts[0] - pos;
                scalar dist = mag(displacement);

                if (dist <= projectDist)
                {
                    return nearestSurfPts[0];
                }
                else
                {
                    return pos + (projectDist/dist) * displacement;
                }
            }
            else
            {
                // No projection. Shortcut
                return pos;
            }
        }
    }
}


// Return true if points using patchPointI have valid faces when using
// points_.
bool Foam::surfaceSnap::validAreas(const label patchPointI) const
{
    const labelList& pFaces = pp_.pointFaces()[patchPointI];

    forAll(pFaces, i)
    {
        label faceI = pFaces[i];

        const face& f = pp_.localFaces()[faceI];

        scalar newArea = f.mag(points_);

        scalar oldArea = f.mag(pp_.localPoints());

        if (newArea < 0.2*oldArea)
        {
            return false;
        }
    }

    return true;
}


bool Foam::surfaceSnap::snapPoint
(
    const label patchPointI,
    const scalar projectDist,
    const point& wantedPos
)
{
    // Save old position and set point to new position
    point oldPos = points_[patchPointI];
    points_[patchPointI] = wantedPos;

    // Check if areas are ok with this new position (usually the case)
    if (validAreas(patchPointI))
    {
        return true;
    }

    // Try with limited movement
    vector disp = 0.5*(wantedPos - oldPos);

    for (label iter = 0; iter < 4; iter++)
    {
        points_[patchPointI] = oldPos + disp;

        if (debug&2)
        {
            Info<< "Failed to move pointI:" << patchPointI << " on to "
                << wantedPos << " from old position " << oldPos << nl
                << "Retrying with " << points_[patchPointI] << endl;
        }

        if (validAreas(patchPointI))
        {
            return true;
        }

        // Retry with smaller displacement
        disp *= 0.5;
    }


    if (debug&2)
    {
        Info<< "Fully failed to move pointI:" << patchPointI << " on to "
            << wantedPos << " from old position " << oldPos << endl;
    }

    // Restore old point position
    points_[patchPointI] = oldPos;

    return false;
}


// Find unsnapped patchPoint that is closest to a feature edge.
Foam::label Foam::surfaceSnap::closestToFeatEdge
(
    const labelList& surfEdges,
    const pointField& surfEdgePoints,
    const boolList& snappedPoint    // whether point has been snapped already
) const
{
    // Find (unsnapped) point nearest to feature edge.
    scalar minDist = GREAT;
    label minPointI = -1;

    forAll(points_, patchPointI)
    {
        if (surfEdges[patchPointI] != -1 && !snappedPoint[patchPointI])
        {
            scalar dist =
                mag
                (
                    points_[patchPointI]
                  - surfEdgePoints[patchPointI]
                );

            if (dist < minDist)
            {
                minDist = dist;
                minPointI = patchPointI;
            }
        }
    }
    return minPointI;
}


// Coming from edgeI, cross fromPointI to next point such that edge max. aligns
// with surface edgeI. Return -1 or edge label.
Foam::label Foam::surfaceSnap::maxAlignedEdge
(
    const labelList& surfEdges,
    const boolList& snappedPoint,   // whether point has been snapped already
    const label fromEdgeI,
    const label fromPointI
) const
{
    const triSurface& surf = querySurf_.surface();

    // Determine vector on surface
    label surfEdgeI = surfEdges[fromPointI];
    vector surfVec = surf.edges()[surfEdgeI].vec(surf.localPoints());
    surfVec /= mag(surfVec);

    // Determine vector of fromEdge
    vector fromVec;
    if (fromEdgeI != -1)
    {
        label otherVertI = pp_.edges()[fromEdgeI].otherVertex(fromPointI);

        fromVec = points_[fromPointI] - points_[otherVertI];
        fromVec /= mag(fromVec);
    }

    //
    // Find connected point that
    // - makes >90 degr angle with fromEdgeI
    // - makes max angle with surfVec.

    const labelList& pEdges = pp_.pointEdges()[fromPointI];

    scalar maxCos = -GREAT;
    label maxEdgeI = -1;

    forAll(pEdges, i)
    {
        label edgeI = pEdges[i];

        if (edgeI != fromEdgeI)
        {
            const edge& e = pp_.edges()[edgeI];

            label otherVertI = e.otherVertex(fromPointI);

            if (surfEdges[otherVertI] != -1 && !snappedPoint[otherVertI])
            {
                // otherVertI close to feature. Check angles.
                vector eVec = points_[otherVertI] - points_[fromPointI];
                eVec /= mag(eVec);

                if (fromEdgeI == -1 || ((eVec & fromVec) > 0.7))
                {
                    // Angle with previous edge ok.
                    // Check angle with surface edge
                    scalar cosAngle = mag(eVec & surfVec);

                    if (cosAngle > maxCos)
                    {
                        maxCos = cosAngle;
                        maxEdgeI = edgeI;
                    }
                }
            }
        }
    }

    return maxEdgeI;
}


// Snap all unsnapped points that are close to feature edges onto
// these feature edges. It tries to stitch whole lines of connected
// mesh edges. Updates snappedPoint for all points snapped.
void Foam::surfaceSnap::stitchToFeatureEdges
(
    const scalarField& snapDist,
    boolList& snappedPoint
)
{
    // Get nearest feature edge
    labelList surfEdges;
    labelList surfEndPoints;
    pointField surfEdgePoints;

    surfFeatures_.nearestSurfFeatures
    (
        points_,
        snapDist,
        surfEdges,
        surfEndPoints,
        surfEdgePoints
    );


    OFstream snapStream("snap.obj");
    label vertI = 0;

    do
    {
        // Find (unsnapped) point nearest to feature edge.
        label pointI =
            closestToFeatEdge
            (
                surfEdges,
                surfEdgePoints,
                snappedPoint
            );

        if (pointI == -1)
        {
            if (debug)
            {
                Info<< "Cannot find unsnapped patchPoint to snap to feature"
                    << " edge. Finishing ..." << nl << endl;
            }
            break;
        }

        // Start off with illegal edge (handled in maxAlignedEdge)
        label edgeI = -1;

        do
        {
            // Stitch point
            if (debug)
            {
                Info<< "Snapping " << pointI << " from pos:"
                    << points_[pointI] << " to " << surfEdgePoints[pointI]
                    << endl;

                meshTools::writeOBJ(snapStream, points_[pointI]);
                vertI++;
                meshTools::writeOBJ(snapStream, surfEdgePoints[pointI]);
                vertI++;
                snapStream << "l " << vertI-1 << ' ' << vertI << nl;
            }

            snapPoint(pointI, snapDist[pointI], surfEdgePoints[pointI]);

            // Mark point off as visited (even when snapping actually has
            // failed)
            snappedPoint[pointI] = true;


            // Find connected point which is in same direction as surface edge 
            // minPointI is snapped to. Note:
            // - minPointI might be snapped to edge end point. In this case get
            //   max correspondence between connected edges and vector between
            //   snapped points.
            // - if not on edge endpoint two solutions:
            //   - get current surface edge vector to compare next point to
            //   - get surface edge belonging to next point to compare vector
            //     to.
            //   (note: is question of where to sample surface - at start or
            //    at end or in middle)


            // Step to next edge.
            edgeI =
                maxAlignedEdge
                (
                    surfEdges,
                    snappedPoint,
                    edgeI,
                    pointI
                );


            if (edgeI == -1)
            {
                Info<< "Did not find point connected to " << pointI
                    << " aligned with surface edge vector "
                    << endl;
                break;
            }

            // Step to next point.
            label nextPointI = pp_.edges()[edgeI].otherVertex(pointI);

            if (debug)
            {
                Info<< "From point:" << pointI << " pos:" << points_[pointI]
                    << " to point:" << nextPointI
                    << " pos:" << points_[nextPointI] << endl;
            }

            pointI = nextPointI;
        }
        while (!snappedPoint[pointI]);

    }
    while (true);
}


void Foam::surfaceSnap::snapToFeatureEdges
(
    const scalarField& snapDist,
    const Map<label>& patchToFeatureEdge,
    boolList& snappedPoint
)
{
    labelList unsnapped(patchToFeatureEdge.toc());
    pointField unsnappedPoints(unsnapped.size());

    forAll(unsnapped, i)
    {
        unsnappedPoints[i] = points_[unsnapped[i]];
    }

    // Get nearest feature edge (so discard patchToFeatureEdge surface edge
    // information since this is the reverse: nearest patchPoint to a given
    // edge; we need to find the edge where the distance is nearest)
    labelList featEdges;
    labelList featEndPoints;
    pointField featEdgePoints;

    surfFeatures_.nearestSurfFeatures
    (
        unsnappedPoints,
        snapDist,
        featEdges,
        featEndPoints,
        featEdgePoints
    );

    forAll(unsnapped, i)
    {
        label patchPointI = unsnapped[i];

        if (featEdges[i] != -1)
        {
            if
            (
                snapPoint
                (
                    patchPointI,
                    snapDist[patchPointI],
                    featEdgePoints[i]
                )
            )
            {
                snappedPoint[patchPointI] = true;
            }
        }
    }
}



Foam::label Foam::surfaceSnap::snapFeaturesPost
(
    const scalarField& snapDist,
    const Map<label>& patchToFeatureEdge,
    boolList& snappedPoint
)
{
    // Build tree out of all pp_ edges.
    labelList edgeLabels(pp_.nEdges());
    forAll(edgeLabels, i)
    {
        edgeLabels[i] = i;
    }
    octreeDataEdges shapes(pp_.edges(), points_, edgeLabels);

    treeBoundBox overallBb(points_);

    octree<octreeDataEdges> ppTree
    (
        overallBb,  // overall search domain
        shapes,     // all information needed to do checks on cells
        1,          // min levels
        20.0,       // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity; n/a since no bounding boxes.
    );

    //
    // Now find nearest pp edge for every sampled feature edge.
    //

    const triSurface& surf = querySurf_.surface();
    const labelHashSet& surfEdgeLabels = surfFeatures_.edgeLabels();
    const edgeList& surfEdges = surf.edges();
    const pointField& surfPoints = surf.localPoints();

    // Local search area
    const scalar maxSearch = max(snapDist);
    const vector span(maxSearch, maxSearch, maxSearch);

    // All unsnapped points.
    DynamicList<label> unsnapped;
    DynamicList<point> unsnappedPoints;

    OFstream ppStream("unsnapped.obj");
    label ppVertI = 0;

    OFstream surfStream("unmapped.obj");
    label surfVertI = 0;

    for
    (
        labelHashSet::const_iterator iter = surfEdgeLabels.begin();
        iter != surfEdgeLabels.end();
        ++iter
    )
    {
        label surfEdgeI = iter.key();

        const edge& e = surfEdges[surfEdgeI];

        // Normalized edge vector
        vector eVec = e.vec(surfPoints);
        scalar eMag = mag(eVec);
        eVec /= eMag;


        //
        // Sample along edge
        //

        // Whether to exit sampling current edge
        bool exit = false;

        // Coordinate along edge (0 .. eMag)
        scalar s = 0.0;
        // Sampling distance along edge
        scalar sOffset = 0.1*eMag;

        while (true)
        {
            point edgePoint(surfPoints[e.start()] + s*eVec);

            treeBoundBox tightest(edgePoint - span, edgePoint + span);
            scalar tightestDist = GREAT;

            label patchEdgeI = ppTree.findNearest
            (
                edgePoint,
                tightest,
                tightestDist
            );

            if (patchEdgeI == -1)
            {
                WarningIn("Foam::surfaceSnap::snapFeaturesPost")
                    << "Did not find patch edge close to point "
                    << edgePoint << " on surface feature edge "
                    << surfEdgeI << endl;
            }
            else
            {
                const edge& patchE = pp_.edges()[patchEdgeI];

                bool v0Snapped = patchToFeatureEdge.found(patchE[0]);
                bool v1Snapped = patchToFeatureEdge.found(patchE[1]);


                if (v0Snapped && v1Snapped)
                {
                    //scalar edgeSampleDist =
                    //    0.5
                    //  * (
                    //        snapDist[patchE[0]]
                    //      + snapDist[patchE[1]]
                    //    );
                    //
                    //sOffset = max(0.1*eMag, edgeSampleDist);
                }
                else if (v0Snapped)
                {
                    unsnapped.append(patchE[1]);
                    unsnappedPoints.append(points_[patchE[1]]);
                }
                else if (v1Snapped)
                {
                    unsnapped.append(patchE[0]);
                    unsnappedPoints.append(points_[patchE[0]]);
                }
                else
                {
                    //Info<< "UnsnappedSurfedge:" << surfEdgeI << " s:" << s
                    //    << " nearest ppEdge:"
                    //    << patchE << " points " << points_[patchE[0]]
                    //    << points_[patchE[1]]
                    //    << " dist:" << tightestDist << endl;

                    // Found edge close to surface feature but one or both
                    // of the edge end points is not snapped.
                    meshTools::writeOBJ(ppStream, points_[patchE[0]]);
                    ppVertI++;
                    meshTools::writeOBJ(ppStream, points_[patchE[1]]);
                    ppVertI++;
                    ppStream << "l " << ppVertI-1 << ' ' << ppVertI << endl;

                    meshTools::writeOBJ(surfStream, surfPoints[e[0]]);
                    surfVertI++;
                    meshTools::writeOBJ(surfStream, surfPoints[e[1]]);
                    surfVertI++;
                    surfStream << "l " << surfVertI-1 << ' ' << surfVertI
                        << endl;
                }
            }

            if (exit)
            {
                break;
            }

            // Step to next sample point
            s += sOffset;

            if (s >= 0.95*eMag)
            {
                // Do one more sample, at edge endpoint
                s = eMag;
                exit = true;
            }
        }
    }

    unsnapped.shrink();
    unsnappedPoints.shrink();

    // Transfer unsnappedPoints to pointField
    pointField allUnsnappedPoints;
    allUnsnappedPoints.transfer(unsnappedPoints);
    unsnappedPoints.clear();

    // Get nearest feature edge
    labelList featEdges;
    labelList featEndPoints;
    pointField featEdgePoints;

    surfFeatures_.nearestSurfFeatures
    (
        allUnsnappedPoints,
        snapDist,
        featEdges,
        featEndPoints,
        featEdgePoints
    );

    label nSnapped = 0;

    forAll(unsnapped, i)
    {
        label patchPointI = unsnapped[i];

        if (featEdges[i] != -1)
        {
            if
            (
                snapPoint
                (
                    patchPointI,
                    snapDist[patchPointI],
                    featEdgePoints[i]
                )
            )
            {
                if (debug)
                {
                    Info<< "Unsnapped PatchPoint " << patchPointI << " pos:"
                        << points_[patchPointI]
                        << " to position " << points_[patchPointI] << endl;
                }
                snappedPoint[patchPointI] = true;
                nSnapped++;
            }
            else
            {
                // Snapping would create illegal face. Don't do anything
                // unless absolutely needed.

                if (forceSnap_)
                {
                    points_[patchPointI] = featEdgePoints[i];
                    snappedPoint[patchPointI] = true;
                    nSnapped++;

                    if (debug)
                    {
                        Info<< "Unsnapped PatchPoint " << patchPointI << " pos:"
                            << points_[patchPointI]
                            << " forced to position " << points_[patchPointI]
                            << endl;
                    }
                }
            }
        }
    }

    if (debug)
    {
        Info<< "*** nSnapped:" << nSnapped << endl;
    }

    return nSnapped;
}
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::surfaceSnap::surfaceSnap
(
    const primitivePatch& pp,
    const triSurfaceSearch& querySurf,
    const features& surfFeatures,
    const dictionary& dict
)
:
    pp_(pp),
    querySurf_(querySurf),
    surfFeatures_(surfFeatures),
    relax_(readScalar(dict.lookup("relax"))),
    snapTol_(readScalar(dict.lookup("snapTol"))),
    projectTol_(readScalar(dict.lookup("projectTol"))),
    forceSnap_(readBool(dict.lookup("forceSnap"))),
    points_(pp.localPoints())
{
    if (relax_ < 0 || relax_ > 1)
    {
        FatalErrorIn
        (
            "surfaceSnap::surfaceSnap(const primitivePatch&"
            "const triSurfaceSearch&, const features&, const dictionary&"
        )   << "Relaxation parameter out of range [0..1] : " << relax_
            << exit(FatalError);
    }

    snap();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfaceSnap::snap()
{
    // Local length scale/max distance to move (per patch point).
    //scalarField localDist = minEdgeLength(pp_.localPoints());
    scalarField localDist = minEdgeLength(points_);

    scalarField snapDist = snapTol_ * localDist;
    scalarField projectDist = projectTol_ * localDist;

    if (debug)
    {
        Info<< "Local length scale:" << endl
            << "    min:" << min(localDist) << endl
            << "    max:" << max(localDist) << endl
            << "    avg:" << sum(localDist)/localDist.size() << endl
            << endl;
    }

    // Find labels of patch point for each surfFeatPoint.
    // Note: only want one patchPoint for each surfFeatPoint since
    // don't want to collapse hexes.
    Map<label> patchToFeaturePoint
    (
        surfFeatures_.nearestSamples
        (
            surfFeatures_.pointLabels(),
            points_,
            snapDist
        )
    );

    if (debug)
    {
        Info<< "Found:\n"
            << "    patch points:" << patchToFeaturePoint.size() << nl
            << "    close to pointFeatures:"
            << surfFeatures_.pointLabels().size() << nl << endl;
    }

    // Find labels of patch points that will be snapped to surface feature edge.
    Map<label> patchToFeatureEdge
    (
        surfFeatures_.nearestSamplesToFeatEdges
        (
            points_,
            0.5*localDist,  // Sampling distance based on OLD position
                            // (is to prevent sampling dist becoming 0)
            snapDist
        )
    );

    if (debug)
    {
        Info<< "Found:\n"
            << "    patch points:" << patchToFeatureEdge.size() << nl
            << "    close to edgeFeatures:"
            << surfFeatures_.edgeLabels().size() << nl << endl;
    }

    boolList snappedPoint(pp_.nPoints(), false);

    //patchToFeatureEdge.clear();
    //stitchToFeatureEdges(snapDist, snappedPoint);

    snapToFeatureEdges(snapDist, patchToFeatureEdge, snappedPoint);


    //// Check faces where all verts are snapped to features since they quite
    //// possibly might all be snapped to the same feature and thus create
    //// a sliver face.
    //labelList nFeatsPerFace(pp_.size(), 0);
    //
    //for
    //(
    //    Map<label>::const_iterator iter = patchToFeatureEdge.begin();
    //    iter != patchToFeatureEdge.end();
    //    ++iter
    //)
    //{
    //    const labelList& pFaces = pp_.pointFaces()[iter.key()];
    //
    //    forAll(pFaces, i)
    //    {
    //        nFeatsPerFace[pFaces[i]]++;
    //    }
    //}
    //
    //forAll(nFeatsPerFace, faceI)
    //{
    //    if (nFeatsPerFace[faceI] == pp_[faceI].size())
    //    {
    //        // All points of face would get snapped.
    //        Info<< "Face:" << faceI << " verts:" << pp_[faceI]
    //            << " nFeats:" << nFeatsPerFace[faceI] << endl;
    //
    //        unsnapFace(faceI, patchToFeatureEdge);
    //    }
    //}


    forAll(points_, patchPointI)
    {
        if (debug && (patchPointI % 1000) == 0)
        {
            Info<< "Snapping point " << patchPointI << endl;
        }

        if (!snappedPoint[patchPointI])
        {
            point wantedPos
            (
                snapPosition
                (
                    patchToFeatureEdge,
                    patchToFeaturePoint,
                    patchPointI,
                    projectDist[patchPointI]
                )
            );

            bool okFaces = snapPoint
            (
                patchPointI,
                projectDist[patchPointI],
                wantedPos
            );

            snappedPoint[patchPointI] = true;

            if (!okFaces)
            {
                Info<< "Did not project patchPointI:" << patchPointI
                    << " pos:" << points_[patchPointI]
                    << " to position " << wantedPos
                    << " since would create illegal face" << endl;
            }
        }
    }


    // Check whether all of the surface feature edges is 'covered' and
    // correct points if nessecary
    label oldSnapped = -1;

    while (true)
    {    
        label n = snapFeaturesPost(snapDist, patchToFeatureEdge, snappedPoint);

        if (n == oldSnapped)
        {
            break;
        }
        oldSnapped = n;
    }

    return true;
}


// ************************************************************************* //
