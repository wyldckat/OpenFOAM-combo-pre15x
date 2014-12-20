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

#include "features.H"
#include "triSurface.H"
#include "octree.H"
#include "octreeDataEdges.H"
#include "octreeDataPoint.H"
#include "meshTools.H"
#include "linePointRef.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::features, 0);


//- Get nearest point on edge and classify position on edge.
Foam::pointIndexHit Foam::features::edgeNearest
(
    const point& start,
    const point& end,
    const point& sample
)
{
    pointHit eHit = linePointRef(start, end).nearestDist(sample);

    // Classification of position on edge.
    label endPoint;

    if (eHit.hit())
    {
        // Nearest point is on edge itself.
        // Note: might be at or very close to endpoint. We should use tolerance
        // here.
        endPoint = -1;
    }
    else
    {
        // Nearest point has to be one of the end points. Find out
        // which one.
        if
        (
            mag(eHit.rawPoint() - start)
          < mag(eHit.rawPoint() - end)
        )
        {
            endPoint = 0;
        }
        else
        {
            endPoint = 1;
        }
    }

    return pointIndexHit(eHit.hit(), eHit.rawPoint(), endPoint);
}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Returns next feature edge connected to pointI with correct value.
Foam::label Foam::features::nextFeatEdge
(
    const boolList& isFeatEdge,
    const labelList& featVisited,
    const label unsetVal,
    const label prevEdgeI,
    const label vertI
) const
{
    const labelList& pEdges = surf_.pointEdges()[vertI];

    label nextEdgeI = -1;

    forAll(pEdges, i)
    {
        label edgeI = pEdges[i];

        if
        (
            edgeI != prevEdgeI
         && isFeatEdge[edgeI]
         && featVisited[edgeI] == unsetVal
        )
        {
            if (nextEdgeI == -1)
            {
                nextEdgeI = edgeI;
            }
            else
            {
                // More than one feature edge to choose from. End of segment.
                return -1;
            }
        }
    }

    return nextEdgeI;
}


// Finds connected feature edges by walking from prevEdgeI in direction of
// prevPointI. Marks feature edges visited in featVisited by assigning them
// the current feature line number. Returns cumulative length of edges walked.
// Works in one of two modes:
// - mark : step to edges with featVisited = -1. 
//          Mark edges visited with currentFeatI.
// - clear : step to edges with featVisited = currentFeatI
//           Mark edges visited with -2 and erase from feature edges.
Foam::scalar Foam::features::walkSegment
(
    const bool mark,
    const boolList& isFeatEdge,
    const label startEdgeI,
    const label startPointI,
    const label currentFeatI,     
    labelList& featVisited
)
{
    label edgeI = startEdgeI;

    label vertI = startPointI;


    //
    // Now we have:
    //    edgeI : first edge on this segment
    //    vertI : one of the endpoints of this segment
    //
    // Start walking, marking off edges (in featVisited)
    // as we go along.
    //

    label unsetVal;

    if (mark)
    {
        unsetVal = -1;
    }
    else
    {
        unsetVal = currentFeatI;
    }


    scalar visitedLength = 0.0;

    label nVisited = 0;

    do
    {
        // Step to next feature edge with value unsetVal
        edgeI = nextFeatEdge(isFeatEdge, featVisited, unsetVal, edgeI, vertI);

        if (edgeI == -1 || edgeI == startEdgeI)
        {
            break;
        }

        // Mark with current value. If in clean mode also remove feature edge.

        if (mark)
        {
            featVisited[edgeI] = currentFeatI;
        }
        else
        {
            featVisited[edgeI] = -2;
            edgeLabels_.erase(edgeI);
        }


        // Step to next vertex on edge
        const edge& e = surf_.edges()[edgeI];

        vertI = e.otherVertex(vertI);


        // Update cumulative length
        visitedLength += e.mag(surf_.localPoints());

        nVisited++;

        if (nVisited > surf_.nEdges())
        {
            Warning<< "walkSegment : reached iteration limit in walking "
                << "feature edges on surface from edge:" << startEdgeI
                << " vertex:" << startPointI << nl
                << "Returning with large length" << endl;

            return GREAT;
        }
    }
    while (true);

    return visitedLength;
}



void Foam::features::filterFeatureEdges()
{
    // Convert feature edges to boolList for speed.
    boolList isFeatureEdge(surf_.nEdges(), false);
    for
    (
        labelHashSet::const_iterator iter = edgeLabels_.begin();
        iter != edgeLabels_.end();
        ++iter
    )
    {
        isFeatureEdge[iter.key()] = true;
    }
    

    // Mark feature edges according to connected lines.
    // -1: unassigned
    // -2: part of too small a feature line
    // >0: feature line number
    labelList featLines(surf_.nEdges(), -1);

    // Current featureline number.
    label featI = 0;

    do
    {
        // Find unset featureline
        label startEdgeI = -1;

        forAll(isFeatureEdge, edgeI)
        {
            if (isFeatureEdge[edgeI] && featLines[edgeI] == -1)
            {
                startEdgeI = edgeI;

                break;
            }
        }


        if (startEdgeI == -1)
        {
            // No unset feature edge found. All feature edges have line number
            // assigned.
            break;
        }

        featLines[startEdgeI] = featI;

        const edge& startEdge = surf_.edges()[startEdgeI];

        // Walk 'left' and 'right' from startEdge.
        scalar visitedLength =
            walkSegment
            (
                true,           // 'mark' mode
                isFeatureEdge,
                startEdgeI,     // edge, not yet assigned to a featureLine
                startEdge[0],   // start of edge
                featI,          // Mark value
                featLines       // Mark for all edges
            );

        visitedLength +=
            walkSegment
            (
                true,
                isFeatureEdge,
                startEdgeI,
                startEdge[1],   // end of edge
                featI,
                featLines
            );

        if (visitedLength < minDim_)
        {
            // Rewalk same route (recognizable by featLines == featI)
            // to reset featLines.

            featLines[startEdgeI] = -2;
            edgeLabels_.erase(startEdgeI);

            walkSegment
            (
                false,          // 'clean' mode
                isFeatureEdge,
                startEdgeI,     // edge, not yet assigned to a featureLine
                startEdge[0],   // start of edge
                featI,          // Unset value
                featLines       // Mark for all edges
            );

            walkSegment
            (
                false,
                isFeatureEdge,
                startEdgeI,
                startEdge[1],   // end of edge
                featI,
                featLines
            );
        }
        else
        {
            featI++;
        }
    }
    while (true);
}


void Foam::features::calculateFeatures()
{
    const labelListList& edgeFaces = surf_.edgeFaces();
    const labelListList& faceEdges = surf_.faceEdges();
    const vectorField& faceAreas = surf_.faceNormals();
    const labelListList& pointEdges = surf_.pointEdges();
    const List<labelledTri>& localFaces = surf_.localFaces();
    const edgeList& edges = surf_.edges();
    const vectorField& pointNormals = surf_.pointNormals();


    //
    // Feature edges.
    //

    if (edgeMinCos_ > (-1 + SMALL))
    {
        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = edgeFaces[edgeI];

            if (eFaces.size() == 2)
            {
                label face0I = eFaces[0];
                label face1I = eFaces[1];

                const labelledTri& f0 = localFaces[face0I];

                const labelledTri& f1 = localFaces[face1I];

                if (f0.region() != f1.region())
                {
                    edgeLabels_.insert(edgeI);
                }
                else
                {
                    const edge& e = edges[edgeI];

                    bool inOrder0 =
                        (f0[0] == e[0] && f0[1] == e[1])
                     || (f0[1] == e[0] && f0[2] == e[1])
                     || (f0[2] == e[0] && f0[0] == e[1]);

                    bool inOrder1 =
                        (f1[0] == e[0] && f1[1] == e[1])
                     || (f1[1] == e[0] && f1[2] == e[1])
                     || (f1[2] == e[0] && f1[0] == e[1]);

                    const vector& n0 = faceAreas[face0I];

                    const vector& n1 = faceAreas[face1I];

                    scalar faceAngle = n0 & n1;

                    if (inOrder0 ^ inOrder1)
                    {
                        // Differing ordering of edge points. Ok.
                    }
                    else
                    {
                        // Same ordering of edge points
                        faceAngle = -faceAngle;
                    }

                    if (faceAngle < edgeMinCos_)
                    {
                        edgeLabels_.insert(edgeI);
                    }
                }
            }
            else
            {
                //What to do here? Edge with more than two faces connected to
                // it.
                //edgeLabels_.insert(edgeI);
            }
        }

        if (minDim_ > SMALL)
        {
            filterFeatureEdges();
        }
    }


    //
    // Feature points.
    //

    if (pointMinCos_ > (-1 + SMALL))
    {
        // Calculate normals on edges as avg of normals on faces.
        vectorField edgeNormals(edges.size(), vector::zero);

        forAll(faceAreas, faceI)
        {
            vector n = faceAreas[faceI];
            n /= mag(n);

            const labelList& fEdges = faceEdges[faceI];

            forAll(fEdges, i)
            {
                edgeNormals[fEdges[i]] += n;
            }
        }
        edgeNormals /= mag(edgeNormals);

        forAll(pointEdges, pointI)
        {
            const labelList& myEdges = pointEdges[pointI];

            //label nFeatEdges = 0;
            //
            //forAll(myEdges, i)
            //{
            //    if (edgeLabels_.found(myEdges[i]))
            //    {
            //        nFeatEdges++;
            //
            //        if (nFeatEdges > 1)
            //        {
            //            break;
            //        }
            //    }
            //}

            //if (nFeatEdges == 1)
            //{
            //    // Only one connected featureEdge -> endpoint of feature line.
            //    pointLabels_.insert(pointI);
            //}
            //else
            //{
                // Check if pointNormal varies too much.
                const vector& n = pointNormals[pointI];

                forAll(myEdges, i)
                {
                    if (mag(n & edgeNormals[myEdges[i]]) < pointMinCos_)
                    {
                        pointLabels_.insert(pointI);

                        break;
                    }
                }
            //}
        }

        if (debug)
        {
            // Dump to obj file
            writeOBJ();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from single angle only.
Foam::features::features(const triSurface& surf, const scalar minCos)
:
    surf_(surf),
    pointMinCos_(minCos),
    edgeMinCos_(minCos),
    minDim_(0),
    pointLabels_(),
    edgeLabels_()
{
    calculateFeatures();
}


// Construct from separate angles for points and edges.
Foam::features::features
(
    const triSurface& surf,
    const scalar pointMinCos,
    const scalar edgeMinCos,
    const scalar minDim
)
:
    surf_(surf),
    pointMinCos_(pointMinCos),
    edgeMinCos_(edgeMinCos),
    minDim_(minDim),
    pointLabels_(surf.nPoints()/100 + 100),
    edgeLabels_(surf.nEdges()/100 + 100)
{
    calculateFeatures();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::features::writeOBJ() const
{
    const pointField& localPoints = surf_.localPoints();
    const edgeList& edges = surf_.edges();

    //
    // Featurepoints : dump to obj file
    //

    Info<< endl << "Dumping surface feature points (angle "
        << 180.0 * Foam::acos(pointMinCos_)/physicalConstant::pi
        << ") to surfFeatPoints.obj\n"
        << "View this Lightwave-OBJ file with e.g. javaview\n"
        << endl;

    OFstream pointStream("surfFeatPoints.obj");

    for
    (
        labelHashSet::const_iterator iter = pointLabels_.begin();
        iter != pointLabels_.end();
        ++iter
    )
    {
        meshTools::writeOBJ(pointStream, localPoints[iter.key()]);
    }


    //
    // Featureedges : dump to obj file
    //

    Info<< endl << "Dumping surface feature edges (angle "
        << 180.0 * Foam::acos(edgeMinCos_)/physicalConstant::pi
        << ") to surfFeatEdges.obj\n"
        << "View this Lightwave-OBJ file with e.g. javaview\n"
        << endl;

    // Write compact.
    label newPointI = 0;
    labelList pointMap(localPoints.size(), -1);

    OFstream edgeStream("surfFeatEdges.obj");

    // Write used points and construct map.
    for
    (
        labelHashSet::const_iterator iter = edgeLabels_.begin();
        iter != edgeLabels_.end();
        ++iter
    )
    {
        const edge& e = edges[iter.key()];

        forAll(e, i)
        {
            label pointI = e[i];

            if (pointMap[pointI] == -1)
            {
                meshTools::writeOBJ(edgeStream, localPoints[pointI]);

                pointMap[pointI] = newPointI++;
            }
        }
    }

    // Write edges using new numbering.
    for
    (
        labelHashSet::const_iterator iter = edgeLabels_.begin();
        iter != edgeLabels_.end();
        ++iter
    )
    {
        const edge& e = edges[iter.key()];

        edgeStream
            << "l " << pointMap[e[0]] + 1
            << ' ' << pointMap[e[1]] + 1
            << endl;
    }
}


// Get nearest vertex on patch for every point of surf in pointSet.
Foam::Map<Foam::label> Foam::features::nearestSamples
(
    const labelHashSet& points,
    const pointField& samples,
    const scalarField& maxDist
) const
{
    // Build tree out of all samples.
    // Note: shapes holds reference!
    octreeDataPoint shapes(samples);

    treeBoundBox overallBb(samples);

    octree<octreeDataPoint> ppTree
    (
        overallBb,  // overall search domain
        shapes,     // all information needed to do checks on cells
        1,          // min levels
        20.0,       // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity; n/a since no bounding boxes.
    );

    // From patch point to surface point
    Map<label> nearest(2*points.size());

    const pointField& surfPoints = surf_.localPoints();

    for
    (
        labelHashSet::const_iterator iter = points.begin();
        iter != points.end();
        ++iter
    )
    {
        label surfPointI = iter.key();

        treeBoundBox tightest = overallBb;
        scalar tightestDist = Foam::GREAT;

        label sampleI = ppTree.findNearest
        (
            surfPoints[surfPointI],
            tightest,
            tightestDist
        );

        if (sampleI == -1)
        {
            FatalErrorIn("features::nearestSamples") << "Problem for point "
                << surfPointI << " in tree " << ppTree.octreeBb()
                << abort(FatalError);
        }

        if (mag(samples[sampleI] - surfPoints[surfPointI]) < maxDist[sampleI])
        {
            nearest.insert(sampleI, surfPointI);
        }
    }


    if (debug)
    {
        //
        // Dump to obj file
        //

        Info<< endl
            << "Dumping nearest surface feature points to nearestSamples.obj"
            << endl
            << "View this Lightwave-OBJ file with e.g. javaview" << endl
            << endl;

        OFstream objStream("nearestSamples.obj");

        label vertI = 0;
        for
        (
            Map<label>::const_iterator iter = nearest.begin();
            iter != nearest.end();
            ++iter
        )
        {
            meshTools::writeOBJ(objStream, samples[iter.key()]); vertI++;
            meshTools::writeOBJ(objStream, surfPoints[iter()]); vertI++;
            objStream<< "l " << vertI-1 << ' ' << vertI << endl;
        }
    }

    return nearest;
}


// Get nearest vertex on patch for regularly sampled points along the feature
// edges. Return map from sample to feature edge.
Foam::Map<Foam::label> Foam::features::nearestSamplesToFeatEdges
(
    const pointField& samples,
    const scalarField& sampleDist,
    const scalarField& maxDist,
    const scalar minSampleDist
) const
{
    const pointField& surfPoints = surf_.localPoints();
    const edgeList& surfEdges = surf_.edges();

    scalar maxSearch = max(maxDist);
    vector span(maxSearch, maxSearch, maxSearch);

    // shapes holds reference!
    octreeDataPoint shapes(samples);

    treeBoundBox overallBb(samples);

    octree<octreeDataPoint> ppTree
    (
        overallBb,  // overall search domain
        shapes,     // all information needed to do checks on cells
        1,          // min levels
        20.0,       // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity; n/a since no bounding boxes.
    );

    // From patch point to surface edge.
    Map<label> nearest(2*edgeLabels_.size());

    label i = 0;

    for
    (
        labelHashSet::const_iterator iter = edgeLabels_.begin();
        iter != edgeLabels_.end();
        ++iter
    )
    {
        label surfEdgeI = iter.key();

        const edge& e = surfEdges[surfEdgeI];

        if (debug && (i % 1000) == 0)
        {
            Info<< "looking at surface feature edge " << surfEdgeI
                << " verts:" << e << " points:" << surfPoints[e[0]]
                << ' ' << surfPoints[e[1]] << endl;
        }

        // Normalized edge vector
        vector eVec = e.vec(surfPoints);
        scalar eMag = mag(eVec);
        eVec /= eMag;


        //
        // Sample along edge
        //

        bool exit = false;

        // Coordinate along edge (0 .. eMag)
        scalar s = 0.0;

        while (true)
        {
            point edgePoint(surfPoints[e.start()] + s*eVec);

            treeBoundBox tightest(edgePoint - span, edgePoint + span);
            scalar tightestDist = Foam::GREAT;

            label sampleI = ppTree.findNearest
            (
                edgePoint,
                tightest,
                tightestDist
            );

            if (sampleI == -1)
            {
                // No point close enough to surface edge.
                break;
            }
            if (tightestDist < maxDist[sampleI])
            {
                nearest.insert(sampleI, surfEdgeI);
            }

            if (exit)
            {
                break;
            }

            // Step to next sample point using local distance.
            // Truncate to max 1/minSampleDist samples per feature edge.
            s += max(minSampleDist*eMag, sampleDist[sampleI]);

            if (s >= (1-minSampleDist)*eMag)
            {
                // Do one more sample, at edge endpoint
                s = eMag;
                exit = true;
            }
        }

        i++;
    }



    if (debug)
    {
        // Dump to obj file

        Info<< "Dumping nearest surface feature edges to nearestEdges.obj\n"
            << "View this Lightwave-OBJ file with e.g. javaview\n" << endl;

        OFstream objStream("nearestEdges.obj");

        label vertI = 0;
        for
        (
            Map<label>::const_iterator iter = nearest.begin();
            iter != nearest.end();
            ++iter
        )
        {
            label sampleI = iter.key();

            meshTools::writeOBJ(objStream, samples[sampleI]); vertI++;

            const edge& e = surfEdges[iter()];

            point nearPt =
                e.line(surfPoints).nearestDist(samples[sampleI]).rawPoint();

            meshTools::writeOBJ(objStream, nearPt); vertI++;

            objStream<< "l " << vertI-1 << ' ' << vertI << endl;
        }
    }

    return nearest;
}


// Get nearest edge on patch for regularly sampled points along the feature
// edges. Return map from patch edge to feature edge.
//
// Q: using point-based sampleDist and maxDist (distance to look around
// each point). Should they be edge-based e.g. by averaging or max()?
Foam::Map<Foam::pointIndexHit> Foam::features::nearestEdgesToFeatEdges
(
    const edgeList& sampleEdges,
    const pointField& samplePoints,
    const scalarField& sampleDist,
    const scalarField& maxDist,
    const scalar minSampleDist
) const
{
    // Build tree out of all sample edges.
    labelList allEdges(sampleEdges.size());
    forAll(allEdges, i)
    {
        allEdges[i] = i;
    }
    octreeDataEdges shapes(sampleEdges, samplePoints, allEdges);

    treeBoundBox overallBb(samplePoints);

    octree<octreeDataEdges> ppTree
    (
        overallBb,  // overall search domain
        shapes,     // all information needed to do checks on cells
        1,          // min levels
        20.0,       // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity
    );

    const pointField& surfPoints = surf_.localPoints();
    const edgeList& surfEdges = surf_.edges();

    scalar maxSearch = max(maxDist);
    vector span(maxSearch, maxSearch, maxSearch);


    Map<pointIndexHit> nearest(2*sampleEdges.size());

    label i = 0;

    //
    // Loop over all feature edges. Sample at regular intervals. Find nearest
    // sampleEdges (using octree)
    //

    for
    (
        labelHashSet::const_iterator iter = edgeLabels_.begin();
        iter != edgeLabels_.end();
        ++iter
    )
    {
        label surfEdgeI = iter.key();

        const edge& e = surfEdges[surfEdgeI];

        if (debug && (i % 1000) == 0)
        {
            Info<< "looking at surface feature edge " << surfEdgeI
                << " verts:" << e << " points:" << surfPoints[e[0]]
                << ' ' << surfPoints[e[1]] << endl;
        }

        // Normalized edge vector
        vector eVec = e.vec(surfPoints);
        scalar eMag = mag(eVec);
        eVec /= eMag;


        //
        // Sample along edge
        //

        bool exit = false;

        // Coordinate along edge (0 .. eMag)
        scalar s = 0.0;

        while (true)
        {
            point edgePoint(surfPoints[e.start()] + s*eVec);

            treeBoundBox tightest(edgePoint - span, edgePoint + span);
            scalar tightestDist = Foam::GREAT;

            label sampleEdgeI = ppTree.findNearest
            (
                edgePoint,
                tightest,
                tightestDist
            );

            if (sampleEdgeI == -1)
            {
                // No edge close enough to surface edge.
                break;
            }

            const edge& e = sampleEdges[sampleEdgeI];

            if (tightestDist < maxDist[e.start()])
            {
                nearest.insert
                (
                    sampleEdgeI,
                    pointIndexHit(true, edgePoint, surfEdgeI)
                );
            }

            if (exit)
            {
                break;
            }

            // Step to next sample point using local distance.
            // Truncate to max 1/minSampleDist samples per feature edge.
//            s += max(minSampleDist*eMag, sampleDist[e.start()]);
s += 0.1*eMag;

            if (s >= (1-minSampleDist)*eMag)
            {
                // Do one more sample, at edge endpoint
                s = eMag;
                exit = true;
            }
        }

        i++;
    }


    if (debug)
    {
        // Dump to obj file

        Info<< "Dumping nearest surface feature edges to nearestEdgeEdges.obj\n"
            << "View this Lightwave-OBJ file with e.g. javaview\n" << endl;

        OFstream objStream("nearestEdgeEdges.obj");

        label vertI = 0;
        for
        (
            Map<pointIndexHit>::const_iterator iter = nearest.begin();
            iter != nearest.end();
            ++iter
        )
        {
            label sampleEdgeI = iter.key();

            const edge& sampleEdge = sampleEdges[sampleEdgeI];

            // Write line from edgeMid to point on feature edge
            const point edgePt =
                0.5
              * (
                    samplePoints[sampleEdge[0]]
                  + samplePoints[sampleEdge[1]]
                );

            meshTools::writeOBJ(objStream, edgePt); vertI++;

            const point& surfPt = iter().rawPoint();

            meshTools::writeOBJ(objStream, surfPt); vertI++;

            objStream<< "l " << vertI-1 << ' ' << vertI << endl;
        }
    }

    return nearest;
}


// Get nearest surface feature edge for every sample. Return in form of
// labelLists giving surfaceEdge label&intersectionpoint.
void Foam::features::nearestSurfFeatures
(
    const pointField& samples,
    const scalarField& maxDist,
    labelList& edgeLabel,
    labelList& edgeEndPoint,
    pointField& edgePoint
) const
{
    edgeLabel.setSize(samples.size());
    edgeEndPoint.setSize(samples.size());
    edgePoint.setSize(samples.size());

    const pointField& points = surf_.localPoints();

    labelList featEdgeLabels(edgeLabels_.toc());

    octreeDataEdges shapes
    (
        surf_.edges(),
        surf_.localPoints(),
        featEdgeLabels
    );

    treeBoundBox overallBb(points);

    octree<octreeDataEdges> ppTree
    (
        overallBb,  // overall search domain
        shapes,     // all information needed to do checks on cells
        1,          // min levels
        20.0,       // maximum ratio of cubes v.s. cells
        100.0       // max. duplicity; n/a since no bounding boxes.
    );


    forAll(samples, i)
    {
        const point& sample = samples[i];

        point maxDistPt(maxDist[i], maxDist[i], maxDist[i]);

        treeBoundBox tightest(sample - maxDistPt, sample + maxDistPt);

        scalar tightestDist = maxDist[i];

        label index =
            ppTree.findNearest
            (
                sample,
                tightest,
                tightestDist
            );


        if (index == -1)
        {
            edgeLabel[i] = -1;
        }
        else
        {
            edgeLabel[i] = featEdgeLabels[index];

            // Unfortunately findNearest does not return nearest point so
            // recalculate
            const edge& e = surf_.edges()[edgeLabel[i]];

            pointIndexHit pHit =
                edgeNearest
                (
                    points[e.start()],
                    points[e.end()],
                    sample
                );

            edgePoint[i] = pHit.rawPoint();
            edgeEndPoint[i] = pHit.index();

        }
    }
}


// ************************************************************************* //
