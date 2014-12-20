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

Class
    features

Description
    Surface features calculation and handling.

SourceFiles
    features.C

\*---------------------------------------------------------------------------*/

#ifndef features_H
#define features_H

#include "scalar.H"
#include "pointField.H"
#include "Map.H"
#include "labelHashSet.H"
#include "pointIndexHit.H"
#include "boolList.H"
#include "edgeList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class features Declaration
\*---------------------------------------------------------------------------*/

class features
{
    // Private data

        //- Reference to surface
        const triSurface& surf_;

        //- Feature angle for points
        const scalar pointMinCos_;

        //- Feature angle for edges.
        const scalar edgeMinCos_;

        //- Min lenth for feature edges to be recognized as such
        const scalar minDim_;

        //- Labels of points that are feature points
        labelHashSet pointLabels_;

        //- Labels of edges that are feature edges.
        labelHashSet edgeLabels_;



    // Private Member Functions

        //- Return nearest point on edge (start..end). Also classify nearest:
        //  index=-1: nearest on mid of edge. index=0:nearest on edge.start()
        //  index=1: nearest on edge.end().
        static pointIndexHit edgeNearest
        (
            const point& start,
            const point& end,
            const point& sample
        );

        //- Choose next unset feature edge.
        label nextFeatEdge
        (
            const boolList& isFeatEdge,
            const labelList& featVisited,
            const label unsetVal,
            const label prevEdgeI,
            const label vertI
        ) const;

        //- Walk connected feature edges. Marks edges in featVisited.
        scalar walkSegment
        (
            const bool mark,
            const boolList& isFeatEdge,
            const label startEdgeI,
            const label startPointI,
            const label currentFeatI,     
            labelList& featVisited
        );

        //- Remove edges < minDim cumulative length
        void filterFeatureEdges();

        //- Do all calculations. Sets pointLabels_, edgeLabels_
        void calculateFeatures();


        //- Disallow default bitwise copy construct
        features(const features&);

        //- Disallow default bitwise assignment
        void operator=(const features&);


public:

    ClassName("features");

    // Constructors

        //- Construct from surface, cos of angle and minimum length.
        //  Set cos to <-1 to have no features
        features(const triSurface&, const scalar minCos);

        //- Construct from surface, separate angles for pointFeatures and
        //  edge features and minimum length
        //  Set cos to <-1 to have no features
        features
        (
            const triSurface&,
            const scalar pointMinCos,
            const scalar edgeMinCos,
            const scalar minDim = 0
        );


    // Member Functions

        // Access

            //- Labels of points on surface that are feature points
            const labelHashSet& pointLabels() const
            {
                return pointLabels_;
            }

            //- Labels of edges on surface that are feature edges
            const labelHashSet& edgeLabels() const
            {
                return edgeLabels_;
            }


        // Write

            //- Dump feature points and edges to .obj file for visualization
            void writeOBJ() const;


        // Find

            //- Find nearest sample for every point (usually the set of
            //  featurePoints). Return map from index in
            //  samples to surface point. Do not include points that are further
            //  than maxDist away (separate maxDist for every sample)
            Map<label> nearestSamples
            (
                const labelHashSet& points,
                const pointField& samples,
                const scalarField& maxDist
            ) const;

            //- Find nearest sample for regularly sampled points along the
            //  feature edges. Return map from sample to feature edge.
            //  maxDist is distance below which gets snapped.
            //  Edge gets sampled at points sampleDist[sampleI] apart.
            //  (with a maximum of 10 samples per edge)
            Map<label> nearestSamplesToFeatEdges
            (
                const pointField& samples,
                const scalarField& sampleDist,
                const scalarField& maxDist,
                const scalar minSampleDist = 0.1
            ) const;

            //- Like nearestSamplesToFeatEdges but now gets nearest point on
            //  sample-edge instead of nearest sample-point itself.
            //  Return map from patch edge to feature edge and map from patch
            Map<pointIndexHit> nearestEdgesToFeatEdges
            (
                const edgeList& sampleEdges,
                const pointField& samplePoints,
                const scalarField& sampleDist,
                const scalarField& maxDist,
                const scalar minSampleDist = 0.1
            ) const;


            //- Find nearest surface feature edge for each sample.
            //  Sets:
            //  - edgeLabel : label of surface edge.
            //  - edgePoint : exact position of nearest point on edge.
            //  - edgeEndPoint : -1, 0, 1 depending on whether nearest is
            //                  on inside/start/end of edge
            void nearestSurfFeatures
            (
                const pointField& samples,
                const scalarField& maxDist,
                labelList& edgeLabel,
                labelList& edgeEndPoint,
                pointField& edgePoint
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
