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

Class
    surfaceSnap

Description
    Takes points on primitivePatch and snaps them to a surface in order of:
    - feature points
    - feature edges
    - nearest point

    Gets points to move as a primitivePatch. Various modes can be set through
    the dictionary (see the code).

    Note: is work in progress!

SourceFiles
    surfaceSnap.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceSnap_H
#define surfaceSnap_H

#include "scalar.H"
#include "pointField.H"
#include "Map.H"
#include "labelHashSet.H"
#include "primitivePatch.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class triSurfaceSearch;
class features;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class surfaceSnap Declaration
\*---------------------------------------------------------------------------*/

class surfaceSnap
{
    // Private data

        //- Reference to patch whose points to move
        const primitivePatch& pp_;

        //- Reference to surface
        const triSurfaceSearch& querySurf_;

        //- Reference to feature points and edges.
        const features& surfFeatures_;

        //- Relaxation factor (0:take point only; 1:take average of neighbours)
        scalar relax_;

        //- Max distance for points to snap to feature
        scalar snapTol_;

        //- Max distance for points to project on to surface
        scalar projectTol_;

        //- Whether ok to create illegal faces when snapping to features
        const bool forceSnap_;

        //- Current patch points.
        pointField points_;


    // Private Member Functions

        // Calculate minimum length of edges on patch
        scalarField minEdgeLength(const pointField&) const;

        //- Calculate average of all neighbours of pointI.
        //  If allPoints=true use all neighbours, if false use only
        //  non-feature points.
        point getAvgPos
        (
            const Map<label>& patchToFeatureEdge,
            const Map<label>& patchToFeaturePoint,
            const bool allPoints,
            const label pointI
        ) const;

        ////- Check if snapping points of face would result in zero sized face
        ////  with current point set. If so remove one face point from 
        ////  patchToFeatureEdge.
        //void unsnapFace
        //(
        //    const label faceI,
        //    Map<label>& patchToFeatureEdge
        //) const;

        //- Calculate snapped position of patchPointI. Limits distance for
        //  projection (not feature snapping) to projectDist.
        point snapPosition
        (
            const Map<label>& patchToFeatureEdge,
            const Map<label>& patchToFeaturePoint,
            const label patchPointI,
            const scalar projectDist
        ) const;

        //- Return true if faces using point have valid areas with current
        //  points_.
        bool validAreas(const label patchPointI) const;

        //- Snap point onto surface feature or nearest. Limits distance
        //  projected such that no zero sized faces are constructed.
        //  Modifies points_.
        bool snapPoint
        (
            const label patchPointI,
            const scalar projectDist,
            const point& wantedPos
        );

        //- Find unsnapped patchPoint that is closest to a feature edge.
        //  Feature edges in
        //  - surfEdges[patchPointI] (-1 or surfEdgeI)
        //  - surfEdgePoints[patchPointI] (nearest position on surfEdgeI)
        label closestToFeatEdge
        (
            const labelList& surfEdges,
            const pointField& surfEdgePoints,
            const boolList& snappedPoint
        ) const;

        //- Determine which edge connected to fromPointI makes best angle
        //  with surface feature edge.
        label maxAlignedEdge
        (
            const labelList& surfEdges,
            const boolList& snappedPoint,
            const label fromEdgeI,
            const label fromPointI
        ) const;

        //- Snap all unsnapped points that are close to feature edges onto
        //  these feature edges. It tries to stitch whole lines of connected
        //  mesh edges.
        void stitchToFeatureEdges
        (
            const scalarField& snapDist,
            boolList& snappedPoint
        );

        //- Snap all unsnapped points that are close to feature edges onto
        //  these feature edges. It tries to stitch whole lines of connected
        //  mesh edges.
        void snapToFeatureEdges
        (
            const scalarField& snapDist,
            const Map<label>& patchToFeatureEdge,
            boolList& snappedPoint
        );


        //- Snap remaining feature edges. If forceSnap even snap points creating
        //  illegal faces.
        label snapFeaturesPost
        (
            const scalarField& snapDist,
            const Map<label>& patchToFeatureEdge,
            boolList& snappedPoint
        );


        //- Disallow default bitwise copy construct
        surfaceSnap(const surfaceSnap&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceSnap&);


public:

        ClassName("surfaceSnap");

    // Constructors

        //- Construct from patch and surface and do one snap iteration.
        surfaceSnap
        (
            const primitivePatch& pp,
            const triSurfaceSearch& querySurf,
            const features& surfFeatures,
            const dictionary& dict
        );


    // Member Functions

        // Access

            scalar relax() const
            {
                return relax_;
            }

            scalar& relax()
            {
                return relax_;
            }

            scalar snapTol() const
            {
                return snapTol_;
            }

            scalar& snapTol()
            {
                return snapTol_;
            }

            scalar projectTol() const
            {
                return projectTol_;
            }

            scalar& projectTol()
            {
                return projectTol_;
            }

            const pointField& points() const
            {
                return points_;
            }

            pointField& points()
            {
                return points_;
            }


        // Edit

            // Move all of points_ on to surface. Return true if successful.
            bool snap();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
