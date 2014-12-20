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
    orientedSurface

Description
    Given outside point flip all faces such that normals point outwards.

SourceFiles
    orientedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef orientedSurface_H
#define orientedSurface_H

#include "triSurface.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class orientedSurface Declaration
\*---------------------------------------------------------------------------*/

class orientedSurface
:
    public triSurface
{
    // Data types

        //- Enumeration listing whether face needs to be flipped.
        enum sideStat
        {
            UNVISITED,
            FLIP,
            NOFLIP
        };


    // Private Member Functions

        //- Return true if face uses edge from start to end.
        static bool edgeOrder(const labelledTri&, const edge&);

        //- Return true if edge is used in opposite order in faces
        static bool consistentEdge
        (
            const edge& e,
            const labelledTri& f0,
            const labelledTri& f1
        );

        //- From changed faces get the changed edges
        labelList faceToEdge(const labelList& changedFaces) const;

        //- From changed edges check the orientation of the connected faces
        //  and flip them. Return changed faces.
        labelList edgeToFace
        (
            const labelList& changedEdges,
            labelList& flip
        ) const;


        static label findEdge(const edgeList&, const labelList&, const edge&);

        //- Given nearest point and face check orientation to nearest face
        //  and flip if nessecary (only marked in flipState) and propagate.
        void propagateOrientation
        (
            const point& outsidePoint,
            const label nearestFaceI,
            const point& nearestPt,
            labelList& flipState
        ) const;

        //- Given flipState reverse triangles of *this
        void flipSurface(const labelList& flipState);

        //- Flip faces such that normals are consistent with outsidePoint
        void orientSurface(const point& outsidePoint);

public:

    ClassName("orientedSurface");


    // Constructors

        //- Construct null
        orientedSurface();

        //- Construct from triSurface and outside point. Uses linear search
        //  to find nearest.
        orientedSurface(const triSurface&, const point& outsidePoint);

        //- Construct from triSurface. Calculates outside point as being
        //  outside the boundaing box of the surface.
        orientedSurface(const triSurface&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
