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

#include "orientedSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "treeBoundBox.H"
#include "octree.H"
#include "octreeDataTriSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::orientedSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Return true if face uses edge from start to end.
bool Foam::orientedSurface::edgeOrder
(
    const labelledTri& f,
    const edge& e
)
{
    if
    (
        (f[0] == e[0] && f[1] == e[1])
     || (f[1] == e[0] && f[2] == e[1])
     || (f[2] == e[0] && f[0] == e[1])
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Return true if edge is used in opposite order in faces
bool Foam::orientedSurface::consistentEdge
(
    const edge& e,
    const labelledTri& f0,
    const labelledTri& f1
)
{
    return edgeOrder(f0, e) ^ edgeOrder(f1, e);
}


Foam::labelList Foam::orientedSurface::faceToEdge
(
    const labelList& changedFaces
) const
{
    labelList changedEdges(3*changedFaces.size());
    label changedI = 0;

    forAll(changedFaces, i)
    {
        const labelList& fEdges = faceEdges()[changedFaces[i]];

        forAll(fEdges, j)
        {
            changedEdges[changedI++] = fEdges[j];
        }
    }
    changedEdges.setSize(changedI);

    return changedEdges;
}


Foam::labelList Foam::orientedSurface::edgeToFace
(
    const labelList& changedEdges,
    labelList& flip
) const
{
    labelList changedFaces(2*changedEdges.size());
    label changedI = 0;

    forAll(changedEdges, i)
    {
        label edgeI = changedEdges[i];

        const labelList& eFaces = edgeFaces()[edgeI];

        if (eFaces.size() < 2)
        {
            // Do nothing, faces was already visited.
        }
        else if (eFaces.size() == 2)
        {
            label face0 = eFaces[0];
            label face1 = eFaces[1];

            const labelledTri& f0 = operator[](face0);
            const labelledTri& f1 = operator[](face1);

            if (flip[face0] == UNVISITED)
            {
                if (flip[face1] == UNVISITED)
                {
                    FatalErrorIn("orientedSurface::edgeToFace") << "Problem"
                        << abort(FatalError);
                }
                else
                {
                    // Face1 has a flip state, face0 hasn't
                    if (consistentEdge(edges()[edgeI], f0, f1))
                    {
                        // Take over flip status
                        flip[face0] = (flip[face1] == FLIP ? FLIP : NOFLIP);
                    }
                    else
                    {
                        // Invert
                        flip[face0] = (flip[face1] == FLIP ? NOFLIP : FLIP);
                    }
                    changedFaces[changedI++] = face0;
                }
            }
            else
            {
                if (flip[face1] == UNVISITED)
                {
                    // Face0 has a flip state, face1 hasn't
                    if (consistentEdge(edges()[edgeI], f0, f1))
                    {
                        flip[face1] = (flip[face0] == FLIP ? FLIP : NOFLIP);
                    }
                    else
                    {
                        flip[face1] = (flip[face0] == FLIP ? NOFLIP : FLIP);
                    }
                    changedFaces[changedI++] = face1;
                }
            }
        }
        else
        {
            // Multiply connected. Do what?
        }
    }
    changedFaces.setSize(changedI);

    return changedFaces;
}


void Foam::orientedSurface::propagateOrientation
(
    const point& samplePoint,
    const bool orientOutside,
    const label nearestFaceI,
    const point& nearestPt,
    labelList& flipState
) const
{
    //
    // Determine orientation to normal on nearest face
    //

    vector n =
        triSurfaceTools::surfaceNormal(*this, nearestFaceI, nearestPt);

    if (debug)
    {
        Pout<< "orientedSurface::propagateOrientation : starting face"
            << " orientation:" << nl
            << "     for samplePoint:" << samplePoint << nl
            << "     starting from point:" << nearestPt << nl
            << "     on face:" << nearestFaceI << nl
            << "     with normal:" << n << endl;
    }

    // Check orientation of samplePt.
    bool sampleOutside = (((samplePoint - nearestPt) & n) > 0);

    if (sampleOutside == orientOutside)
    {
        // Normals on surface outwards pointing. No need to flip normals
        flipState[nearestFaceI] = NOFLIP;
    }
    else
    {
        // Need to flip normals.
        flipState[nearestFaceI] = FLIP;
    }


    // List of faces that were changed in the last iteration.
    labelList changedFaces(1, nearestFaceI);
    // List of edges that were changed in the last iteration.
    labelList changedEdges;
    
    while(true)
    {
        changedEdges = faceToEdge(changedFaces);

        if (debug)
        {
            Pout<< "From changedFaces:" << changedFaces.size()
                << " to changedEdges:" << changedEdges.size()
                << endl;
        }

        if (changedEdges.size() == 0)
        {
            break;
        }

        changedFaces = edgeToFace(changedEdges, flipState);

        if (debug)
        {
            Pout<< "From changedEdges:" << changedEdges.size()
                << " to changedFaces:" << changedFaces.size()
                << endl;
        }

        if (changedFaces.size() == 0)
        {
            break;
        }
    }
}


void Foam::orientedSurface::flipSurface(const labelList& flipState)
{
    // Flip tris in *this
    forAll(flipState, faceI)
    {
        if (flipState[faceI] == UNVISITED)
        {
            FatalErrorIn
            (
                "orientSurface(const point&, const label, const point&)"
            )   << "unvisited face " << faceI
                << abort(FatalError);
        }
        else if (flipState[faceI] == FLIP)
        {
            labelledTri& tri = operator[](faceI);

            label tmp = tri[0];

            tri[0] = tri[1];
            tri[1] = tmp;
        }
    }
    // Recalculate normals
    clearOut();
}


void Foam::orientedSurface::orientSurface
(
    const point& samplePoint,
    const bool orientOutside
)
{
    // Whether face has to be flipped.
    //      UNVISITED: unvisited
    //      NOFLIP: no need to flip
    //      FLIP: need to flip
    labelList flipState(size(), UNVISITED);


    while (true)
    {
        // Linear search for nearest unvisited point on surface.

        scalar minDist = GREAT;
        point minPoint;
        label minFaceI = -1;

        forAll(*this, faceI)
        {
            if (flipState[faceI] == UNVISITED)
            {
                const labelledTri& f = operator[](faceI);

                pointHit curHit =
                    triPointRef
                    (
                        points()[f[0]],
                        points()[f[1]],
                        points()[f[2]]
                    ).nearestPoint(samplePoint);

                if (curHit.distance() < minDist)
                {
                    minDist = curHit.distance();
                    minPoint = curHit.rawPoint();
                    minFaceI = faceI;
                }
            }
        }

        // Did we find anything?
        if (minFaceI == -1)
        {
            break;
        }

        // From this nearest face see if needs to be flipped and then
        // go outwards.
        propagateOrientation
        (
            samplePoint,
            orientOutside,
            minFaceI,
            minPoint,
            flipState
        );
    }

    // Now finally flip triangles according to flipState.
    flipSurface(flipState);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::orientedSurface::orientedSurface()
:
    triSurface()
{}


// Construct from surface and point which defines outside
Foam::orientedSurface::orientedSurface
(
    const triSurface& surf,
    const point& samplePoint,
    const bool orientOutside
)
:
    triSurface(surf)
{
    orientSurface(samplePoint, orientOutside);
}


// Construct from surface. Calculate outside point.
Foam::orientedSurface::orientedSurface
(
    const triSurface& surf,
    const bool orientOutside
)
:
    triSurface(surf)
{
    treeBoundBox bb(localPoints());

    point outsidePoint = 2 * bb.max() - bb.min();

    orientSurface(outsidePoint, orientOutside);
}


// ************************************************************************* //
