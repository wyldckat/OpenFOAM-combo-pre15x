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
    Routines to check for point-nearness.

\*---------------------------------------------------------------------------*/

#include "checkCoords.H"
#include "pointSet.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static const scalar edgeRelTol_ = 0.1;

// Min of (length of) connected edges.
scalar minEdgeLen(const primitiveMesh& mesh, const label vertI)
{
    const edgeList& edges = mesh.edges();

    const pointField& points = mesh.points();

    const labelList& pEdges = mesh.pointEdges()[vertI];

    scalar minLen = GREAT;

    forAll(pEdges, pEdgeI)
    {
        minLen = min(minLen, edges[pEdges[pEdgeI]].mag(points));
    }
    return minLen;
}


// Checks if points are close. Uses edge-relative tolerance
bool closePoints
(
    const primitiveMesh& mesh,
    const scalar tol,
    const label point0,
    const label point1
)
{
    const point& p0 = mesh.points()[point0];
    const point& p1 = mesh.points()[point1];

    // Determine minimum edge.
    scalar minLen = min(minEdgeLen(mesh, point0), minEdgeLen(mesh, point1));

    if (mag(p0 - p1) < tol*minLen)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Merge points and create mapping array. Return mapping.
bool checkCoords
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar mergeTol,
    labelHashSet* setPtr
)
{
    const pointField& pts = mesh.points();

    // Sort points
    SortableList<scalar> sortedMag(mag(pts));

    label nClose = 0;

    for (label i = 1; i < sortedMag.size(); i++)
    {
        label ptI = sortedMag.indices()[i];

        label prevPtI = sortedMag.indices()[i-1];

        if (closePoints(mesh, edgeRelTol_, ptI, prevPtI))
        {
            nClose++;

            if (report)
            {
                Pout<< "checkCoords : points "
                    << ptI << " and " << prevPtI
                    << " with coordinates:" << pts[ptI]
                    << " and " << pts[prevPtI] << endl
                    << " are relatively close to each other."
                    << " This might be correct in case of e.g. baffles."
                    << endl;
            }

            if (setPtr)
            {
                setPtr->insert(ptI);
                setPtr->insert(prevPtI);
            }
        }
    }

    reduce(nClose, sumOp<label>());

    if (nClose > 0)
    {
        if (Pstream::master())
        {
            WarningIn
            (
                "bool checkCoords(const primitiveMesh& mesh, "
                "const bool report, const scalar mergeTol, "
                "labelHashSet* setPtr)"
            ) << nClose  << " points that are close together found." << nl
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
