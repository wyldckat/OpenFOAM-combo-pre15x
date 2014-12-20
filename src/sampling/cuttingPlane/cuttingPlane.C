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
    Creates a primitivePatch from a plane and a mesh

\*---------------------------------------------------------------------------*/

#include "cuttingPlane.H"
#include "primitiveMesh.H"
#include "linePointRef.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Find cut cells
void cuttingPlane::calcCutCells
(
    const primitiveMesh& mesh,
    const scalarField& dotProducts
)
{
    const labelListList& cellEdges = mesh.cellEdges();
    const edgeList& edges = mesh.edges();

    cutCells_.setSize(cellEdges.size());
    label cutcellI(0);

    // Find the cut cells by detecting any cell that uses points with
    // opposing dotProducts.
    forAll(cellEdges, cellI)
    {
        const labelList& cEdges = cellEdges[cellI];

        label nCutEdges = 0;

        forAll(cEdges, i)
        {
            const edge& e = edges[cEdges[i]];

            if
            (
                (dotProducts[e[0]] < 0 && dotProducts[e[1]] > 0)
             || (dotProducts[e[1]] < 0 && dotProducts[e[0]] > 0)
            )
            {
                nCutEdges++;

                if (nCutEdges > 2)
                {
                    cutCells_[cutcellI++] = cellI;

                    break;
                }
            }
        }
    }

    // Set correct list size
    cutCells_.setSize(cutcellI);
}


// Determine for each edge the intersection point. Calculates
// - cutPoints_ : coordinates of all intersection points
// - edgePoint  : per edge -1 or the index into cutPoints_
labelList cuttingPlane::intersectEdges
(
    const primitiveMesh& mesh,
    const scalarField& dotProducts
)
{
    // Use the dotProducts to find out the cut edges.
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    // Per edge -1 or the label of the intersection point
    labelList edgePoint(edges.size(), -1);

    DynamicList<point> dynCuttingPoints(4*cutCells_.size());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if
        (
            (dotProducts[e[0]] < 0 && dotProducts[e[1]] > 0)
         || (dotProducts[e[1]] < 0 && dotProducts[e[0]] > 0)
        )
        {
            // Edge is cut.
            const point& p0 = points[e[0]];
            const point& p1 = points[e[1]];

            scalar alpha = lineIntersect(linePointRef(p0, p1));

            dynCuttingPoints.append((1-alpha)*p0 + alpha*p1);
            edgePoint[edgeI] = dynCuttingPoints.size() - 1;
        }
    }

    dynCuttingPoints.shrink();
    cutPoints_.transfer(dynCuttingPoints);
    dynCuttingPoints.clear();

    return edgePoint;
}


// Coming from startEdgeI cross the edge to the other face
// across to the next cut edge.
bool cuttingPlane::walkCell
(
    const primitiveMesh& mesh,
    const labelList& edgePoint,
    const label cellI,
    const label startEdgeI,
    DynamicList<label>& faceVerts
)
{
    label faceI = -1;
    label edgeI = startEdgeI;

    label nIter = 0;

    do
    {
        faceVerts.append(edgePoint[edgeI]);

        // Cross edge to other face
        faceI = meshTools::otherFace(mesh, cellI, faceI, edgeI);

        // Find next cut edge on face.
        const labelList& fEdges = mesh.faceEdges()[faceI];

        label nextEdgeI = -1;

        //Note: here is where we should check for whether there are more
        // than 2 intersections with the face (warped/non-convex face).
        // If so should e.g. decompose the cells on both faces and redo
        // the calculation.

        forAll(fEdges, i)
        {
            label edge2I = fEdges[i];

            if (edge2I != edgeI && edgePoint[edge2I] != -1)
            {
                nextEdgeI = edge2I;
                break;
            }
        }

        if (nextEdgeI == -1)
        {
            // Did not find another cut edge on faceI. Do what?
            WarningIn("cuttingPlane::walkCell")
                << "Did not find closed walk along surface of cell " << cellI
                << " starting from edge " << startEdgeI
                << " in " << nIter << " iterations." << nl
                << "Collected cutPoints so far:" << faceVerts
                << endl;

            return false;
        }

        edgeI = nextEdgeI;

        nIter++;

        if (nIter > 1000)
        {
            WarningIn("cuttingPlane::walkCell")
                << "Did not find closed walk along surface of cell " << cellI
                << " starting from edge " << startEdgeI
                << " in " << nIter << " iterations." << nl
                << "Collected cutPoints so far:" << faceVerts
                << endl;
            return false;
        }

    } while (edgeI != startEdgeI);


    if (faceVerts.size() >= 3)
    {
        return true;
    }
    else
    {
        WarningIn("cuttingPlane::walkCell")
            << "Did not find closed walk along surface of cell " << cellI
            << " starting from edge " << startEdgeI << nl
            << "Collected cutPoints so far:" << faceVerts
            << endl;

        return false;
    }
}


// For every cut cell determine a walk through all? its cuts.
void cuttingPlane::walkCellCuts
(
    const primitiveMesh& mesh,
    const labelList& edgePoint
)
{
    cutFaces_.setSize(cutCells_.size());
    label cutFaceI = 0;

    forAll(cutCells_, i)
    {
        label cellI = cutCells_[i];

        // Find the starting edge to walk from.
        const labelList& cEdges = mesh.cellEdges()[cellI];

        label startEdgeI = -1;

        forAll(cEdges, cEdgeI)
        {
            label edgeI = cEdges[cEdgeI];

            if (edgePoint[edgeI] != -1)
            {
                startEdgeI = edgeI;
                break;
            }
        }

        // Check for the unexpected ...
        if (startEdgeI == -1)
        {
            FatalErrorIn("cuttingPlane::walkCellCuts") << "Cannot find cut edge"
                << " for cut cell " << cellI << abort(FatalError);
        }

        // Walk from starting edge around the circumference of the cell.
        DynamicList<label> faceVerts(2*mesh.faces()[cellI].size());
        bool okCut = walkCell
        (
            mesh,
            edgePoint,
            cellI,
            startEdgeI,
            faceVerts
        );

        if (okCut)
        {
            faceVerts.shrink();

            face cutFace(faceVerts);

            // Orient face.
            if ((cutFace.normal(cutPoints_) && normal()) < 0)
            {
                cutFace = cutFace.reverseFace();
            }

            cutFaces_[cutFaceI++] = cutFace;
        }
    }

    cutFaces_.setSize(cutFaceI);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cuttingPlane::cuttingPlane(const primitiveMesh& mesh, const plane& newPlane)
:
    plane(newPlane)
{
    scalarField dotProducts = (mesh.points() - refPoint()) & normal();

    //// Perturb points cuts so edges are cut properly.
    //const tmp<scalarField> tdotProducts = stabilise(rawDotProducts, SMALL);
    //const scalarField& dotProducts = tdotProducts();

    // Determine cells that are (probably) cut.
    calcCutCells(mesh, dotProducts);

    // Determine cutPoints and return list of edge cuts. (per edge -1 or
    // the label of the intersection point (in cutPoints_)
    labelList edgePoint(intersectEdges(mesh, dotProducts));

    // Do topological walk around cell to find closed loop. 
    walkCellCuts(mesh, edgePoint);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return vectorField of cutting points
const pointField& cuttingPlane::points() const
{
    return cutPoints_;
}


// Return unallocFaceList of points in cells
const faceList& cuttingPlane::faces() const
{
    return cutFaces_;
}


// Return labelList of cut cells
const labelList& cuttingPlane::cells() const
{
    return cutCells_;
}


bool cuttingPlane::cut()
{
    if(cutCells_.size() > 0)
    {
        return true;
    }
    else 
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void cuttingPlane::operator=(const cuttingPlane& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("cuttingPlane::operator=(const cuttingPlane&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    cutCells_ = rhs.cells();
    cutPoints_ = rhs.points();
    cutFaces_ = rhs.faces();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
