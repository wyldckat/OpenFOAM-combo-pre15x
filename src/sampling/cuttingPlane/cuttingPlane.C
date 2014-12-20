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
    Creates a primitivePatch from a plane and a mesh

\*---------------------------------------------------------------------------*/

#include "cuttingPlane.H"
#include "primitiveMesh.H"
#include "linePointRef.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Find cut cells
void cuttingPlane::cutCells
(
    const scalarField& dotProducts, 
    const cellList& CellFaces,
    const faceList& Faces
)
{
    cutCells_.setSize(CellFaces.size());
    label cutcellI(0);

    // Find the cut cells
    forAll(CellFaces, cellI)
    {
        labelList cellPoints = CellFaces[cellI].labels(Faces);
        bool plusPlane = false;
        bool minusPlane = false;

        forAll(cellPoints, pointi)
        {
            // Cuts included faces of negative dotProducts cells
            if(dotProducts[cellPoints[pointi]] >= 0)
            {
                plusPlane = true;
            }
            else
            {
                minusPlane = true;
            }

            //If cell is cut
            if(plusPlane && minusPlane)
            {
                cutCells_[cutcellI] = cellI;
                cutcellI++;

                break;
            }
        }
    }

    // Set correct list size
    cutCells_.setSize(cutcellI);
}

labelListList cuttingPlane::cuttingPoints
(
    const scalarField& dotProducts,
    const labelListList& CellEdges,
    const pointField& Points
)
{
    // Edge field references
    const edgeList& Edges = mesh_.edges();

    // Checklist to make sure points and edges are recorded only once 
    // if they are cut
    labelListList intersections(2);

    {
        labelList edgeIntersections(Edges.size(), -1);
        intersections[0] = edgeIntersections;
        labelList pointIntersections(Points.size(), -1);
        intersections[1] = pointIntersections;
    }

    // Set private data list sizes
    cuttingPoints_.setSize(12*cutCells_.size());
    cutCellCuttingPoints_.setSize(cutCells_.size());

    // Extra loop labels
    label cuttingPointI = 0;
    label cellCuttingPointI = 0;

    // Make a pointList and set size for cutCellCuttingPoint labelListList 
    // of cutting points
    forAll(cutCells_, cellI)
    {
        // Determine number of cutting points and record those points
        cellCuttingPointI = 0;
        forAll(CellEdges[cutCells_[cellI]], edgeI)
        {
            if(intersections[0][CellEdges[cutCells_[cellI]][edgeI]] == -1)
            {
                labelList EdgePoints(2);
                EdgePoints[0] =
                    Edges[CellEdges[cutCells_[cellI]][edgeI]].start();

                EdgePoints[1] =
                    Edges[CellEdges[cutCells_[cellI]][edgeI]].end();
    

                if
                (
                    (
                        dotProducts[EdgePoints[0]] > 0
                     && dotProducts[EdgePoints[1]] < 0
                    )
                 || (
                        dotProducts[EdgePoints[0]] < 0
                     && dotProducts[EdgePoints[1]] > 0
                    )
                 || (
                        dotProducts[EdgePoints[0]] == 0
                     && dotProducts[EdgePoints[1]] != 0
                    )
                 || (
                     dotProducts[EdgePoints[0]] != 0
                  && dotProducts[EdgePoints[1]] == 0
                    )
                )
                {
                   
                    scalar alpha =
                        lineIntersect
                        (
                            linePointRef
                            (
                                Points[EdgePoints[0]],
                                Points[EdgePoints[1]]
                            )
                        );

                    // Check for cutting point proximity to vertex
                    if(mag(Points[EdgePoints[0]] - Points[EdgePoints[1]]) == 0)
                    {
                        FatalErrorIn
                        (
                            "labelListList cuttingPlane::cuttingPoints\n"
                            "(\n"
                                "const scalarField& dotProducts,\n"
                                "const labelListList& CellEdges,\n"
                                "const pointIOField& Points\n"
                            ")\n"
                        )   << "Zero edge length."
                            << abort(FatalError);
                    }

                    scalar CutToPnt0 = mag(alpha);
                    scalar CutToPnt1 = mag(1 - alpha);

                    scalar cuttingPointTolerance = SMALL;
    
                    if
                    (
                        CutToPnt0 < cuttingPointTolerance
                     || CutToPnt1 < cuttingPointTolerance
                    )
                    {
                        if(CutToPnt0 < cuttingPointTolerance)
                        {
                            if(intersections[1][EdgePoints[0]] == -1)
                            {
                                intersections[1][EdgePoints[0]] = cuttingPointI;
                                cuttingPoints_[cuttingPointI] =
                                    Points[EdgePoints[0]];

                                cellCuttingPointI++;
                                cuttingPointI++;
                            }
                            else
                            {
                                cellCuttingPointI++;
                            }
                        }
                        else
                        {
                            if(intersections[1][EdgePoints[1]] == -1)
                            {
                                intersections[1][EdgePoints[1]] = cuttingPointI;
                                cuttingPoints_[cuttingPointI] =
                                    Points[EdgePoints[1]];

                                cellCuttingPointI++;
                                cuttingPointI++;
                            }
                            else
                            {
                                cellCuttingPointI++;
                            }
                        }
                    }
                    else
                    {
                        cuttingPoints_[cuttingPointI] = 
                            Points[EdgePoints[0]]
                          + alpha*
                            (
                                Points[EdgePoints[1]]
                              - Points[EdgePoints[0]]
                            );

                        intersections[0][CellEdges[cutCells_[cellI]][edgeI]] =
                            cuttingPointI;

                        cellCuttingPointI++;
                        cuttingPointI++;
                    }
                }
            }
            else
            {
                cellCuttingPointI++;
            }
        }

        // Set the labelListList size for current cell
        cutCellCuttingPoints_[cellI].setSize(cellCuttingPointI);
    }

    // Refine pointList size
    cuttingPoints_.setSize(cuttingPointI);
    
    return intersections;
}

void cuttingPlane::cutCellCuttingPoints
(
    const labelListList intersections,
    const labelListList& CellEdges,
    const cellList& CellFaces,
    const faceList& Faces
)
{
    // Labels for deleting zero area faces
    label toCell = 0;
    label cellSkip = 0;
    label cellCuttingPointI = 0;

    // Make list of labels for cutting points for each cell 
    // and delete cells with fewer than 3 cutting points. !!!
    forAll(cutCells_, cellI)
    {
        toCell = cellI - cellSkip;

        // Resize toCell for writing list
        if(cellSkip)
        {
            cutCellCuttingPoints_[toCell].setSize
                (cutCellCuttingPoints_[cellI].size());
        }

        cellCuttingPointI = 0;

        // Add cutEdges to cutCellCuttingPoints
        forAll(CellEdges[cutCells_[cellI]], edgeI)
        {
            if(intersections[0][CellEdges[cutCells_[cellI]][edgeI]] != -1)
            {
                cutCellCuttingPoints_[toCell][cellCuttingPointI] =
                    intersections[0][CellEdges[cutCells_[cellI]][edgeI]];
                
                cellCuttingPointI++;
            }
        }

        // Add cutPoints to cutCellCuttingPoints
        const labelList cellLabels =
            CellFaces[cutCells_[cellI]].labels(Faces);

        forAll(cellLabels, pointI)
        {
            if
            (
                intersections[1][cellLabels[pointI]] != -1
            )
            {
                cutCellCuttingPoints_[toCell][cellCuttingPointI] =
                    intersections[1][cellLabels[pointI]];
                
                cellCuttingPointI++;
            }
        }

        cutCells_[toCell] = cutCells_[cellI];

        // Resize cell intersection size
        if(cellCuttingPointI < cutCellCuttingPoints_[toCell].size())
        {
            cutCellCuttingPoints_[toCell].setSize(cellCuttingPointI);
        }

        
        // Change read-write offset
        if(cutCellCuttingPoints_[toCell].size() < 3)
        {
            cellSkip++;
        }
    }

    // Resize transfer vars
    cutCells_.setSize(cutCells_.size() - cellSkip);
    cutCellCuttingPoints_.setSize(cutCells_.size());
}

void cuttingPlane::orderFacePoints()
{
    cutFaces_.setSize(cutCells_.size());

    forAll(cutCells_, cellI)
    {

        const labelList& cellPoints(cutCellCuttingPoints_[cellI]);

        // Calculate face centroid
        vector faceCentre(0, 0, 0);

        forAll(cellPoints, pointI)
        {
            faceCentre += cuttingPoints_[cellPoints[pointI]];
        }

        faceCentre /= cellPoints.size();

        // Calculate right hand angle between point-to-centroid lines
        // and baseLine 
        vector baseLine(cuttingPoints_[cellPoints[0]] - faceCentre);
        slList includedAngle(cellPoints.size());
        includedAngle[0].s_ = 0;

        forAll(cellPoints, pointI)
        {   
            includedAngle[pointI].l_ = cellPoints[pointI];

            if(pointI > 0)
            {
                vector lineCto
                (
                    cuttingPoints_[cellPoints[pointI]] - faceCentre
                );
            
                vector crossProduct = (baseLine ^ lineCto);
                scalar dotProduct = (baseLine & lineCto);
                scalar magProduct = mag(baseLine)*mag(lineCto);

                scalar angle = acos(dotProduct/(magProduct + SMALL));

                if(mag(crossProduct) < VSMALL)
                {
                    includedAngle[pointI].s_ = physicalConstant::pi;
                }
                else
                {
                    scalar alignment =
                        mag(crossProduct/mag(crossProduct) - normal());
                
                    if (alignment < 1)
                    {
                        includedAngle[pointI].s_ = angle;
                    }  
                    else
                    {
                        includedAngle[pointI].s_ =
                            2*physicalConstant::pi - angle;  
                    }
                }
            }
        }      

        sort(includedAngle);

        forAll(cellPoints, pointI)
        {
            cutCellCuttingPoints_[cellI][pointI] = includedAngle[pointI].l_; 
        }

        cutFaces_[cellI] = face(cutCellCuttingPoints_[cellI]);
    }
}


void cuttingPlane::showPoints()
{

    Info << "Cutting points:" << endl;
    forAll(cutCells_, cellI)
    {
        Info << "Cut cell #" << cellI << " (" << cutCells_[cellI]
             << ") has " << cutCellCuttingPoints_[cellI].size()
             << " cutting points."
             << endl;

        forAll(cutCellCuttingPoints_[cellI], pointI)
        {
            Info << " #"<< pointI << ": " 
                 << cuttingPoints_[cutCellCuttingPoints_[cellI][pointI]]
                 << " (" 
                 << cutCellCuttingPoints_[cellI][pointI]
                 << ")" << endl;
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from components
cuttingPlane::cuttingPlane(const polyMesh& mesh, const plane& newPlane)
:
    plane(newPlane),
    mesh_(mesh)
{
    const faceList& Faces = mesh_.faces();
    const cellList& CellFaces = mesh_.cells();
    const pointField& Points = mesh_.points();
    const labelListList& CellEdges = mesh_.cellEdges();

    scalarField dotProducts = (Points - refPoint()) & normal();

    cutCells
    (
        dotProducts,
        CellFaces,
        Faces
    );

    labelListList intersections
    (    
        cuttingPoints
        (
            dotProducts,
            CellEdges,
            Points
        )
    );

    cutCellCuttingPoints
    (
        intersections, 
        CellEdges, 
        CellFaces,
        Faces
    );

    orderFacePoints();

    //showPoints();
    
    cutCellCuttingPoints_.setSize(0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return vectorField of cutting points
const pointField& cuttingPlane::points() const
{
    return cuttingPoints_;
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
