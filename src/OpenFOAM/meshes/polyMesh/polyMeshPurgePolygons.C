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
    Purge polygonal face descriptions from extra points.
    Additionally, all point in point zones are kept.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "pointZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMesh::purgePolygons()
{
    // Algorithm:
    // Create pointEdges and mark all points with the number of edges
    // using them.  Remove all points used by less than 3 edges; two should
    // be minimum.  pointZone points are kept because of the increment.

    Info<< "void polyMesh::purgePolygons() : "
        << "purging polygons" << endl;

    labelList pointsUsage(points().size(), 0);

    const labelListList& pe = pointEdges();

    forAll (pointsUsage, pointI)
    {
        pointsUsage[pointI] = pe[pointI].size();
    }

    forAll (pointZones_, pzI)
    {
        const labelList& pointsToKeep = pointZones_[pzI].addressing();

        forAll (pointsToKeep, keepI)
        {
            pointsUsage[pointsToKeep[keepI]]++;
        }
    }

    bool meshChanged = false;

    // Check if the mesh has changed before renumbering
    forAll (pointsUsage, pointI)
    {
        if (pointsUsage[pointI] < 3)
        {
            meshChanged = true;

            break;
        }
    }

    if (!meshChanged)
    {
        Info<< "void polyMesh::purgePolygons() : "
            << "mesh unchanged" << endl;

        return;
    }
    else
    {
        label nNewPoints = 0;

        // Sort out the primitive points using point maps
        labelList pointMap(pointsUsage.size(), -1);

        // Memory management
        {
            // Create the new points list
            const pointField oldPoints = allPoints();

            forAll (pointsUsage, pointI)
            {
                if (pointsUsage[pointI] >= 3)
                {
                    // The point is kept. Add it to the list
                    points_[nNewPoints] = oldPoints[pointI];

                    pointMap[pointI] = nNewPoints;
                    nNewPoints++;
                }
            }

            // Add trailing points to the list with no change
            for (label p = pointsUsage.size(); p < oldPoints.size(); p++)
            {
                // The point is kept. Add it to the list
                points_[nNewPoints] = oldPoints[p];

                pointMap[p] = nNewPoints;
                nNewPoints++;
            }

            points_.setSize(nNewPoints);
        }

        // Create a new face list with renumbered points

        forAll (faces_, faceI)
        {
            const labelList oldLabels = faces_[faceI];

            labelList& newFace = faces_[faceI];
            label nNewLabels = 0;

            forAll (oldLabels, i)
            {
                if (pointMap[oldLabels[i]] > -1)
                {
                    // Point is kept
                    newFace[nNewLabels] = pointMap[oldLabels[i]];
                    nNewLabels++;
                }
            }

            // If less than three point in face, something has gone wrong
            if (nNewLabels < 3)
            {
                FatalErrorIn
                (
                    "void polyMesh::purgePolygons(const labelList&)"
                )   << "number of points in face " << faceI
                    << " less than three: " << nNewLabels
                    << "purge failed!"
                    << abort(FatalError);
            }

            newFace.setSize(nNewLabels);
        }

        List<pointZone*> newPointZones(pointZones_.size());

        // Renumber the point zones
        forAll (pointZones_, pzI)
        {
            const labelList& oldPointLabels = pointZones_[pzI].addressing();

            labelList newPointLabels(oldPointLabels.size());

            forAll (oldPointLabels, i)
            {
                newPointLabels = pointMap[oldPointLabels[i]];
            }

            newPointZones[pzI] =
                pointZones_[pzI].clone
                (
                    pointZones_, 
                    pzI,
                    newPointLabels
                ).ptr();
        }

        // Re-hook the point zones
        pointZones_.clear();
        
        forAll (newPointZones, pzI)
        {
            pointZones_.hook(newPointZones[pzI]);
        }

        // Re-create the basic mesh
        primitiveMesh::reset
        (
            nPoints(),
            nInternalFaces(),
            nFaces(),
            nCells(),
            allPoints(),
            allFaces(),
            allOwner(),
            allNeighbour()
        );

        // Flags the mesh files as being changed so they get written.
        faces_.writeOpt() = IOobject::AUTO_WRITE;
        pointZones_.writeOpt() = IOobject::AUTO_WRITE;

        clearOut();

        // Done. Check new mesh
        checkMesh();

        Info<< "void polyMesh::purgePolygons() : "
            << "finished purging polygons" << endl;

        return;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
