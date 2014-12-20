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
    write function : write points, cells and boundary mesh

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set instance for mesh files
void polyMesh::setInstance(const fileName& inst)
{
    if (debug)
    {
        Info<< "void polyMesh::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    points_.instance() = inst;
    faces_.instance() = inst;
    cells_.instance() = inst;
    boundary_.instance() = inst;

    pointZones_.instance() = inst;
    faceZones_.instance() = inst;
    cellZones_.instance() = inst;

    if (morphEnginePtr_)
    {
        morphEnginePtr_->instance() = inst;
    }
}


polyMesh::readUpdateState polyMesh::readUpdate()
{
    if (debug)
    {
        Info<< "polyMesh::readUpdateState polyMesh::readUpdate() : "
            << "Updating mesh based on saved data." << endl;
    }

    // Find the point and cell instance
    fileName pointsInst(time().findInstance(meshDir(), "points"));
    fileName cellsInst(time().findInstance(meshDir(), "cells"));

    if (debug)
    {
        Info<< "Cells instance: old = " << cellsInstance()
            << " new = " << cellsInst << nl
            << "Points instance: old = " << pointsInstance()
            << " new = " << pointsInst << endl;
    }

    if (cellsInst != cellsInstance())
    {
        // Topological change
        if (debug)
        {
            Info << "Topological change" << endl;
        }

        clearOut();

        // Set instance to new instance
        setInstance(cellsInst);

        points_ =
            pointIOField
            (
                IOobject
                (
                    "points",
                    cellsInst,
                    meshSubDir,
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        faces_ =
            faceIOList
            (
                IOobject
                (
                    "faces",
                    cellsInst,
                    meshSubDir,
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        cells_ =
            cellIOList
            (
                IOobject
                (
                    "cells",
                    cellsInst,
                    meshSubDir,
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

        // Recalculate the owner/neighbour addressing and reset the
        // primitiveMesh
        clearFaceCells();
        calcFaceCells();

        // Reset the boundary patches
        polyBoundaryMesh newBoundary
        (
            IOobject
            (
                "boundary",
                cellsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        // Check that patch types and names are unchanged
        bool boundaryChanged = false;

        if (newBoundary.size() != boundary_.size())
        {
            boundaryChanged = true;
        }
        else
        {
            wordList newTypes = newBoundary.types();
            wordList newNames = newBoundary.names();

            wordList oldTypes = boundary_.types();
            wordList oldNames = boundary_.names();

            forAll (oldTypes, patchI)
            {
                if
                (
                    oldTypes[patchI] != newTypes[patchI]
                 || oldNames[patchI] != newNames[patchI]
                )
                {
                    boundaryChanged = true;
                    break;
                }
            }
        }

        if (boundaryChanged)
        {
            Warning
                << "polyMesh::readUpdateState polyMesh::readUpdate() : "
                << "Number of patches has changed.  This may have "
                << "unexpected consequences.  Proceed with care." << endl;

            boundary_.clear();
            boundary_.setSize(newBoundary.size());

            forAll (newBoundary, patchI)
            {
                boundary_.hook(newBoundary[patchI].clone(boundary_));
            }
        }
        else
        {
            forAll (boundary_, patchI)
            {
                boundary_[patchI] = polyPatch
                (
                    newBoundary[patchI].name(),
                    newBoundary[patchI].size(),
                    newBoundary[patchI].start(),
                    patchI,
                    boundary_
                );
            }
        }

        pointZoneMesh newPointZones
        (
            IOobject
            (
                "pointZones",
                cellsInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        forAll (pointZones_, pzI)
        {
            pointZones_[pzI].resetAddressing(newPointZones[pzI].addressing());
        }

        faceZoneMesh newFaceZones
        (
            IOobject
            (
                "faceZones",
                cellsInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        forAll (faceZones_, fzI)
        {
            faceZones_[fzI].resetAddressing
            (
                newFaceZones[fzI].addressing(),
                newFaceZones[fzI].flipMap()
            );
        }

        cellZoneMesh newCellZones
        (
            IOobject
            (
                "cellZones",
                cellsInst,
                meshSubDir,
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        forAll (cellZones_, czI)
        {
            cellZones_[czI].resetAddressing(newCellZones[czI].addressing());
        }

        if (boundaryChanged)
        {
            return polyMesh::TOPO_PATCH_CHANGE;
        }
        else
        {
            return polyMesh::TOPO_CHANGE;
        }
    }
    else if (pointsInst != pointsInstance())
    {
        // Points moved
        if (debug)
        {
            Info << "Point motion" << endl;
        }

        clearGeom();

        points_.instance() = pointsInst;

        points_ =
            pointIOField
            (
                IOobject
                (
                    "points",
                    pointsInst,
                    meshSubDir,
                    *this,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

        return polyMesh::POINTS_MOVED;
    }
    else
    {
        if (debug)
        {
            Info << "No change" << endl;
        }

        return polyMesh::UNCHANGED;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


