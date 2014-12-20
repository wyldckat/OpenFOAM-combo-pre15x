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

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "parallelInfo.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMesh::clearFaceCells() const
{
    if (debug)
    {
        Info<< "void polyMesh::clearFaceCells() const: "
            << "removing mesh modifiers"
            << endl;
    }

    // Delete addressing used by primitiveMesh
    deleteDemandDrivenData(allOwnerPtr_);
    deleteDemandDrivenData(allNeighbourPtr_);
}


void polyMesh::removeBoundary()
{
    if (debug)
    {
        Info<< "void polyMesh::removeBoundary(): "
            << "Removing boundary patches."
            << endl;
    }

    // Remove the point zones
    boundary_.clear();
    boundary_.setSize(0);

    clearOut();
}


void polyMesh::removePointZones()
{
    if (debug)
    {
        Info<< "void polyMesh::removePointZones(): "
            << "Removing point zones."
            << endl;
    }

    // Remove the point zones
    pointZones_.clear();

    // Reset the size of the points list
    points_.setSize(nPoints());

    // Reset the primitive mesh
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

    clearOut();
}


void polyMesh::removeFaceZones()
{
    if (debug)
    {
        Info<< "void polyMesh::removeFaceZones(): "
            << "Removing face zones."
            << endl;
    }

    // Remove the face zones
    faceZones_.clear();

    // Reset the size of the faces list
    faces_.setSize(nFaces());

    // Re-do the boundary patches with the new face list
    List<polyPatch*> newPatches(boundary_.size());

    forAll (boundary_, patchI)
    {
        newPatches[patchI] =
            boundary_[patchI].clone
            (
                boundary_
            ).ptr();
    }

    boundary_.clear();
    addPatches(newPatches);

    // Reset the primitive mesh
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

    clearOut();
}


void polyMesh::removeCellZones()
{
    if (debug)
    {
        Info<< "void polyMesh::removeCellZones(): "
            << "Removing cell zones."
            << endl;
    }

    // Remove the cell zones
    cellZones_.clear();

    // Reset the size of the cells list
    cells_.setSize(nCells());

    // Reset the primitive mesh
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

    clearOut();
}


void polyMesh::removeMeshModifiers()
{
    if (debug)
    {
        Info<< "void polyMesh::removeMorphEngine(): "
            << "removing mesh modifiers"
            << endl;
    }

    deleteDemandDrivenData(morphEnginePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMesh::printAllocated() const
{
    Info<< "polyMesh allocated :" << endl;

    if (allOwnerPtr_)
    {
        Info<< "    All owner" << endl;
    }

    if (allNeighbourPtr_)
    {
        Info<< "    All neighbour" << endl;
    }

    if (oldPointsPtr_)
    {
        Info<< "    Old motion points" << endl;
    }

    if (morphMap_)
    {
        Info<< "    Topology morph map" << endl;
    }
}


void polyMesh::clearGeom()
{
    if (debug)
    {
        Info<< "void polyMesh::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    primitiveMesh::clearGeom();

    forAll (boundary_, patchI)
    {
        boundary_[patchI].clearGeom();
    }
}


void polyMesh::clearAddressing()
{
    if (debug)
    {
        Info<< "void polyMesh::clearAddressing() : "
            << "clearing topology"
            << endl;
    }

    primitiveMesh::clearAddressing();

    // parallelData depends on the processorPatch ordering so force
    // recalculation
    deleteDemandDrivenData(parallelDataPtr_);
}


// Clear primitive data (points and cells)
void polyMesh::clearPrimitives()
{
    resetMotion();
    resetMorph();

    points_.setSize(0);
    faces_.setSize(0);
    cells_.setSize(0);
}


void polyMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
