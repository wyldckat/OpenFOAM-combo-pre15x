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
#include "Time.H"
#include "primitiveMesh.H"
#include "cellIOList.H"
#include "SubList.H"
#include "polyPatch.H"
#include "parallelInfo.H"
#include "processorPolyPatch.H"
#include "OSspecific.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyMesh, 0);

word polyMesh::defaultRegion = "region0";
word polyMesh::meshSubDir = "polyMesh";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOobject
// Search back in time for points, cells and boundary file
// if not found read from constant directory
// WARNING: ASSUMES CORRECT ORDERING OF DATA. 
polyMesh::polyMesh(const IOobject& io)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    faces_
    (
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    allOwner_
    (
        IOobject
        (
            "owner",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    allNeighbour_
    (
        IOobject
        (
            "neighbour",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "boundary",
            time().findInstance(meshDir(), "boundary"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            time().findInstance
            (
                meshDir(),
                "pointZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            time().findInstance
            (
                meshDir(),
                "faceZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            time().findInstance
            (
                meshDir(),
                "cellZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    parallelDataPtr_(NULL),
    moving_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    if (allOwner_.size())
    {
        initMesh();
    }
    else
    {
        cellIOList c
        (
            IOobject
            (
                "cells",
                time().findInstance(meshDir(), "cells"),
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Set the primitive mesh
        initMesh(c);

        allOwner_.write();
        allNeighbour_.write();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// Construct from components without boundary.
// Boundary is added using addPatches() member function
// WARNING: ASSUMES CORRECT ORDERING OF DATA. 
polyMesh::polyMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    objectRegistry(io),
    primitiveMesh(),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faces
    ),
    allOwner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    allNeighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        0
    ),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    parallelDataPtr_(NULL),
    moving_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldPointsPtr_(NULL)
{
    // Check if the faces and cells are valid
    forAll (faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject& io,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const cellList& cells\n"
                ")\n"
            )   << "Face " << faceI << "contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }

    // Check if the faces and cells are valid
    forAll (cells, cellI)
    {
        const cell& curCell = cells[cellI];

        if (min(curCell) < 0 || max(curCell) > faces_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject& io,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const cellList& cells\n"
                ")\n"
            )   << "Cell " << cellI << "contains face labels out of range: "
                << curCell << " Max face index = " << faces_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh(const_cast<cellList&>(cells));
}



// Reset mesh primitive data
// WARNING: ASSUMES CORRECT ORDERING OF DATA.
void polyMesh::resetPrimitives
(
    const label nUsedFaces,
    const pointField& points,
    const faceList& faces,
    const labelList& owner,
    const labelList& neighbour,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const bool validBoundary
)
{
    // Clear everything (copied from ~polyMesh)
    clearOut();
    resetMotion();

    // Take over new primitive data. Note extra optimization to prevent
    // assignment to self.
    if (&points_ != &points)
    {
        points_ = points;
    }
    if (&faces_ != &faces)
    {
        faces_ = faces;
    }
    if (&allOwner_ != &owner)
    {
        allOwner_ = owner;
    }
    if (&allNeighbour_ != &neighbour)
    {
        allNeighbour_ = neighbour;
    }

    // Reset patch sizes and starts
    forAll(boundary_, patchI)
    {
        boundary_[patchI] = polyPatch
        (
            boundary_[patchI].name(),
            patchSizes[patchI],
            patchStarts[patchI],
            patchI,
            boundary_
        );
    }


    // Flags the mesh files as being changed
    setInstance(time().timeName());

    // Check if the faces and cells are valid
    forAll (faces_, faceI)
    {
        const face& curFace = faces_[faceI];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh::resetPrimitives\n"
                "(\n"
                "    const label nUsedFaces,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const labelList& owner,\n"
                "    const labelList& neighbour,\n"
                "    const labelList& patchSizes,\n"
                "    const labelList& patchStarts\n"
                ")\n"
            )   << "Face " << faceI << "contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }


    // Set the primitive mesh from the allOwner_, allNeighbour_. Works
    // out from patch end where the active faces stop.
    initMesh();


    if (validBoundary)
    {
        // Note that we assume that all the patches stay the same and are
        // correct etc. so we can already use the patches to do
        // processor-processor comms.

        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMesh::~polyMesh()
{
    clearOut();
    resetMotion();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const fileName& polyMesh::dbDir() const
{
    if (objectRegistry::dbDir() == defaultRegion)
    {
        return parent().dbDir();
    }
    else
    {
        return objectRegistry::dbDir();
    }
}


fileName polyMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const fileName& polyMesh::pointsInstance() const
{
    return points_.instance();
}


const fileName& polyMesh::facesInstance() const
{
    return faces_.instance();
}


// Add boundary patches. Constructor helper
void polyMesh::addPatches(const List<polyPatch*> & p)
{
    if (boundaryMesh().size() > 0)
    {
        FatalErrorIn("void polyMesh::addPatches(const List<polyPatch*>& p)")
            << "boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    // Copy the patch pointers
    forAll (p, pI)
    {
        boundary_.hook(p[pI]);
    }

    // parallelData depends on the processorPatch ordering so force
    // recalculation. Problem: should really be done in removeBoundary but
    // there is some info in parallelData which might be interesting inbetween
    // removeBoundary and addPatches.
    deleteDemandDrivenData(parallelDataPtr_);

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    boundary_.checkDefinition();
}


// Add mesh zones. Constructor helper
void polyMesh::addZones
(
    const List<pointZone*>& pz,
    const List<faceZone*>& fz,
    const List<cellZone*>& cz
)
{
    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        FatalErrorIn
        (
            "void addZones\n"
            "(\n"
            "    const List<pointZone*>& pz,\n"
            "    const List<faceZone*>& fz,\n"
            "    const List<cellZone*>& cz\n"
            ")"
        )   << "point, face or cell zone already exists"
            << abort(FatalError);
    }

    // Point zones
    if (pz.size())
    {
        pointZones_.setSize(pz.size());

        // Copy the zone pointers
        forAll (pz, pI)
        {
            pointZones_.hook(pz[pI]);
        }

        pointZones_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // Face zones
    if (fz.size())
    {
        faceZones_.setSize(fz.size());

        // Copy the zone pointers
        forAll (fz, fI)
        {
            faceZones_.hook(fz[fI]);
        }

        faceZones_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // Cell zones
    if (cz.size())
    {
        cellZones_.setSize(cz.size());

        // Copy the zone pointers
        forAll (cz, cI)
        {
            cellZones_.hook(cz[cI]);
        }

        cellZones_.writeOpt() = IOobject::AUTO_WRITE;
    }
}


const pointField& polyMesh::allPoints() const
{
    if (points_.size() == 0)
    {
        FatalErrorIn("const pointField& polyMesh::allPoints() const")
            << "points deallocated"
            << abort(FatalError);
    }

    return points_;
}


const faceList& polyMesh::allFaces() const
{
    if (faces_.size() == 0)
    {
        FatalErrorIn("const faceList& polyMesh::allFaces() const")
            << "faces deallocated"
            << abort(FatalError);
    }

    return faces_;
}


const labelList& polyMesh::allOwner() const
{
    return allOwner_;
}


const labelList& polyMesh::allNeighbour() const
{
    return allNeighbour_;
}


// Return old mesh motion points
const pointField& polyMesh::oldAllPoints() const
{
    if (!oldPointsPtr_)
    {
        if (debug)
        {
            WarningIn("const pointField& polyMesh::oldAllPoints() const")
                << "Old points not available.  Forcing storage of old points"
                << endl;
        }

        oldPointsPtr_ = new pointField(points_);
        curMotionTimeIndex_ = time().timeIndex();
    }

    return *oldPointsPtr_;
}


// Move points
tmp<scalarField> polyMesh::movePoints(const pointField& newPoints)
{
    if (debug)
    {
        Info<< "tmp<scalarField> polyMesh::movePoints(const pointField&) : "
            << " Moving points for time " << time().value()
            << " index " << time().timeIndex() << endl;
    }

    moving_ = true;

    // Pick up old points
    if (curMotionTimeIndex_ != time().timeIndex())
    {
        // Mesh motion in the new time step
        deleteDemandDrivenData(oldPointsPtr_);
        oldPointsPtr_ = new pointField(points_);
        curMotionTimeIndex_ = time().timeIndex();
    }

    points_ = newPoints;

    // Adjust parallel shared points
    if (parallelDataPtr_)
    {
        parallelDataPtr_->movePoints(points_);
    }

    if (debug)
    {
        // Check mesh motion
        if (primitiveMesh::checkMeshMotion(points_, true))
        {
            Info<< "tmp<scalarField> polyMesh::movePoints"
                << "(const pointField&) : "
                << "Moving the mesh with given points will "
                << "invalidate the mesh." << nl
                << "Mesh motion should not be executed." << endl;
        }
    }

    points_.writeOpt() = IOobject::AUTO_WRITE;
    points_.instance() = time().timeName();

    // Force recalculation of all geometric data with new points
    boundary_.movePoints(allPoints());

    pointZones_.movePoints(allPoints());
    faceZones_.movePoints(allPoints());
    cellZones_.movePoints(allPoints());

    primitiveMesh::clearGeom();
    return primitiveMesh::movePoints(allPoints(), oldAllPoints());
}


// Reset motion by deleting old points
void polyMesh::resetMotion() const
{
    curMotionTimeIndex_ = 0;
    deleteDemandDrivenData(oldPointsPtr_);
}


// Return parallel info
const parallelInfo& polyMesh::parallelData() const
{
    if (!parallelDataPtr_)
    {
        if (debug)
        {
            Info<< "polyMesh::parallelData() const : "
                << "Constructing parallelData from processor topology" << nl
                << "This needs the patch faces to be correctly matched"
                << endl;
        }
        // Construct parallelInfo using processorPatch information only.
        parallelDataPtr_ = new parallelInfo(*this);
    }

    return *parallelDataPtr_;
}


// Remove all files
void polyMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = db().path()/instanceDir/meshSubDir;

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"owner");
    rm(meshFilesPath/"neighbour");
    rm(meshFilesPath/"cells");
    rm(meshFilesPath/"boundary");
    rm(meshFilesPath/"pointZones");
    rm(meshFilesPath/"faceZones");
    rm(meshFilesPath/"cellZones");
    rm(meshFilesPath/"meshModifiers");
    rm(meshFilesPath/"parallelData");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
