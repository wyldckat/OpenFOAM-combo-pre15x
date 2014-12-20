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

\*---------------------------------------------------------------------------*/

#include "parallelInfo.H"
#include "Time.H"
#include "Pstream.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"
#include "globalPoints.H"
#include "labelIOList.H"
#include "PackedList.H"
#include "mergePoints.H"
#include "matchPoints.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parallelInfo, 0);
}

// Geometric matching tolerance. Factor of mesh bounding box.
const Foam::scalar Foam::parallelInfo::matchTol_ = 1E-10;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::parallelInfo::initProcAddr()
{
    processorPatchIndices_.setSize(mesh_.boundaryMesh().size());
    processorPatchIndices_ = -1;

    processorPatchNeighbours_.setSize(mesh_.boundaryMesh().size());
    processorPatchNeighbours_ = -1;

    // Construct processor patch indexing. processorPatchNeighbours_ only
    // set if running in parallel!
    processorPatches_.setSize(mesh_.boundaryMesh().size());

    label nNeighbours = 0;

    forAll (mesh_.boundaryMesh(), patchi)
    {
        if
        (
            typeid(mesh_.boundaryMesh()[patchi])
         == typeid(processorPolyPatch)
        )
        {
            processorPatches_[nNeighbours] = patchi;
            processorPatchIndices_[patchi] = nNeighbours++;
        }
    }
    processorPatches_.setSize(nNeighbours);


    if (Pstream::parRun())
    {
        // Send indices of my processor patches to my neighbours
        forAll (processorPatches_, i)
        {
            label patchi = processorPatches_[i];

            OPstream toNeighbour
            (
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo()
            );

            toNeighbour << processorPatchIndices_[patchi];
        }

        forAll(processorPatches_, i)
        {
            label patchi = processorPatches_[i];
            
            IPstream fromNeighbour
            (
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo()
            );
            
            fromNeighbour >> processorPatchNeighbours_[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components. For testing only!
Foam::parallelInfo::parallelInfo
(
    const polyMesh& mesh,
    const bool parallel,
    const bool cyclicParallel,
    const label nTotalPoints,
    const label nTotalFaces,
    const label nTotalCells,
    const label nGlobalPoints,
    const labelList& sharedPointLabels,
    const labelList& sharedPointAddr,
    const labelList& sharedPointGlobalLabels
)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    cyclicParallel_(cyclicParallel),
    bb_(vector::zero, vector::zero),
    nTotalPoints_(nTotalPoints),
    nTotalFaces_(nTotalFaces),
    nTotalCells_(nTotalCells),
    processorPatches_(0),
    processorPatchIndices_(0),
    processorPatchNeighbours_(0),
    nGlobalPoints_(nGlobalPoints),
    sharedPointLabels_(sharedPointLabels),
    sharedPointAddr_(sharedPointAddr),
    sharedPointGlobalLabelsPtr_(new labelList(sharedPointGlobalLabels))
{
    initProcAddr();
}


// Construct from components
Foam::parallelInfo::parallelInfo(const polyMesh& mesh)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    cyclicParallel_(false),
    bb_(vector::zero, vector::zero),
    nTotalPoints_(-1),
    nTotalFaces_(-1),
    nTotalCells_(-1),
    processorPatches_(0),
    processorPatchIndices_(0),
    processorPatchNeighbours_(0),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0),
    sharedPointGlobalLabelsPtr_(NULL)
{
    // Dummy map - not used.
    mapPolyMesh dummyMap
    (
        mesh,
        0,
        0,
        0,
        labelList(0),
        labelList(0),
        List<objectMap>(0),
        List<objectMap>(0),
        labelList(0),
        List<objectMap>(0),
        List<objectMap>(0),
        List<objectMap>(0),
        labelList(0),
        labelList(0),
        labelList(0),
        labelHashSet(0),
        labelListList(0),
        labelListList(0),
        labelListList(0),
        labelListList(0),
        labelListList(0),
        pointField(0),
        labelList(1, 0),
        labelList(0)
    );

    updateTopology(dummyMap);
}


// Read constructor given IOobject and a polyMesh reference
Foam::parallelInfo::parallelInfo(const IOobject& io, const polyMesh& mesh)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    cyclicParallel_(false),
    bb_(mesh.points()),
    sharedPointGlobalLabelsPtr_(NULL)
{
    initProcAddr();

    IOdictionary dict(io);

    dict.lookup("cyclicParallel") >> cyclicParallel_;
    dict.lookup("nTotalPoints") >> nTotalPoints_;
    dict.lookup("nTotalFaces") >> nTotalFaces_;
    dict.lookup("nTotalCells") >> nTotalCells_;
    dict.lookup("nGlobalPoints") >> nGlobalPoints_;
    dict.lookup("sharedPointLabels") >> sharedPointLabels_;
    dict.lookup("sharedPointAddr") >> sharedPointAddr_;
    labelList sharedPointGlobalLabels(dict.lookup("sharedPointGlobalLabels"));

    sharedPointGlobalLabelsPtr_ = new labelList(sharedPointGlobalLabels);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parallelInfo::~parallelInfo()
{
    clearOut();
}


void Foam::parallelInfo::clearOut()
{
    deleteDemandDrivenData(sharedPointGlobalLabelsPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return shared point global labels.
const Foam::labelList& Foam::parallelInfo::sharedPointGlobalLabels() const
{
    if (!sharedPointGlobalLabelsPtr_)
    {
        sharedPointGlobalLabelsPtr_ = new labelList(sharedPointLabels_.size());
        labelList& sharedPointGlobalLabels = *sharedPointGlobalLabelsPtr_;

        IOobject addrHeader
        (
            "pointProcAddressing",
            mesh_.cellsInstance()/mesh_.meshSubDir,
            mesh_,
            IOobject::MUST_READ
        );

        if (addrHeader.headerOk())
        {
            // There is a pointProcAddressing file so use it to get labels
            // on the original mesh
            Sout<< "parallelInfo::sharedPointGlobalLabels : "
                << "Reading pointProcAddressing" << endl;

            labelIOList pointProcAddressing(addrHeader);

            forAll(sharedPointLabels_, i)
            {
                // Get my mesh point
                label pointI = sharedPointLabels_[i];

                // Map to mesh point of original mesh
                sharedPointGlobalLabels[i] = pointProcAddressing[pointI];
            }
        }
        else
        {
            Sout<< "parallelInfo::sharedPointGlobalLabels :"
                << " Setting pointProcAddressing to -1" << endl;

            sharedPointGlobalLabels = -1;
        }
    }
    return *sharedPointGlobalLabelsPtr_;
}


// Collect coordinates of shared points. (does parallel communication!)
Foam::pointField Foam::parallelInfo::sharedPoints() const
{
    //Get all processors to send their shared points to master.
    // (not very efficient)

    pointField sharedPoints(nGlobalPoints_);

    if (Pstream::master())
    {
        // Master:
        // insert my own data first
        forAll(sharedPointLabels_, i)
        {
            label sharedPointI = sharedPointAddr_[i];

            sharedPoints[sharedPointI] = mesh_.points()[sharedPointLabels_[i]];
        }

        // Receive data from slaves and insert
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave(slave);

            labelList nbrSharedPointAddr;
            pointField nbrSharedPoints;
            fromSlave >> nbrSharedPointAddr >> nbrSharedPoints;

            forAll(nbrSharedPointAddr, i)
            {
                label sharedPointI = nbrSharedPointAddr[i];

                sharedPoints[sharedPointI] = nbrSharedPoints[i];
            }
        }

        // Send back
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave(slave, sharedPoints.size()*sizeof(vector::zero));
            toSlave << sharedPoints;
        }
    }
    else
    {
        // Slave:
        // send points
        {
            OPstream toMaster(Pstream::masterNo());

            toMaster
                << sharedPointAddr_ 
                << IndirectList<point>(mesh_.points(), sharedPointLabels_)();
        }

        // Receive sharedPoints
        {
            IPstream fromMaster(Pstream::masterNo());
            fromMaster >> sharedPoints;
        }
    }

    return sharedPoints;
}


// Collect coordinates of shared points. (does parallel communication!)
Foam::pointField Foam::parallelInfo::geometricSharedPoints() const
{
    // Get coords of my shared points
    pointField sharedPoints(sharedPointLabels_.size());

    forAll(sharedPointLabels_, i)
    {
        label meshPointI = sharedPointLabels_[i];

        sharedPoints[i] = mesh_.points()[meshPointI];
    }

    // Append from all processors
    combineReduce(sharedPoints, plusEqOp<pointField>());

    // Merge tolerance
    scalar tolDim = matchTol_*mag(bb_.max() - bb_.min());

    // And see how many are unique
    labelList pMap;
    pointField mergedPoints;

    mergePoints
    (
        sharedPoints,   // coordinates to merge
        tolDim,         // tolerance
        false,          // verbosity
        pMap,
        mergedPoints
    );

    return mergedPoints;
}


void Foam::parallelInfo::movePoints(const pointField& newPoints)
{}


// Update all data after morph
void Foam::parallelInfo::updateTopology(const mapPolyMesh& map)
{
    // Clear out old data
    clearOut();

    // Do processor patch addressing
    initProcAddr();

    // CyclicParallel. Is equivalent to having 'separated' processor patches.
    bool separatedPatches = false;

    forAll(processorPatches_, i)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[processorPatches_[i]];

        separatedPatches |= refCast<const processorPolyPatch>(pp).separated();
    }

    reduce(separatedPatches, orOp<bool>());

    cyclicParallel_ = separatedPatches;

    if (debug)
    {
        Sout<< "parallelInfo : cyclicParallel_:" << cyclicParallel_ << endl;
    }

    {
        // Calculate all shared points. This does all the hard work.
        globalPoints parallelPoints(mesh_);

        // Copy data out.
        nGlobalPoints_ = parallelPoints.nGlobalPoints();
        sharedPointLabels_ = parallelPoints.sharedPointLabels();
        sharedPointAddr_ = parallelPoints.sharedPointAddr();
    }


    // Bounding box (does communication)
    bb_ = boundBox(mesh_.points());

    scalar tolDim = matchTol_*mag(bb_.max() - bb_.min());

    if (debug)
    {
        Sout<< "parallelInfo : bb_:" << bb_ << " merge dist:" << tolDim << endl;
    }


    // Total number of faces. Start off from all faces. Remove coincident
    // processor faces (on highest numbered processor) before summing.
    nTotalFaces_ = mesh_.nFaces();

    // Do not count processorpatch faces that are coincident.
    forAll(processorPatches_, i)
    {
        label patchI = processorPatches_[i];

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

        if (Pstream::myProcNo() > procPatch.neighbProcNo())
        {
            // Uncount my faces. Handle cyclics separately.

            if (procPatch.separated())
            {
                const vectorField& separationDist = procPatch.separation();

                forAll(separationDist, faceI)
                {
                    if (mag(separationDist[faceI]) < tolDim)
                    {
                        // Faces coincide -> are proper processor patch
                        nTotalFaces_--;
                    }
                    else
                    {
                        // Distance between faces -> originate from cyclic.
                        // -> keep.
                    }
                }
            }
            else
            {
                // Normal, unseparated processor patch. Remove duplicates.
                nTotalFaces_ -= procPatch.size();
            }
        }
    }
    reduce(nTotalFaces_, sumOp<label>());

    if (debug)
    {
        Sout<< "parallelInfo : nTotalFaces_:" << nTotalFaces_ << endl;
    }


    nTotalCells_ = mesh_.nCells();
    reduce(nTotalCells_, sumOp<label>());

    if (debug)
    {
        Sout<< "parallelInfo : nTotalCells_:" << nTotalCells_ << endl;
    }

    nTotalPoints_ = mesh_.nPoints();

    // Correct points for duplicate ones. We have
    // - points shared between 2 processor patches only. Count only on
    //   lower numbered processor. Make sure to count only once since points
    //   can be on multiple patches on the same processor.
    // - globally shared points.

    if (Pstream::parRun())
    {
        const label UNSET = 0;      // not set
        const label SHARED = 1;     // globally shared
        const label VISITED = 2;    // corrected for

        // Mark globally shared points
        PackedList<2> pointStatus(mesh_.nPoints(), UNSET);

        forAll(sharedPointLabels_, i)
        {
            label meshPointI = sharedPointLabels_[i];

            pointStatus.set(meshPointI, SHARED);
        }

        // Send patch local points
        forAll(processorPatches_, i)
        {
            label patchI = processorPatches_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            OPstream toNeighbour(procPatch.neighbProcNo());

            toNeighbour << procPatch.localPoints();
        }

        // Receive patch local points and uncount if coincident (and not shared)
        forAll(processorPatches_, i)
        {
            label patchI = processorPatches_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            IPstream fromNeighbour(procPatch.neighbProcNo());

            pointField nbrPoints(fromNeighbour);

            if (Pstream::myProcNo() > procPatch.neighbProcNo())
            {
                labelList pMap;
                matchPoints
                (
                    procPatch.localPoints(),
                    nbrPoints,
                    scalarField(procPatch.nPoints(), tolDim),   // tolerance
                    false,      // verbosity
                    pMap        // map from my points to nbrPoints
                );

                forAll(pMap, patchPointI)
                {
                    label meshPointI = procPatch.meshPoints()[patchPointI];

                    label stat = pointStatus.get(meshPointI);

                    if (stat == UNSET)
                    {
                        // Mark point as visited so if point is on multiple proc
                        // patches it only gets uncounted once.
                        pointStatus.set(meshPointI, VISITED);

                        if (pMap[patchPointI] != -1)
                        {
                            // Points share same coordinate so uncount.
                            nTotalPoints_--;
                        }
                    }
                }
            }
        }
        // Sum all points
        reduce(nTotalPoints_, sumOp<label>());
    }

    // nTotalPoints has not been corrected yet for shared points. For these
    // just collect all their coordinates and count unique ones.

    label mySharedPoints = sharedPointLabels_.size();
    reduce(mySharedPoints, sumOp<label>());

    // Collect and merge shared points (does parallel communication)
    pointField geomSharedPoints(geometricSharedPoints());
    label nGeomSharedPoints = geomSharedPoints.size();

    // Shared points merged down to mergedPoints size.
    nTotalPoints_ -= mySharedPoints - nGeomSharedPoints;

    if (debug)
    {
        Sout<< "parallelInfo : nTotalPoints_:" << nTotalPoints_ << endl;
    }

    //
    // Now we have all info we wanted.
    // Do some checking (if debug is set)
    //

    if (debug)
    {
        if (Pstream::master())
        {
            // We have the geometricSharedPoints already so write them.
            // Ideally would like to write the networks of connected points as
            // well but this is harder. (Todo)
            Sout<< "parallelInfo : writing geometrically separated shared"
                << " points to geomSharedPoints.obj" << endl;

            OFstream str("geomSharedPoints.obj");

            forAll(geomSharedPoints, i)
            {
                const point& pt = geomSharedPoints[i];

                str << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
                    << nl;
            }
        }

        // Compare my global data against what I can read from parallelData
        // file.
        parallelInfo oldData
        (
            IOobject
            (
                "parallelData",
                mesh_.time().findInstance
                (
                    mesh_.dbDir()/mesh_.meshSubDir,
                    "parallelData"
                ),
                mesh_.meshSubDir,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        if (oldData.cyclicParallel() != cyclicParallel_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "cyclicParallel : old:" << oldData.cyclicParallel()
                << " new:" << cyclicParallel_ << abort(FatalError);
        }
        if (oldData.nTotalPoints() != nTotalPoints_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "nTotalPoints : old:" << oldData.nTotalPoints()
                << " new:" << nTotalPoints_ << abort(FatalError);
        }
        if (oldData.nTotalFaces() != nTotalFaces_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "nTotalFaces : old:" << oldData.nTotalFaces()
                << " new:" << nTotalFaces_ << abort(FatalError);
        }
        if (oldData.nTotalCells() != nTotalCells_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "nTotalPoints : old:" << oldData.nTotalCells()
                << " new:" << nTotalCells_ << abort(FatalError);
        }
        if (oldData.nGlobalPoints() != nGlobalPoints_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "nGlobalPoints : old:" << oldData.nGlobalPoints()
                << " new:" << nGlobalPoints_ << abort(FatalError);
        }
        if (oldData.sharedPointLabels() != sharedPointLabels_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "sharedPointLabels : old:" << oldData.sharedPointLabels()
                << " new:" << sharedPointLabels_ << abort(FatalError);
        }
        if (oldData.sharedPointAddr() != sharedPointAddr_)
        {
            FatalErrorIn("parallelInfo::updateTopology(const mapPolyMesh& map)")
                << "sharedPointAddr : old:" << oldData.sharedPointAddr()
                << " new:" << sharedPointAddr_ << abort(FatalError);
        }
        // Cannot compare sharedPointGlobalLabels easily.

    }
}


// Write data
bool Foam::parallelInfo::write() const
{
    IOdictionary dict
    (
        IOobject
        (
            "parallelData",
            mesh_.cellsInstance(),
            mesh_.meshSubDir,
            mesh_
        )
    );

    dict.add("cyclicParallel", cyclicParallel_);
    dict.add("nTotalPoints", nTotalPoints());
    dict.add("nTotalFaces", nTotalFaces());
    dict.add("nTotalCells", nTotalCells());

    dict.add("nGlobalPoints", nGlobalPoints());
    dict.add("sharedPointLabels", sharedPointLabels());
    dict.add("sharedPointAddr", sharedPointAddr());
    dict.add("sharedPointGlobalLabels", sharedPointGlobalLabels());

    return dict.regIOobject::write
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const parallelInfo& p)
{
    os  << "cyclicParallel " << p.cyclicParallel_ << token::END_STATEMENT << nl
        << "nTotalPoints " << p.nTotalPoints() << token::END_STATEMENT << nl
        << "nTotalFaces " << p.nTotalFaces() << token::END_STATEMENT << nl
        << "nTotalCells " << p.nTotalCells() << token::END_STATEMENT << nl
        << "nGlobalPoints " << p.nGlobalPoints() << token::END_STATEMENT << nl
        << "sharedPointLabels " << p.sharedPointLabels()
        << token::END_STATEMENT << nl
        << "sharedPointAddr " << p.sharedPointAddr()
        << token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //
