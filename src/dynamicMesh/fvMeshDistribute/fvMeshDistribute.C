/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    fvMeshDistribute

\*----------------------------------------------------------------------------*/

#include "fvMeshDistribute.H"
#include "ProcessorTopology.H"
#include "commSchedule.H"
#include "PstreamCombineReduceOps.H"
#include "fvMeshAdder.H"
#include "faceCoupleInfo.H"
#include "processorFvPatchField.H"
#include "directPolyTopoChange.H"
#include "directRemoveCells.H"
#include "polyModifyFace.H"
#include "mapDistributePolyMesh.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fvMeshDistribute, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::List<Foam::labelPair> Foam::fvMeshDistribute::getSchedule
//(
//    const labelList& distribution
//)
//{
//    labelList nCellsPerProc(countCells(distribution));
//
//    if (debug)
//    {
//        Pout<< "getSchedule : Wanted distribution:" << nCellsPerProc << endl;
//    }
//
//    // Processors I need to send data to
//    labelListList mySendProcs(Pstream::nProcs());
//
//    // Count
//    label nSendProcs = 0;
//    forAll(nCellsPerProc, sendProcI)
//    {
//        if (sendProcI != Pstream::myProcNo() && nCellsPerProc[sendProcI] > 0)
//        {
//            nSendProcs++;
//        }
//    }
//
//    // Fill 
//    mySendProcs[Pstream::myProcNo()].setSize(nSendProcs);
//    nSendProcs = 0;
//    forAll(nCellsPerProc, sendProcI)
//    {
//        if (sendProcI != Pstream::myProcNo() && nCellsPerProc[sendProcI] > 0)
//        {
//            mySendProcs[Pstream::myProcNo()][nSendProcs++] = sendProcI;
//        }
//    }
//
//    // Synchronise
//    Pstream::gatherList(mySendProcs);
//    Pstream::scatterList(mySendProcs);
//
//    // Combine into list (same on all procs) giving sending and receiving
//    // processor
//    label nComms = 0;
//    forAll(mySendProcs, procI)
//    {
//        nComms += mySendProcs[procI].size();
//    }
//
//    List<labelPair> schedule(nComms);
//    nComms = 0;
//
//    forAll(mySendProcs, procI)
//    {
//        const labelList& sendProcs = mySendProcs[procI];
//
//        forAll(sendProcs, i)
//        {
//            schedule[nComms++] = labelPair(procI, sendProcs[i]);
//        }
//    }
//
//    return schedule;
//}


Foam::labelList Foam::fvMeshDistribute::select
(
    const bool selectEqual,
    const labelList& values,
    const label value
)
{
    label n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            n++;
        }
    }

    labelList indices(n);
    n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            indices[n++] = i;
        }
    }
    return indices;
}


// Check all procs have same names and in exactly same order.
void Foam::fvMeshDistribute::checkEqualWordList(const wordList& lst)
{
    wordList myObjects(lst);

    // Check that all meshes have same objects
    Pstream::listCombineGather(myObjects, checkEqualType());
    // Below scatter only needed to balance sends and receives.
    Pstream::listCombineScatter(myObjects);
}


// Print some info on mesh.
void Foam::fvMeshDistribute::printMeshInfo(const fvMesh& mesh)
{
    Pout<< "Primitives:" << nl
        << "    points       :" << mesh.nPoints() << nl
        << "    internalFaces:" << mesh.nInternalFaces() << nl
        << "    faces        :" << mesh.nFaces() << nl
        << "    cells        :" << mesh.nCells() << nl;

    const fvBoundaryMesh& patches = mesh.boundary();

    Pout<< "Patches:" << endl;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI].patch();

        Pout<< "    " << patchI << " name:" << pp.name()
            << " size:" << pp.size()
            << " start:" << pp.start()
            << " type:" << pp.type()
            << endl;
    }
}


void Foam::fvMeshDistribute::printCoupleInfo
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourceNewProc
)
{
    Pout<< nl
        << "Current coupling info:"
        << endl;

    forAll(sourceFace, bFaceI)
    {
        label meshFaceI = mesh.nInternalFaces() + bFaceI;

        Pout<< "    meshFace:" << meshFaceI
            << " fc:" << mesh.faceCentres()[meshFaceI]
            << " connects to proc:" << sourceProc[bFaceI]
            << "/face:" << sourceFace[bFaceI]
            << " which will move to proc:" << sourceNewProc[bFaceI]
            << endl;
    }
}


// Finds (non-empty) patch that exposed internal and proc faces can be put into.
Foam::label Foam::fvMeshDistribute::findNonEmptyPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nonEmptyPatchI = -1;

    forAllReverse(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!isA<emptyPolyPatch>(pp) && !isA<processorPolyPatch>(pp))
        {
            nonEmptyPatchI = patchI;
            break;
        }
    }

    if (nonEmptyPatchI == -1)
    {
        FatalErrorIn("findNonEmptyPatch() const")
            << "Cannot find a patch which is neither of type empty nor"
            << " of type processor in patches " << patches.names() << endl
            << "There has to be at least one such patch for"
            << " distribution to work" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "findNonEmptyPatch : using patch " << nonEmptyPatchI
            << " name:" << patches[nonEmptyPatchI].name()
            << " type:" << patches[nonEmptyPatchI].type()
            << " to put exposed faces into." << endl;
    }


    // Do additional test for processor patches intermingled with non-proc
    // patches.
    label procPatchI = -1;

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            procPatchI = patchI;
        }
        else if (procPatchI != -1)
        {
            FatalErrorIn("findNonEmptyPatch() const")
                << "Processor patches should be at end of patch list."
                << endl
                << "Have processor patch " << procPatchI
                << " followed by non-processor patch " << patchI
                << " in patches " << patches.names()
                << abort(FatalError);
        }
    }

    return nonEmptyPatchI;
}


// Appends processorPolyPatch. Returns patchID.
Foam::label Foam::fvMeshDistribute::addProcPatch
(
    const word& patchName,
    const label nbrProc
)
{
    // Clear local fields and e.g. polyMesh globalMeshData.
    mesh_.clearOut();


    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh_.boundary());

    if (polyPatches.findPatchID(patchName) != -1)
    {
        FatalErrorIn("addProcPatch(const word&, const label)")
            << "Cannot create patch " << patchName << " since already exists."
            << nl
            << "Current patch names:" << polyPatches.names()
            << exit(FatalError);
    }



    // Add the patch
    // ~~~~~~~~~~~~~

    label sz = polyPatches.size();

    // Add polyPatch
    polyPatches.setSize(sz+1);
    polyPatches.set
    (
        sz,
        new processorPolyPatch
        (
            patchName,
            0,              // size
            mesh_.nFaces(),
            sz,
            mesh_.boundaryMesh(),
            Pstream::myProcNo(),
            nbrProc
        )
    );
    fvPatches.setSize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh_.boundary()
        )
    );

    return sz;
}


// Deletes last patch
void Foam::fvMeshDistribute::deleteTrailingPatch()
{
    // Clear local fields and e.g. polyMesh globalMeshData.
    mesh_.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh_.boundary());

    if (polyPatches.size() == 0)
    {
        FatalErrorIn("deleteTrailingPatch(fvMesh&)")
            << "No patches in mesh"
            << abort(FatalError);
    }

    label sz = polyPatches.size();

    label nFaces = polyPatches[sz-1].size();

    // Remove last polyPatch
    if (debug)
    {
        Pout<< "deleteTrailingPatch : Removing patch " << sz-1
            << " : " << polyPatches[sz-1].name() << " size:" << nFaces << endl;
    }

    if (nFaces != 0)
    {
        FatalErrorIn("deleteTrailingPatch()")
            << "There are still " << nFaces << " faces in patch to be deleted "
            << sz-1 << ' ' << polyPatches[sz-1].name()
            << abort(FatalError);
    }

    // Remove actual patch
    polyPatches.setSize(sz-1);
    fvPatches.setSize(sz-1);
}


// Delete all processor patches. Move any processor faces into the last
// non-processor patch.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::deleteProcPatches
(
    const label destinationPatch
)
{
    // New patchID per boundary faces to be repatched. Is -1 (no change)
    // or new patchID
    labelList newPatchID(mesh_.nFaces() - mesh_.nInternalFaces(), -1);

    label nProcPatches = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            if (debug)
            {
                Pout<< "Moving all faces of patch " << pp.name()
                    << " into patch " << destinationPatch
                    << endl;
            }

            label offset = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                newPatchID[offset+i] = destinationPatch;
            }

            nProcPatches++;
        }
    }

    // Note: order of boundary faces been kept the same since the
    // destinationPatch is at the end and we have visited the patches in
    // incremental order.
    labelListList dummyFaceMaps;
    autoPtr<mapPolyMesh> map = repatch(newPatchID, dummyFaceMaps);


    // Delete (now empty) processor patches.
    forAllReverse(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            deleteTrailingPatch();
            deleteTrailingPatchFields<scalar, volMesh>();
            deleteTrailingPatchFields<vector, volMesh>();
            deleteTrailingPatchFields<tensor, volMesh>();
            deleteTrailingPatchFields<scalar, surfaceMesh>();
            deleteTrailingPatchFields<vector, surfaceMesh>();
            deleteTrailingPatchFields<tensor, surfaceMesh>();
        }
    }

    return map;
}


// Repatch the mesh.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::repatch
(
    const labelList& newPatchID,         // per boundary face -1 or new patchID
    labelListList& constructFaceMap
)
{
    directPolyTopoChange meshMod(mesh_);
    
    forAll(newPatchID, bFaceI)
    {
        if (newPatchID[bFaceI] != -1)
        {
            label faceI = mesh_.nInternalFaces() + bFaceI;

            label zoneID = mesh_.faceZones().whichZone(faceI);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[faceI],       // modified face
                    faceI,                      // label of face
                    mesh_.faceOwner()[faceI],   // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    newPatchID[bFaceI],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );
        }
    }


    // Do mapping of fields from one patchField to the other ourselves since
    // is currently not supported by updateMesh.

    // Store boundary fields (we only do this for surfaceFields)
    PtrList<FieldField<fvPatchField, scalar> > sFlds;
    saveBoundaryFields<scalar, surfaceMesh>(sFlds);
    PtrList<FieldField<fvPatchField, vector> > vFlds;
    saveBoundaryFields<vector, surfaceMesh>(vFlds);
    PtrList<FieldField<fvPatchField, tensor> > tFlds;
    saveBoundaryFields<tensor, surfaceMesh>(tFlds);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Map patch fields using stored boundary fields. Note: assumes order
    // of fields has not changed in object registry!
    mapBoundaryFields<scalar, surfaceMesh>(map, sFlds);
    mapBoundaryFields<vector, surfaceMesh>(map, vFlds);
    mapBoundaryFields<tensor, surfaceMesh>(map, tFlds);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        mesh_.clearOut();
    }

    // Adapt constructMaps.

    if (debug)
    {
        label index = findIndex(map().reverseFaceMap(), -1);

        if (index != -1)
        {
            FatalErrorIn
            (
                "fvMeshDistribute::repatch(const labelList&, labelListList&)"
            )   << "reverseFaceMap contains -1 at index:"
                << index << endl
                << "This means that the repatch operation was not just a shuffle?"
                << abort(FatalError);
        }
    }

    forAll(constructFaceMap, procI)
    {
        inplaceRenumber(map().reverseFaceMap(), constructFaceMap[procI]);
    }


    return map;
}


// Construct the local environment of all boundary faces.
void Foam::fvMeshDistribute::getNeighbourData
(
    const labelList& distribution,
    labelList& sourceFace,
    labelList& sourceProc,
    labelList& sourceNewProc
) const
{
    sourceFace.setSize(mesh_.nFaces() - mesh_.nInternalFaces());
    sourceProc.setSize(sourceFace.size());
    sourceNewProc.setSize(sourceFace.size());

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Send meshFace labels of processor patches and destination processor.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Labels of faces on this side
            labelList meshFaceLabels(pp.size());
            forAll(meshFaceLabels, i)
            {
                meshFaceLabels[i] = pp.start()+i;
            }

            // Which processor they will end up on
            const labelList newProc
            (
                IndirectList<label>(distribution, pp.faceCells())
            );

            OPstream toNeighbour(procPatch.neighbProcNo());

            toNeighbour << meshFaceLabels << newProc;
        }
    }

    // Receive meshFace labels and destination processors of processor faces.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        label offset = pp.start() - mesh_.nInternalFaces();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Receive the data
            IPstream fromNeighbour(procPatch.neighbProcNo());
            labelList nbrFaces(fromNeighbour);
            labelList nbrNewProc(fromNeighbour);

            // Check which of the two faces we store.

            if (Pstream::myProcNo() < procPatch.neighbProcNo())
            {
                // Use my local face labels
                forAll(pp, i)
                {
                    sourceFace[offset + i] = pp.start()+i;
                    sourceProc[offset + i] = Pstream::myProcNo();
                    sourceNewProc[offset + i] = nbrNewProc[i];
                }
            }
            else
            {
                // Use my neighbours face labels
                forAll(pp, i)
                {
                    sourceFace[offset + i] = nbrFaces[i];
                    sourceProc[offset + i] = procPatch.neighbProcNo();
                    sourceNewProc[offset + i] = nbrNewProc[i];
                }
            }
        }
        else
        {
            // Normal (physical) boundary
            forAll(pp, i)
            {
                sourceFace[offset + i] = patchI;
                sourceProc[offset + i] = -1;
                sourceNewProc[offset + i] = -1;
            }
        }
    }
}


// Subset the neighbourCell/neighbourProc fields
void Foam::fvMeshDistribute::subsetBoundaryData
(
    const fvMesh& mesh,
    const labelList& faceMap,
    const labelList& cellMap,

    const labelList& oldDistribution,
    const labelList& oldFaceOwner,
    const labelList& oldFaceNeighbour,
    const label oldInternalFaces,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourceNewProc,

    labelList& subFace,
    labelList& subProc,
    labelList& subNewProc
)
{
    subFace.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subProc.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subNewProc.setSize(mesh.nFaces() - mesh.nInternalFaces());

    forAll(subFace, newBFaceI)
    {
        label newFaceI = newBFaceI + mesh.nInternalFaces();

        label oldFaceI = faceMap[newFaceI];

        // Was oldFaceI internal face? If so which side did we get.
        if (oldFaceI < oldInternalFaces)
        {
            subFace[newBFaceI] = oldFaceI;
            subProc[newBFaceI] = Pstream::myProcNo();

            label oldOwn = oldFaceOwner[oldFaceI];
            label oldNei = oldFaceNeighbour[oldFaceI];

            if (oldOwn == cellMap[mesh.faceOwner()[newFaceI]])
            {
                // We kept the owner side. Where does the neighbour move to?
                subNewProc[newBFaceI] = oldDistribution[oldNei];
            }
            else
            {
                // We kept the neighbour side.
                subNewProc[newBFaceI] = oldDistribution[oldOwn];
            }
        }
        else
        {
            // Was boundary face. Take over boundary information
            label oldBFaceI = oldFaceI - oldInternalFaces;

            subFace[newBFaceI] = sourceFace[oldBFaceI];
            subProc[newBFaceI] = sourceProc[oldBFaceI];
            subNewProc[newBFaceI] = sourceNewProc[oldBFaceI];
        }
    }
}


// Find cells on mesh whose faceID/procID match the neighbour cell/proc of
// domainMesh. Store the matching face.
void Foam::fvMeshDistribute::findCouples
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,

    const label domain,
    const primitiveMesh& domainMesh,
    const labelList& domainFace,
    const labelList& domainProc,

    labelList& masterCoupledFaces,
    labelList& slaveCoupledFaces
)
{
    // Store domain neighbour as map so we can easily look for pair
    // with same face+proc.
    HashTable<label, labelPair, labelPairHash> map(domainFace.size());

    forAll(domainFace, bFaceI)
    {
        map.insert(labelPair(domainFace[bFaceI], domainProc[bFaceI]), bFaceI);
    }


    // Try to match mesh data.

    masterCoupledFaces.setSize(domainFace.size());
    slaveCoupledFaces.setSize(domainFace.size());
    label coupledI = 0;

    forAll(sourceFace, bFaceI)
    {
        if (sourceProc[bFaceI] != -1)
        {
            labelPair myData(sourceFace[bFaceI], sourceProc[bFaceI]);

            HashTable<label, labelPair, labelPairHash>::const_iterator iter =
                map.find(myData);

            if (iter != map.end())
            {
                label nbrBFaceI = iter();

                masterCoupledFaces[coupledI] = mesh.nInternalFaces() + bFaceI;
                slaveCoupledFaces[coupledI] =
                    domainMesh.nInternalFaces()
                  + nbrBFaceI;

                coupledI++;
            }
        }
    }

    masterCoupledFaces.setSize(coupledI);
    slaveCoupledFaces.setSize(coupledI);

    if (debug)
    {
        Pout<< "findCouples : found " << coupledI
            << " faces that will be stitched" << nl << endl;
    }
}


// Map data on boundary faces to new mesh (resulting from adding two meshes)
Foam::labelList Foam::fvMeshDistribute::mapBoundaryData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // mesh before adding
    const label nInternalFaces1,
    const labelList& boundaryData1  // added mesh
)
{
    labelList newBoundaryData(mesh.nFaces() - mesh.nInternalFaces());

    forAll(boundaryData0, oldBFaceI)
    {
        label newFaceI = map.oldFaceMap()[oldBFaceI + map.nOldInternalFaces()];

        // Face still exists (is nessecary?) and still boundary face
        if (newFaceI >= 0 && newFaceI >= mesh.nInternalFaces())
        {
            newBoundaryData[newFaceI - mesh.nInternalFaces()] =
                boundaryData0[oldBFaceI];
        }
    }

    forAll(boundaryData1, addedBFaceI)
    {
        label newFaceI = map.addedFaceMap()[addedBFaceI + nInternalFaces1];

        if (newFaceI >= 0 && newFaceI >= mesh.nInternalFaces())
        {
            newBoundaryData[newFaceI - mesh.nInternalFaces()] =
                boundaryData1[addedBFaceI];
        }
    }

    return newBoundaryData;
}


// Remove cells. Add all exposed faces to patch oldInternalPatchI
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::removeCells
(
    const labelList& cellsToRemove,
    const label oldInternalPatchI
)
{
    // Mesh change engine
    directPolyTopoChange meshMod(mesh_);

    // Cell removal topo engine. Do NOT synchronize parallel since
    // we are doing a local cell removal.
    directRemoveCells cellRemover(mesh_, false);

    // Get all exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    // Insert the topo changes
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        labelList(exposedFaces.size(), oldInternalPatchI),  // patch for exposed
                                                            // faces.
        meshMod
    );

    // Change the mesh. No inflation. Note: no parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, false);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        mesh_.clearOut();
    }

    return map;
}


// Delete and add processor patches. Changes mesh and returns per neighbour proc
// the processor patchID.
void Foam::fvMeshDistribute::addProcPatches
(
    const labelList& neighbourNewProc,   // processor that neighbour is on
    labelList& procPatchID              
)
{
    // Now use the neighbourFace/Proc to repatch the mesh. These two lists
    // contain for all current boundary faces the global patchID (for non-proc
    // patch) or the processor.

    labelList procPatchSizes(Pstream::nProcs(), 0);

    forAll(neighbourNewProc, bFaceI)
    {
        if (neighbourNewProc[bFaceI] != -1)
        {
            procPatchSizes[neighbourNewProc[bFaceI]]++;
        }
    }

    // Per neighbour processor the label of the processor patch
    procPatchID.setSize(Pstream::nProcs());

    forAll(procPatchSizes, procI)
    {
        if (procPatchSizes[procI] > 0)
        {
            const word patchName =
                "procBoundary"
              + name(Pstream::myProcNo())
              + "to"
              + name(procI);


            procPatchID[procI] = addProcPatch(patchName, procI);
            addPatchFields<scalar, volMesh>
            (
                processorFvPatchField<scalar>::typeName
            );
            addPatchFields<vector, volMesh>
            (
                processorFvPatchField<vector>::typeName
            );
            addPatchFields<tensor, volMesh>
            (
                processorFvPatchField<tensor>::typeName
            );

            addPatchFields<scalar, surfaceMesh>
            (
                processorFvPatchField<scalar>::typeName
            );
            addPatchFields<vector, surfaceMesh>
            (
                processorFvPatchField<vector>::typeName
            );
            addPatchFields<tensor, surfaceMesh>
            (
                processorFvPatchField<tensor>::typeName
            );
        }
        else
        {
            procPatchID[procI] = -1;
        }
    }
}


// Get boundary faces to be repatched. Is -1 or new patchID
Foam::labelList Foam::fvMeshDistribute::getProcBoundaryPatch
(
    const labelList& neighbourNewProc,  // new processor per boundary face
    const labelList& procPatchID        // patchID
)
{
    labelList patchIDs(neighbourNewProc);

    forAll(neighbourNewProc, bFaceI)
    {
        if (neighbourNewProc[bFaceI] != -1)
        {
            label nbrProc = neighbourNewProc[bFaceI];

            patchIDs[bFaceI] = procPatchID[nbrProc];
        }
        else
        {
            patchIDs[bFaceI] = -1;
        }
    }
    return patchIDs;
}


// Send mesh and coupling data.
void Foam::fvMeshDistribute::sendMesh
(
    const label domain,
    const fvMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourceNewProc
)
{
    if (debug)
    {
        Pout<< "Sending to domain " << domain << nl
            << "    nPoints:" << mesh.nPoints() << nl
            << "    nFaces:" << mesh.nFaces() << nl
            << "    nCells:" << mesh.nCells() << nl
            << "    nPatches:" << mesh.boundaryMesh().size() << nl
            << endl;
    }

    // Send
    OPstream toDomain(domain);
    toDomain
        << mesh.points()
        << mesh.faces()
        << mesh.allOwner()
        << mesh.allNeighbour()
        << mesh.boundaryMesh()
        << sourceFace
        << sourceProc
        << sourceNewProc;
}


// Receive mesh. Opposite of sendMesh
Foam::autoPtr<Foam::fvMesh> Foam::fvMeshDistribute::receiveMesh
(
    const label domain,
    Time& runTime,
    labelList& domainSourceFace,
    labelList& domainSourceProc,
    labelList& domainSourceNewProc
)
{
    IPstream fromNbr(domain);

    pointField domainPoints(fromNbr);
    faceList domainFaces(fromNbr);
    labelList domainAllOwner(fromNbr);
    labelList domainAllNeighbour(fromNbr);
    PtrList<entry> patchEntries(fromNbr);

    fromNbr
        >> domainSourceFace
        >> domainSourceProc
        >> domainSourceNewProc;

    // Construct fvMesh
    autoPtr<fvMesh> domainMeshPtr
    (
        new fvMesh 
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            domainPoints,
            domainFaces,
            domainAllOwner,
            domainAllNeighbour
        )
    );
    List<polyPatch*> patches(patchEntries.size());

    forAll(patchEntries, patchI)
    {
        patches[patchI] = polyPatch::New
        (
            patchEntries[patchI].keyword(),
            patchEntries[patchI].dict(),
            patchI,
            domainMeshPtr().boundaryMesh()
        ).ptr();
    }
    domainMeshPtr().addFvPatches(patches);

    return domainMeshPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fvMeshDistribute::fvMeshDistribute(fvMesh& mesh, const scalar mergeTol)
:
    mesh_(mesh),
    mergeTol_(mergeTol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvMeshDistribute::countCells
(
    const labelList& distribution
)
{
    labelList nCells(Pstream::nProcs(), 0);
    forAll(distribution, cellI)
    {
        label newProc = distribution[cellI];

        if (newProc < 0 || newProc >= Pstream::nProcs())
        {
            FatalErrorIn("fvMeshDistribute::distribute(const labelList&)")
                << "Distribution should be in range 0.." << Pstream::nProcs()-1
                << endl
                << "At index " << cellI << " distribution:" << newProc
                << abort(FatalError);
        }
        nCells[newProc]++;
    }
    return nCells;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::fvMeshDistribute::distribute
(
    const labelList& distribution
)
{
    // Some checks on distribution
    if (distribution.size() != mesh_.nCells())
    {
        FatalErrorIn("fvMeshDistribute::distribute(const labelList&)")
            << "Size of distribution:"
            << distribution.size() << " mesh nCells:" << mesh_.nCells()
            << abort(FatalError);
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Check all processors have same non-proc patches in same order.
    if (patches.checkParallelSync(true))
    {
        FatalErrorIn("fvMeshDistribute::distribute(const labelList&)")
            << "This application requires all non-processor patches"
            << " to be present in the same order on all patches" << nl
            << "followed by the processor patches (which of course are unique)."
            << nl
            << "Local patches:" << mesh_.boundaryMesh().names()
            << abort(FatalError);
    }

    // Save some data for mapping later on
    const label nOldPoints(mesh_.nPoints());
    const label nOldFaces(mesh_.nFaces());
    const label nOldCells(mesh_.nCells());
    labelList oldPatchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    forAll(patches, patchI)
    {
        oldPatchStarts[patchI] = patches[patchI].start();
        oldPatchNMeshPoints[patchI] = patches[patchI].nPoints();
    }



    // Short circuit trivial case.
    if (!Pstream::parRun())
    {
        // Collect all maps and return
        return autoPtr<mapDistributePolyMesh>
        (
            new mapDistributePolyMesh
            (
                mesh_,

                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts,
                oldPatchNMeshPoints,

                labelListList(1, identity(mesh_.nPoints())),//subPointMap
                labelListList(1, identity(mesh_.nFaces())), //subFaceMap
                labelListList(1, identity(mesh_.nCells())), //subCellMap
                labelListList(1, identity(patches.size())), //subPatchMap

                labelListList(1, identity(mesh_.nPoints())),//constructPointMap
                labelListList(1, identity(mesh_.nFaces())), //constructFaceMap
                labelListList(1, identity(mesh_.nCells())), //constructCellMap
                labelListList(1, identity(patches.size()))  //constructPatchMap
            )
        );
    }



    // Local environment of all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // A face is uniquely defined by
    //  - proc
    //  - local face no
    //
    // To glue the parts of meshes which can get sent from anywhere we
    // need to know on boundary faces what the above tuple on both sides is.
    // So we need to maintain:
    //  - original face
    //  - original processor id (= trivial)
    // For coupled boundaries (where the faces are 'duplicate') we take the
    // lowest numbered processor as the data to store.
    //
    // Additionally to create the procboundaries we need to know where the owner
    // cell on the other side now is: newNeighbourProc.
    //

    // physical boundary:
    //     sourceProc = -1
    //     sourceNewProc = -1
    //     sourceFace = patchID
    // coupled boundary:
    //     sourceProc = proc
    //     sourceNewProc = distribution[cell on proc]
    //     sourceFace = face
    labelList sourceFace;
    labelList sourceProc;
    labelList sourceNewProc;
    getNeighbourData(distribution, sourceFace, sourceProc, sourceNewProc);


    // Remove meshPhi. Since this would otherwise dissappear anyway
    // during topo changes and we have to guarantee that all the fields
    // can be sent.
    mesh_.clearOut();

    const wordList volScalars(mesh_.names(volScalarField::typeName));
    checkEqualWordList(volScalars);
    const wordList volVectors(mesh_.names(volVectorField::typeName));
    checkEqualWordList(volVectors);
    const wordList volTensors(mesh_.names(volTensorField::typeName));
    checkEqualWordList(volTensors);

    const wordList surfScalars(mesh_.names(surfaceScalarField::typeName));
    checkEqualWordList(surfScalars);
    const wordList surfVectors(mesh_.names(surfaceVectorField::typeName));
    checkEqualWordList(surfVectors);
    const wordList surfTensors(mesh_.names(surfaceTensorField::typeName));
    checkEqualWordList(surfTensors);


    // Find patch to temporarily put exposed and processor faces into.
    label oldInternalPatchI = findNonEmptyPatch();



    // Delete processor patches, starting from the back. Move all faces into
    // oldInternalPatchI.
    labelList repatchFaceMap;
    {
        autoPtr<mapPolyMesh> repatchMap = deleteProcPatches(oldInternalPatchI);

        // Store face map (only face ordering that changed)
        repatchFaceMap = repatchMap().faceMap();

        // Reorder all boundary face data (sourceProc, sourceFace etc.)
        labelList bFaceMap
        (
            SubList<label>
            (
                repatchMap().reverseFaceMap(),
                mesh_.nFaces() - mesh_.nInternalFaces(),
                mesh_.nInternalFaces()
            )
          - mesh_.nInternalFaces()
        );

        inplaceReorder(bFaceMap, sourceFace);
        inplaceReorder(bFaceMap, sourceProc);
        inplaceReorder(bFaceMap, sourceNewProc);
    }



    // Print a bit.
    if (debug)
    {
        Pout<< nl << "MESH WITH PROC PATCHES DELETED:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<scalar, volMesh>(mesh_);
        printFieldInfo<vector, volMesh>(mesh_);
        printFieldInfo<tensor, volMesh>(mesh_);
        printFieldInfo<scalar, surfaceMesh>(mesh_);
        printFieldInfo<vector, surfaceMesh>(mesh_);
        printFieldInfo<tensor, surfaceMesh>(mesh_);
        Pout<< nl << endl;
    }



    // Maps from subsetted mesh (that is sent) back to original maps
    labelListList subCellMap(Pstream::nProcs());
    labelListList subFaceMap(Pstream::nProcs());
    labelListList subPointMap(Pstream::nProcs());
    labelListList subPatchMap(Pstream::nProcs());
    // Maps from subsetted mesh to reconstructed mesh
    labelListList constructCellMap(Pstream::nProcs());
    labelListList constructFaceMap(Pstream::nProcs());
    labelListList constructPointMap(Pstream::nProcs());
    labelListList constructPatchMap(Pstream::nProcs());




    // Find out schedule
    // ~~~~~~~~~~~~~~~~~

    labelListList nSendCells(Pstream::nProcs());
    nSendCells[Pstream::myProcNo()] = countCells(distribution);
    Pstream::gatherList(nSendCells);
    Pstream::scatterList(nSendCells);

    // What to send to neighbouring domains
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(nSendCells[Pstream::myProcNo()], recvProc)
    {
        if
        (
            recvProc != Pstream::myProcNo()
         && nSendCells[Pstream::myProcNo()][recvProc] > 0
        )
        {
            // Send to recvProc

            if (debug)
            {
                Pout<< nl
                    << "SUBSETTING FOR DOMAIN " << recvProc
                    << " cells to send:"
                    << nSendCells[Pstream::myProcNo()][recvProc]
                    << nl << endl;
            }

            // Mesh subsetting engine
            fvMeshSubset subsetter(mesh_);

            // Subset the cells of the current domain.
            subsetter.setLargeCellSubset
            (
                distribution,
                recvProc,
                oldInternalPatchI,  // oldInternalFaces patch
                false               // no parallel sync
            );

            subCellMap[recvProc] = subsetter.cellMap();
            subFaceMap[recvProc] = renumber
            (
                repatchFaceMap,
                subsetter.faceMap()
            );
            subPointMap[recvProc] = subsetter.pointMap();
            subPatchMap[recvProc] = subsetter.patchMap();


            // Subset the boundary fields (owner/neighbour/processor)
            labelList procSourceFace;
            labelList procSourceProc;
            labelList procSourceNewProc;

            subsetBoundaryData
            (
                subsetter.subMesh(),
                subsetter.faceMap(),        // from subMesh to mesh
                subsetter.cellMap(),        //      ,,      ,,

                distribution,               // old mesh distribution
                mesh_.faceOwner(),          // old owner
                mesh_.faceNeighbour(),
                mesh_.nInternalFaces(),

                sourceFace,
                sourceProc,
                sourceNewProc,

                procSourceFace,
                procSourceProc,
                procSourceNewProc
            );

            // Send to neighbour
            sendMesh
            (
                recvProc,
                subsetter.subMesh(),
                procSourceFace,
                procSourceProc,
                procSourceNewProc
            );
            sendFields<scalar, volMesh>(recvProc, volScalars, subsetter);
            sendFields<vector, volMesh>(recvProc, volVectors, subsetter);
            sendFields<tensor, volMesh>(recvProc, volTensors, subsetter);

            sendFields<scalar, surfaceMesh>(recvProc, surfScalars, subsetter);
            sendFields<vector, surfaceMesh>(recvProc, surfVectors, subsetter);
            sendFields<tensor, surfaceMesh>(recvProc, surfTensors, subsetter);
        }
    }



    // Subset the part that stays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Save old mesh maps before changing mesh
        const labelList oldFaceOwner(mesh_.faceOwner());
        const labelList oldFaceNeighbour(mesh_.faceNeighbour());
        const label oldInternalFaces = mesh_.nInternalFaces();

        // Remove cells.
        autoPtr<mapPolyMesh> subMap
        (
            removeCells
            (
                select(false, distribution, Pstream::myProcNo()),
                oldInternalPatchI
            )
        );

        // Addressing from subsetted mesh
        subCellMap[Pstream::myProcNo()] = subMap().cellMap();
        subFaceMap[Pstream::myProcNo()] = renumber
        (
            repatchFaceMap,
            subMap().faceMap()
        );
        subPointMap[Pstream::myProcNo()] = subMap().pointMap();
        subPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Initialize all addressing into current mesh
        constructCellMap[Pstream::myProcNo()] = identity(mesh_.nCells());
        constructFaceMap[Pstream::myProcNo()] = identity(mesh_.nFaces());
        constructPointMap[Pstream::myProcNo()] = identity(mesh_.nPoints());
        constructPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Subset the mesh data: neighbourCell/neighbourProc
        // fields
        labelList domainSourceFace;
        labelList domainSourceProc;
        labelList domainSourceNewProc;

        subsetBoundaryData
        (
            mesh_,                          // new mesh
            subMap().faceMap(),             // from new to original mesh
            subMap().cellMap(),

            distribution,                   // distribution before subsetting
            oldFaceOwner,                   // owner before subsetting
            oldFaceNeighbour,               // neighbour        ,,
            oldInternalFaces,               // nInternalFaces   ,,

            sourceFace,
            sourceProc,
            sourceNewProc,

            domainSourceFace,
            domainSourceProc,
            domainSourceNewProc
        );

        sourceFace.transfer(domainSourceFace);
        sourceProc.transfer(domainSourceProc);
        sourceNewProc.transfer(domainSourceNewProc);
    }


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "STARTING MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<scalar, volMesh>(mesh_);
        printFieldInfo<vector, volMesh>(mesh_);
        printFieldInfo<tensor, volMesh>(mesh_);
        printFieldInfo<scalar, surfaceMesh>(mesh_);
        printFieldInfo<vector, surfaceMesh>(mesh_);
        printFieldInfo<tensor, surfaceMesh>(mesh_);
        Pout<< nl << endl;
    }



    // Receive and add what was sent
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(nSendCells, sendProc)
    {
        // Did processor sendProc send anything to me?
        if
        (
            sendProc != Pstream::myProcNo()
         && nSendCells[sendProc][Pstream::myProcNo()] > 0
        )
        {
            if (debug)
            {
                Pout<< nl
                    << "RECEIVING FROM DOMAIN " << sendProc
                    << " cells to receive:"
                    << nSendCells[sendProc][Pstream::myProcNo()]
                    << nl << endl;
            }

            // Receive from sendProc
            labelList domainSourceFace;
            labelList domainSourceProc;
            labelList domainSourceNewProc;

            // Opposite of sendMesh
            autoPtr<fvMesh> domainMeshPtr = receiveMesh
            (
                sendProc,
                const_cast<Time&>(mesh_.time()),
                domainSourceFace,
                domainSourceProc,
                domainSourceNewProc
            );
            fvMesh& domainMesh = domainMeshPtr();

            // Receive fields
            PtrList<volScalarField> vsf;
            receiveFields<scalar, volMesh>
            (
                sendProc,
                volScalars,
                domainMesh,
                vsf
            );

            PtrList<volVectorField> vvf;
            receiveFields<vector, volMesh>
            (
                sendProc,
                volVectors,
                domainMesh,
                vvf
            );
            PtrList<volTensorField> vtf;
            receiveFields<tensor, volMesh>
            (
                sendProc,
                volTensors,
                domainMesh,
                vtf
            );

            PtrList<surfaceScalarField> ssf;
            receiveFields<scalar, surfaceMesh>
            (
                sendProc,
                surfScalars,
                domainMesh,
                ssf
            );
            PtrList<surfaceVectorField> svf;
            receiveFields<vector, surfaceMesh>
            (
                sendProc,
                surfVectors,
                domainMesh,
                svf
            );
            PtrList<surfaceTensorField> stf;
            receiveFields<tensor, surfaceMesh>
            (
                sendProc,
                surfTensors,
                domainMesh,
                stf
            );


            constructCellMap[sendProc] = identity(domainMesh.nCells());
            constructFaceMap[sendProc] = identity(domainMesh.nFaces());
            constructPointMap[sendProc] = identity(domainMesh.nPoints());
            constructPatchMap[sendProc] =
                identity(domainMesh.boundaryMesh().size());


            // Print a bit.
            if (debug)
            {
                Pout<< nl << "RECEIVED MESH FROM:" << sendProc << endl;
                printMeshInfo(domainMesh);
                printFieldInfo<scalar, volMesh>(domainMesh);
                printFieldInfo<vector, volMesh>(domainMesh);
                printFieldInfo<tensor, volMesh>(domainMesh);
                printFieldInfo<scalar, surfaceMesh>(domainMesh);
                printFieldInfo<vector, surfaceMesh>(domainMesh);
                printFieldInfo<tensor, surfaceMesh>(domainMesh);
            }


            // Now this mesh we received (from sendProc) needs to be merged
            // with the current mesh. On the current mesh we have for all
            // boundaryfaces the original face and processor. See if we can
            // match these up to the received domainSourceFace and
            // domainSourceProc.
            labelList masterCoupledFaces;
            labelList slaveCoupledFaces;
            findCouples
            (
                mesh_,

                sourceFace,
                sourceProc,

                sendProc,
                domainMesh,
                domainSourceFace,
                domainSourceProc,

                masterCoupledFaces,
                slaveCoupledFaces
            );

            // Generate additional coupling info (points, edges) from
            // faces-that-match
            faceCoupleInfo couples
            (
                mesh_,
                masterCoupledFaces,
                domainMesh,
                slaveCoupledFaces,
                mergeTol_,              // merge tolerance
                true,                   // faces align
                true,                   // couples are ordered already
                false
            );


            // Add domainMesh to mesh
            // ~~~~~~~~~~~~~~~~~~~~~~

            autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
            (
                mesh_,
                domainMesh,
                couples,
                false           // no parallel comms
            );

            // Update mesh data: sourceFace,sourceProc for added
            // mesh.

            sourceFace = 
                mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceFace,
                    domainMesh.nInternalFaces(),
                    domainSourceFace
                );
            sourceProc = 
                mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceProc,
                    domainMesh.nInternalFaces(),
                    domainSourceProc
                );
            sourceNewProc = 
                mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceNewProc,
                    domainMesh.nInternalFaces(),
                    domainSourceNewProc
                );

            // Update all addressing so xxProcAddressing points to correct item
            // in masterMesh.
            const labelList& oldCellMap = map().oldCellMap();
            const labelList& oldFaceMap = map().oldFaceMap();
            const labelList& oldPointMap = map().oldPointMap();
            const labelList& oldPatchMap = map().oldPatchMap();

            forAll(constructPatchMap, procI)
            {
                if (procI != sendProc && constructPatchMap[procI].size() > 0)
                {
                    // Processor already in mesh (either myProcNo or received)
                    inplaceRenumber(oldCellMap, constructCellMap[procI]);
                    inplaceRenumber(oldFaceMap, constructFaceMap[procI]);
                    inplaceRenumber(oldPointMap, constructPointMap[procI]);
                    inplaceRenumber(oldPatchMap, constructPatchMap[procI]);
                }
            }

            // Added processor
            inplaceRenumber(map().addedCellMap(), constructCellMap[sendProc]);
            inplaceRenumber(map().addedFaceMap(), constructFaceMap[sendProc]);
            inplaceRenumber(map().addedPointMap(), constructPointMap[sendProc]);
            inplaceRenumber(map().addedPatchMap(), constructPatchMap[sendProc]);

            if (debug)
            {
                Pout<< nl << "MERGED MESH FROM:" << sendProc << endl;
                printMeshInfo(mesh_);
                printFieldInfo<scalar, volMesh>(mesh_);
                printFieldInfo<vector, volMesh>(mesh_);
                printFieldInfo<tensor, volMesh>(mesh_);
                printFieldInfo<scalar, surfaceMesh>(mesh_);
                printFieldInfo<vector, surfaceMesh>(mesh_);
                printFieldInfo<tensor, surfaceMesh>(mesh_);
                Pout<< nl << endl;
            }
        }
    }


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "REDISTRIBUTED MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<scalar, volMesh>(mesh_);
        printFieldInfo<vector, volMesh>(mesh_);
        printFieldInfo<tensor, volMesh>(mesh_);
        printFieldInfo<scalar, surfaceMesh>(mesh_);
        printFieldInfo<vector, surfaceMesh>(mesh_);
        printFieldInfo<tensor, surfaceMesh>(mesh_);
        Pout<< nl << endl;
    }



    // Add processorPatches
    // ~~~~~~~~~~~~~~~~~~~~

    // Per neighbour processor the patchID to it (or -1).
    labelList procPatchID;

    // Add processor patches.
    addProcPatches(sourceNewProc, procPatchID);

    // Put faces into correct patch. Note that we now have proper
    // processorPolyPatches again so repatching will take care of coupled face
    // ordering.

    // Get boundary faces to be repatched. Is -1 or new patchID
    labelList newPatchID
    (
        getProcBoundaryPatch
        (
            sourceNewProc,
            procPatchID
        )
    );

    // Change patches. Since this might change ordering of coupled faces
    // we also need to adapt our constructMaps.
    repatch(newPatchID, constructFaceMap);

    // Bit of hack: processorFvPatchField does not get reset since created
    // from nothing so explicitly reset.
    initPatchFields<scalar, volMesh>
    (
        processorFvPatchField<scalar>::typeName,
        pTraits<scalar>::zero
    );
    initPatchFields<vector, volMesh>
    (
        processorFvPatchField<vector>::typeName,
        pTraits<vector>::zero
    );
    initPatchFields<tensor, volMesh>
    (
        processorFvPatchField<tensor>::typeName,
        pTraits<tensor>::zero
    );
    initPatchFields<scalar, surfaceMesh>
    (
        processorFvPatchField<scalar>::typeName,
        pTraits<scalar>::zero
    );
    initPatchFields<vector, surfaceMesh>
    (
        processorFvPatchField<vector>::typeName,
        pTraits<vector>::zero
    );
    initPatchFields<tensor, surfaceMesh>
    (
        processorFvPatchField<tensor>::typeName,
        pTraits<tensor>::zero
    );


    mesh_.setInstance(mesh_.time().timeName());


    // Print a bit
    if (debug)
    {
        Pout<< nl << "FINAL MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<scalar, volMesh>(mesh_);
        printFieldInfo<vector, volMesh>(mesh_);
        printFieldInfo<tensor, volMesh>(mesh_);
        printFieldInfo<scalar, surfaceMesh>(mesh_);
        printFieldInfo<vector, surfaceMesh>(mesh_);
        printFieldInfo<tensor, surfaceMesh>(mesh_);
        Pout<< nl << endl;
    }

    // Collect all maps and return
    return autoPtr<mapDistributePolyMesh>
    (
        new mapDistributePolyMesh
        (
            mesh_,

            nOldPoints,
            nOldFaces,
            nOldCells,
            oldPatchStarts,
            oldPatchNMeshPoints,

            subPointMap,
            subFaceMap,
            subCellMap,
            subPatchMap,

            constructPointMap,
            constructFaceMap,
            constructCellMap,
            constructPatchMap,
            true                // reuse storage
        )
    );
}


// ************************************************************************* //
