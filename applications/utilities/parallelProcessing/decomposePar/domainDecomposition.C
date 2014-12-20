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
    Automatic domain decomposition class for FOAM meshes

\*---------------------------------------------------------------------------*/

#include "domainDecomposition.H"
#include "Time.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "parallelInfo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
domainDecomposition::domainDecomposition(const IOobject& io)
:
    fvMesh(io),
    decompositionDict_
    (
        IOobject
        (
            "decomposeParDict",
            time().system(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nProcs_(readInt(decompositionDict_.lookup("numberOfSubdomains"))),
    cellToProc_(nCells()),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procCellAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    globallySharedPoints_(0),
    cyclicParallel_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

domainDecomposition::~domainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool domainDecomposition::writeDecomposition()
{
    Info<< "\nConstructing processor meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(2*globallySharedPoints_.size());

    forAll (globallySharedPoints_, pointI)
    {
        sharedPointLookup.insert(globallySharedPoints_[pointI], pointI);
    }

    // Write out the meshes
    for (label procI = 0; procI < nProcs_; procI++)
    {
        // Create processor points
        const labelList& curPointLabels = procPointAddressing_[procI];

        const pointField& meshPoints = points();

        labelList pointLookup(nPoints(), -1);

        pointField procPoints(curPointLabels.size());

        forAll (curPointLabels, pointI)
        {
            procPoints[pointI] = meshPoints[curPointLabels[pointI]];

            pointLookup[curPointLabels[pointI]] = pointI;
        }

        // Create processor faces
        const labelList& curFaceLabels = procFaceAddressing_[procI];

        const faceList& meshFaces = faces();

        labelList faceLookup(nFaces(), -1);

        faceList procFaces(curFaceLabels.size());

        forAll (curFaceLabels, faceI)
        {
            // Mark the original face as used
            // Remember to decrement the index by one (turning index)
            // 
            label curF = mag(curFaceLabels[faceI]) - 1;

            faceLookup[curF] = faceI;

            // get the original face
            labelList origFaceLabels;

            if (curFaceLabels[faceI] >= 0)
            {
                // face not turned
                origFaceLabels = meshFaces[curF];
            }
            else
            {
                origFaceLabels = meshFaces[curF].reverseFace();
            }

            // translate face labels into local point list
            face& procFaceLabels = procFaces[faceI];

            procFaceLabels.setSize(origFaceLabels.size());

            forAll (origFaceLabels, pointI)
            {
                procFaceLabels[pointI] = pointLookup[origFaceLabels[pointI]];
            }
        }

        // Create processor cells
        const labelList& curCellLabels = procCellAddressing_[procI];

        const cellList& meshCells = cells();

        cellList procCells(curCellLabels.size());

        forAll (curCellLabels, cellI)
        {
            const labelList& origCellLabels = meshCells[curCellLabels[cellI]];

            cell& curCell = procCells[cellI];

            curCell.setSize(origCellLabels.size());

            forAll (origCellLabels, cellFaceI)
            {
                curCell[cellFaceI] = faceLookup[origCellLabels[cellFaceI]];
            }
        }

        // Create processor mesh without a boundary

        fileName processorCasePath
        (
            time().caseName()/fileName(word("processor") + Foam::name(procI))
        );

        // make the processor directory
        mkDir(time().rootPath()/processorCasePath);

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            processorCasePath
        );

        // create the mesh
        polyMesh procMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                "constant",
                processorDb
            ),
            procPoints,
            procFaces,
            procCells
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

        const labelList& curPatchSizes = procPatchSize_[procI];

        const labelList& curPatchStarts = procPatchStartIndex_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        const labelList& curProcessorPatchStarts =
            procProcessorPatchStartIndex_[procI];

        const polyPatchList& meshPatches = boundaryMesh();

        List<polyPatch*> procPatches
        (
            curPatchSizes.size()
          + curProcessorPatchSizes.size(),
            reinterpret_cast<polyPatch*>(NULL)
        );

        label nPatches = 0;

        forAll (curPatchSizes, patchI)
        {
            procPatches[nPatches] =
                meshPatches[curBoundaryAddressing[patchI]].clone
                (
                    procMesh.boundaryMesh(),
                    nPatches,
                    curPatchSizes[patchI],
                    curPatchStarts[patchI]
                ).ptr();

            nPatches++;
        }

        forAll (curProcessorPatchSizes, procPatchI)
        {
            procPatches[nPatches] =
                new processorPolyPatch
                (
                    word("procBoundary") + Foam::name(procI)
                  + word("to")
                  + Foam::name(curNeighbourProcessors[procPatchI]),
                    curProcessorPatchSizes[procPatchI],
                    curProcessorPatchStarts[procPatchI],
                    nPatches,
                    procMesh.boundaryMesh(),
                    procI,
                    curNeighbourProcessors[procPatchI]
            );

            nPatches++;
        }

        // Add boundary patches
        procMesh.addPatches(procPatches);

        // Create and add zones

        // Point zones

        List<pointZone*> procPointZones(0);

        const pointZoneMesh& pz = pointZones();

        labelList createdPointZones(pz.size());
        List<SLList<label> > zonePoints(pz.size());
        label nCreatedPointZones = 0;

        // If there are zoned elements, map them zone by zone
        if (pz.zoneMap().size() > 0)
        {
            // Go through all the zoned points and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll (curPointLabels, pointI)
            {
                label zoneForPoint = pz.whichZone(curPointLabels[pointI]);

                if (zoneForPoint >= 0)
                {
                    // Point belongs to a zone.  Find out if this zone exists
                    // and if so add the point.
                    bool foundZone = false;

                    for (label zoneI = 0; zoneI < nCreatedPointZones; zoneI++)
                    {
                        if (createdPointZones[zoneI] == zoneForPoint)
                        {
                            // Zone exists.  Add the face
                            zonePoints[zoneI].append(pointI);

                            foundZone = true;
                            break;
                        }
                    }

                    if (!foundZone)
                    {
                        // Zone not created yet.  Add it
                        createdPointZones[nCreatedPointZones] = zoneForPoint;
                        zonePoints[nCreatedPointZones].append(pointI);

                        nCreatedPointZones++;
                    }
                }
            }

            // Reset the sizes of lists
            createdPointZones.setSize(nCreatedPointZones);
            zonePoints.setSize(nCreatedPointZones);

            // Hook the point zones
            procPointZones.setSize(createdPointZones.size());

            forAll (createdPointZones, zoneI)
            {
                procPointZones[zoneI] =
                    pz[createdPointZones[zoneI]].clone
                    (
                        procMesh.pointZones(),
                        zoneI,
                        zonePoints[zoneI]
                    ).ptr();
            }
        }


        // Do face zones

        List<faceZone*> procFaceZones(0);

        const faceZoneMesh& fz = faceZones();

        labelList createdFaceZones(fz.size());
        List<SLList<label> > zoneFaces(fz.size());
        List<SLList<bool> > zoneFaceFlips(fz.size());
        label nCreatedFaceZones = 0;

        // If there are zoned elements, map them zone by zone
        if (fz.zoneMap().size() > 0)
        {
            // Go through all the zoned faces and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll (curFaceLabels, faceI)
            {
                // Remember to decrement the index by one (turning index)
                // 
                label curF = mag(curFaceLabels[faceI]) - 1;

                label zoneForFace = fz.whichZone(curF);

                if (zoneForFace >= 0)
                {
                    // Face belongs to a zone.  Find out if this zone exists
                    // and if so add the face.
                    bool foundZone = false;

                    for (label zoneI = 0; zoneI < nCreatedFaceZones; zoneI++)
                    {
                        if (createdFaceZones[zoneI] == zoneForFace)
                        {
                            // Zone exists.  Add the face
                            zoneFaces[zoneI].append(faceI);

                            bool flip = 
                                fz[zoneForFace].flipMap()
                                [
                                    fz[zoneForFace].whichFace(curF)
                                ];

                            if (curFaceLabels[faceI] < 0)
                            {
                                flip = !flip;
                            }

                            zoneFaceFlips[zoneI].append(flip);

                            foundZone = true;
                            break;
                        }
                    }

                    if (!foundZone)
                    {
                        // Zone not created yet.  Add it
                        createdFaceZones[nCreatedFaceZones] = zoneForFace;
                        zoneFaces[nCreatedFaceZones].append(faceI);

                        bool flip = 
                            fz[zoneForFace].flipMap()
                            [
                                fz[zoneForFace].whichFace(curF)
                            ];

                        if (curFaceLabels[faceI] < 0)
                        {
                            flip = !flip;
                        }

                        zoneFaceFlips[nCreatedFaceZones].append(flip);

                        nCreatedFaceZones++;
                    }
                }
            }

            // Reset the sizes of lists
            createdFaceZones.setSize(nCreatedFaceZones);
            zoneFaces.setSize(nCreatedFaceZones);
            zoneFaceFlips.setSize(nCreatedFaceZones);

            // Hook the face zones
            procFaceZones.setSize(createdFaceZones.size());

            forAll (createdFaceZones, zoneI)
            {
                procFaceZones[zoneI] =
                    fz[createdFaceZones[zoneI]].clone
                    (
                        zoneFaces[zoneI],
                        zoneFaceFlips[zoneI],
                        zoneI,
                        procMesh.faceZones()
                    ).ptr();
            }
        }

        // Do cell zones

        List<cellZone*> procCellZones(0);

        const cellZoneMesh& cz = cellZones();

        labelList createdCellZones(cz.size());
        List<SLList<label> > zoneCells(cz.size());
        label nCreatedCellZones = 0;

        // If there are zoned elements, map them zone by zone
        if (cz.zoneMap().size() > 0)
        {
            // Go through all the zoned cells and find out if they
            // belong to a zone.  If so, add it to the zone as
            // necessary
            forAll (curCellLabels, cellI)
            {
                label zoneForCell = cz.whichZone(curCellLabels[cellI]);

                if (zoneForCell >= 0)
                {
                    // Cell belongs to a zone.  Find out if this zone exists
                    // and if so add the cell.
                    bool foundZone = false;

                    for (label zoneI = 0; zoneI < nCreatedCellZones; zoneI++)
                    {
                        if (createdCellZones[zoneI] == zoneForCell)
                        {
                            // Zone exists.  Add the cell
                            zoneCells[zoneI].append(cellI);

                            foundZone = true;
                            break;
                        }
                    }

                    if (!foundZone)
                    {
                        // Zone not created yet.  Add it
                        createdCellZones[nCreatedCellZones] = zoneForCell;
                        zoneCells[nCreatedCellZones].append(cellI);

                        nCreatedCellZones++;
                    }
                }
            }

            // Reset the sizes of lists
            createdCellZones.setSize(nCreatedCellZones);
            zoneCells.setSize(nCreatedCellZones);

            // Hook the cell zones
            procCellZones.setSize(createdCellZones.size());

            forAll (createdCellZones, zoneI)
            {
                procCellZones[zoneI] =
                    cz[createdCellZones[zoneI]].clone
                    (
                        zoneCells[zoneI],
                        zoneI,
                        procMesh.cellZones()
                    ).ptr();
            }
        }


        // Hook the zones onto the mesh
        procMesh.addZones(procPointZones, procFaceZones, procCellZones);

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        procMesh.write();

        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << procMesh.nCells()
            << endl;

        label nBoundaryFaces = 0;

        forAll (procMesh.boundaryMesh(), patchi)
        {
            if
            (
                procMesh.boundaryMesh()[patchi].type()
             == processorPolyPatch::typeName
            )
            {
                const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>
                (
                    procMesh.boundaryMesh()[patchi]
                );

                Info<< "    Number of faces shared with processor " 
                    << ppp.neighbProcNo() << " = " << ppp.size() << endl;
            }
            else
            {
                nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
            }
        }

        Info<< "    Number of boundary faces = " << nBoundaryFaces << endl;


        // create and write the addressing information
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procPointAddressing_[procI]
        );
        pointProcAddressing.write();

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procFaceAddressing_[procI]
        );
        faceProcAddressing.write();

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procCellAddressing_[procI]
        );
        cellProcAddressing.write();

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procBoundaryAddressing_[procI]
        );
        boundaryProcAddressing.write();
    }

    return true;
}


// ************************************************************************* //
