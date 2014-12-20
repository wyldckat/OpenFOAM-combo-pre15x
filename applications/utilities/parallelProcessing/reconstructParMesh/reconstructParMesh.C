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

Application
    reconstructParMesh

Description
    Reconstructs a mesh using geometric information only. Writes
    point/face/cell procAddressing so afterwards reconstructPar can be used to
    reconstruct fields.

    Note:
    - parallelData is incorrect for shared cyclics.
    - uses geometric matching tolerance.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOobjectList.H"
#include "mergePolyMesh.H"
#include "polyTopoChanger.H"
#include "perfectInterface.H"
#include "treeBoundBox.H"
#include "matchPoints.H"
#include "perfectInterface.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "mapPolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tolerance to use for matching points. Is fraction of minimum bounding
// box dimension.
static const scalar matchTol = 1E-9;


// Filter out empty patches and converts processorPolyPatch ones into normal
// polyPatches.
void filterPatches(mergePolyMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    List<polyPatch*> newPatches(patches.size());
    label newPatchI = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.size() != 0)
        {
            if (pp.type() == processorPolyPatch::typeName)
            {
                Info<< "Converting processorPolyPatch " << pp.name()
                    << " to polyPatch" << endl;

                newPatches[newPatchI] = new polyPatch
                (
                    pp.name(),
                    pp.size(),
                    pp.start(),
                    newPatchI,
                    patches
                );

                newPatchI++;
            }
            else
            {
                newPatches[newPatchI++] = pp.clone(patches).ptr();
            }
        }
        else
        {
            Info<< "Removing zero sized patch " << pp.name() << endl;
        }
    }
    newPatches.setSize(newPatchI);

    mesh.removeBoundary();
    mesh.addPatches(newPatches);
}


int main(int argc, char *argv[])
{
    argList::noParallel();

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "This is an experimental tool which tries to merge individual"
        << " processor meshes back into one master mesh" << nl
        << "Use it if the original master mesh has been deleted or if the"
        << " processor meshes have been modified (topology change)" << nl
        << "This tool will write the resulting mesh to a new time step and"
        << " construct xxxxProcAddressing files in the processor meshes so"
        << " reconstructPar can be used to regenerate the fields on the master"
        << " mesh" << nl
        << "Not tested & use at your own risk" << nl << endl;

    int nProcs = 0;

    while
    (
        exists
        (
            args.rootPath()
          / args.caseName()
          / fileName(word("processor") + name(nProcs))
        )
    )
    {
        nProcs++;
    }

    Info<< "Found " << nProcs << " processor directories" << nl << endl;


    // Read all databases.
    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        Info<< "Reading database "
            << args.caseName()/fileName(word("processor") + name(procI))
            << endl;

        databases.hook
        (
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );
    }

    // Read processor0
    Info<< "Reading master mesh from " << databases[0].caseName()
        << " for time = " << databases[0].value()
        << endl;

    autoPtr<polyMesh> proc0MeshPtr
    (
        new polyMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                databases[0].timeName(),
                databases[0],
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );    

    mergePolyMesh masterMesh
    (
        runTime,
        proc0MeshPtr(),
        IOobject::AUTO_WRITE
    );

    proc0MeshPtr.reset(NULL);

    filterPatches(masterMesh);


    for (label procI = 1; procI < nProcs; procI++)
    {
        {
            Info<< "Reading mesh to add from "
                << databases[procI].caseName()
                << " for time = " << databases[procI].value()
                << endl;

            polyMesh meshToAdd
            (
                IOobject
                (
                    polyMesh::defaultRegion,
                    databases[procI].timeName(),
                    databases[procI]
                )
            );

            // Add elements to mesh
            Info<< "Adding to master mesh" << endl;
            masterMesh.addMesh(meshToAdd);
        }

        runTime++;

        masterMesh.merge();
        //masterMesh.resetMorph();
        Info<< "Added mesh from processor " << procI << endl;

        // Stitch this processor with all other processors already there
        for (label masterProcI = 0; masterProcI < procI; masterProcI++)
        {
            Info<< endl;

            // Construct some names
            word masterPatchName
            (
                "procBoundary"
              + Foam::name(masterProcI)
              + "to"
              + Foam::name(procI)
            );
            word slavePatchName
            (
                "procBoundary"
              + Foam::name(procI)
              + "to"
              + Foam::name(masterProcI)
            );


            label masterId =
                masterMesh.boundaryMesh().findPatchID(masterPatchName);

            label slaveId =
                masterMesh.boundaryMesh().findPatchID(slavePatchName);

            Info<< "Master patch name:" << masterPatchName
                << " at index:" << masterId << endl
                << "Slave patch name:" << slavePatchName
                << " at index:" << slaveId << endl;


            if
            (
                (masterId == -1 && slaveId != -1)
             || (masterId != -1 && slaveId == -1)
            )
            {
                FatalErrorIn(args.executable())
                    << "Problem : found only one of corresponding processor"
                    << " patches between processor " << masterProcI
                    << " and " << procI << endl
                    << "Master patch name:" << masterPatchName
                    << " at index:" << masterId << endl
                    << "Slave patch name:" << slavePatchName
                    << " at index:" << slaveId << abort(FatalError);
            }

            if (masterId != -1 && slaveId != -1)
            {
                Info<< "Trying to stitch master " << masterPatchName
                    << " and slave " << slavePatchName << endl;

                // Zone to put resulting internal faces in.
                word cutZoneName
                (
                    "cutFaceZone"
                  + Foam::name(masterProcI)
                  + "and"
                  + Foam::name(procI)
                );


                // Master patch
                const polyPatch& masterPatch =
                    masterMesh.boundaryMesh()[masterId];

                // Make list of masterPatch faces
                labelList isf(masterPatch.size());

                forAll (isf, i)
                {
                    isf[i] = masterPatch.start() + i;
                }


                List<polyMeshModifier*> tm(1);

                DynamicList<pointZone*> pz;
                DynamicList<faceZone*> fz;
                DynamicList<cellZone*> cz;

                // Add empty zone for resulting internal faces
                fz.append
                (
                    new faceZone
                    (
                        cutZoneName,
                        isf,
                        boolList(masterPatch.size(), false),
                        0,
                        masterMesh.faceZones()
                    )
                );

                // Note: make sure to add the zones BEFORE constructing
                // polyMeshModifier
                // (since looks up various zones at construction time)
                Info << "Adding point and face zones" << endl;

                masterMesh.removeZones();
                masterMesh.addZones(pz.shrink(), fz.shrink(), cz.shrink());

                polyTopoChanger coupler(masterMesh);
                coupler.setSize(1);

                // Add the perfect interface mesh modifier
                coupler.hook
                (
                    new perfectInterface
                    (
                        "couple",
                        0,
                        coupler,
                        cutZoneName,
                        masterPatchName,
                        slavePatchName
                    )
                );

                runTime++;

                // Execute all polyMeshModifiers
                Info << "Executing topology modifiers" << endl;
                coupler.changeMesh();

                masterMesh.removeZones();
            }
        }

        Info<< endl;
    }


    // Remove zero size patches
    filterPatches(masterMesh);


    runTime++;

    Info<< "\nWriting mesh to " << runTime.timeName() << endl;

    if (!masterMesh.polyMesh::write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


    // Determine addressing. Use geometric matching since I'm lazy.
    treeBoundBox bb(masterMesh.points());

    scalar typDim = matchTol * bb.minDim();

    Info<< "Bounding box of reconstructed mesh : " << bb << nl
        << "Minimum dimension of box           : " << bb.minDim() << nl
        << "Dimension used for matching        : " << typDim << nl
        << endl;


    // All processor to mastermesh
    PtrList<labelIOList> pointProcAddressing(nProcs);


    forAll(databases, procI)
    {
        Info<< "Reading processor " << procI << " mesh from "
            << databases[procI].caseName() << endl;

        polyMesh procMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                databases[procI].timeName(),
                databases[procI]
            )
        );

        // Match points
        {
            pointProcAddressing.hook
            (
                new labelIOList
                (
                    IOobject
                    (
                        "pointProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false                       // do not register
                    )
                )
            );

            bool matchOk =
                matchPoints
                (
                    procMesh.points(),
                    masterMesh.points(),
                    scalarField(procMesh.nPoints(), typDim),    // tolerance
                    false,                                      // verbose
                    pointProcAddressing[procI]
                );

            if (!matchOk)
            {
                FatalErrorIn(args.executable())
                    << "Points on processor " << procI
                    << " mesh do not match up "
                    << " with points on reconstructed mesh to within tolerance "
                    << typDim << exit(FatalError);
            }

            Info<< "Writing addressing from processor " << procI
                << " mesh points to merged mesh points to "
                << pointProcAddressing[procI].objectPath() << endl;

            pointProcAddressing[procI].write();
        }


        // Match cells
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                polyMesh::meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false                       // do not register
            )
        );

        bool matchOk =
            matchPoints
            (
                procMesh.cellCentres(),
                masterMesh.cellCentres(),
                scalarField(procMesh.nCells(), typDim),     // tolerance
                false,                                      // verbose
                cellProcAddressing
            );

        if (!matchOk)
        {
            FatalErrorIn(args.executable())
                << "Cellcentres on processor " << procI
                << " mesh do not match up "
                << " with points on reconstructed mesh to within tolerance "
                << typDim << exit(FatalError);
        }

        Info<< "Writing addressing from processor " << procI
            << " mesh cells to merged mesh cells to "
            << cellProcAddressing.objectPath() << endl;

        cellProcAddressing.write();


        // Match faces
        {
            labelIOList faceProcAddressing
            (
                IOobject
                (
                    "faceProcAddressing",
                    procMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false                       // do not register
                ),
                procMesh.nFaces()
            );

            bool matchOk =
                matchPoints
                (
                    procMesh.faceCentres(),
                    masterMesh.faceCentres(),
                    scalarField(procMesh.nFaces(), typDim),     // tolerance
                    false,                                      // verbose
                    faceProcAddressing
                );

            if (!matchOk)
            {
                FatalErrorIn(args.executable())
                    << "Facecentres on processor " << procI
                    << " mesh do not match up "
                    << " with points on reconstructed mesh to within tolerance "
                    << typDim << exit(FatalError);
            }

            // Now add turning index to faceProcAddressing.
            // See reconstrurPar for meaning of turning index.
            forAll(faceProcAddressing, procFaceI)
            {
                label masterFaceI = faceProcAddressing[procFaceI];

                if
                (
                   !procMesh.isInternalFace(procFaceI)
                 && masterMesh.isInternalFace(masterFaceI)
                )
                {
                    // proc face is now external but used to be internal face.
                    // Check if we have owner or neighbour.

                    label procOwn = procMesh.faceOwner()[procFaceI];
                    label masterOwn = masterMesh.faceOwner()[masterFaceI];

                    if (cellProcAddressing[procOwn] == masterOwn)
                    {
                        // No turning. Offset by 1.
                        faceProcAddressing[procFaceI]++;
                    }
                    else
                    {
                        // Turned face.
                        faceProcAddressing[procFaceI] =
                            -1 - faceProcAddressing[procFaceI];
                    }
                }
            }

            Info<< "Writing addressing from processor " << procI
                << " mesh faces to merged mesh faces to "
                << faceProcAddressing.objectPath() << endl;

            faceProcAddressing.write();
        }


        // Match boundaries
        {
            const polyBoundaryMesh& patches = procMesh.boundaryMesh();

            labelIOList boundaryProcAddressing
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    procMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false                   // Do not register
                ),
                labelList(patches.size(), -1)
            );

            const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();
            forAll(patches, patchI)
            {
                boundaryProcAddressing[patchI] =
                    masterPatches.findPatchID(patches[patchI].name());
            }

            Info<< "Writing addressing from processor " << procI
                << " mesh patches to merged mesh patches to "
                << boundaryProcAddressing.objectPath() << endl;

            boundaryProcAddressing.write();
        }
        Info<< nl << endl;
    }


    // Count number of times global points are used.
    labelList pointCount(masterMesh.nPoints(), 0);

    forAll(pointProcAddressing, procI)
    {
        const labelList& procToMaster = pointProcAddressing[procI];

        forAll(procToMaster, procPointI)
        {
            pointCount[procToMaster[procPointI]]++;
        }
    }

    // All points used by more than 2 processors are shared
    // Note: does not handle cyclics correctly.
    label nGlobalPoints = 0;

    forAll(pointCount, pointI)
    {
        if (pointCount[pointI] > 2)
        {
            nGlobalPoints++;
        }
    }

    Info<< "nGlobalPoints:" << nGlobalPoints << endl;

    // From mesh point to shared point index (or -1)
    labelList globalPointToShared(masterMesh.nPoints(), -1);

    nGlobalPoints = 0;

    forAll(pointCount, pointI)
    {
        if (pointCount[pointI] > 2)
        {
            globalPointToShared[pointI] = nGlobalPoints++;
        }
    }


    // Construct shared point info per processor.
    forAll(pointProcAddressing, procI)
    {
        const labelList& procToMaster = pointProcAddressing[procI];

        DynamicList<label> localSharedPoints(nGlobalPoints);
        DynamicList<label> localSharedPointAddr(nGlobalPoints);
        DynamicList<label> localSharedPointGlobalLabels(nGlobalPoints);

        forAll(procToMaster, procPointI)
        {
            label masterPointI = procToMaster[procPointI];

            if (pointCount[masterPointI] > 2)
            {
                // Shared points
                localSharedPoints.append(procPointI);
                localSharedPointAddr.append(globalPointToShared[masterPointI]);
                localSharedPointGlobalLabels.append(masterPointI);
            }
        }

        localSharedPoints.shrink();
        localSharedPointAddr.shrink();
        localSharedPointGlobalLabels.shrink();


        //Hack: write without constructing parallelInfo
        fileName facesInst
        (
            databases[procI].findInstance
            (
                polyMesh::meshSubDir, "faces"
            )
        );

        IOdictionary dict
        (
            IOobject
            (
                "parallelData",
                facesInst,
                polyMesh::meshSubDir,
                databases[procI],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false                   // Do not register
            )
        );

        dict.add("cyclicParallel", false);
        dict.add("nTotalPoints", masterMesh.nPoints());
        dict.add("nTotalFaces", masterMesh.nFaces());
        dict.add("nTotalCells", masterMesh.nCells());

        dict.add("nGlobalPoints", nGlobalPoints);
        dict.add("sharedPointLabels", localSharedPoints);
        dict.add("sharedPointAddr", localSharedPointAddr);
        dict.add("sharedPointGlobalLabels", localSharedPointGlobalLabels);


        Info<< "Writing addressing from processor " << procI
            << " mesh shared points to merged mesh points to "
            << dict.instance()
             / dict.name()
            << endl;


        dict.regIOobject::write
        (
            databases[procI].writeFormat(),
            IOstream::currentVersion,
            databases[procI].writeCompression()
        );
    }

    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
