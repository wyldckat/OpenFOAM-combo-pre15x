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
    A multi-block mesh generator.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyMesh.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
#   include "checkOptions.H"

    Info<< nl << "Reading block mesh description dictionary" << endl;

    IOobject meshDescriptionIOobject
    (
        "blockMeshDict",
        runTime.constant(),
        "polyMesh",
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!meshDescriptionIOobject.headerOk())
    {
        meshDescriptionIOobject = IOobject
        (
            "meshDescription",
            runTime.constant(),
            "polyMesh",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
    }

    if (!meshDescriptionIOobject.headerOk())
    {
        meshDescriptionIOobject = IOobject
        (
            "meshDescription",
            runTime.constant(),
            "mesh",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
    }

    if (!meshDescriptionIOobject.headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot find mesh description file " << nl
            << runTime.constant()/"polyMesh"/"blockMeshDict" << " or " << nl
            << runTime.constant()/"polyMesh"/"meshDescription" << " or " << nl
            << runTime.constant()/"mesh"/"meshDescription"
            << exit(FatalError);
    }

    IOdictionary meshDescription(meshDescriptionIOobject);

    Info<< nl << "Creating block mesh" << endl;

    blockMesh blocks(meshDescription);


    if (writeTopo)
    {
        word objMeshFile("blockTopology.obj");

        Info<< nl << "Dumping block structure as Lightwave obj format"
            << " to " << objMeshFile << endl;

        // Write mesh as edges.

        OFstream objStream(objMeshFile);

        blocks.writeTopology(objStream);

        Info<< nl << "end" << endl;

        return 0;
    }



    Info<< nl << "Creating mesh from block mesh" << endl;

    wordList patchNames = blocks.patchNames();
    wordList patchTypes = blocks.patchTypes();
    word defaultFacesType = emptyPolyPatch::typeName;
    wordList patchPhysicalTypes = blocks.patchPhysicalTypes();

    preservePatchTypes
    (
        runTime,
        runTime.constant(),
        polyMesh::meshSubDir,
        patchNames,
        patchTypes,
        defaultFacesType,
        patchPhysicalTypes
    );

    attachPolyMesh pMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        blocks.points(),
        blocks.cells(),
        blocks.patches(),
        patchNames,
        patchTypes,
        defaultFacesType,
        patchPhysicalTypes
    );

    // Read in a list of dictionaries for the merge patch pairs
    if (meshDescription.found("mergePatchPairs"))
    {
        Info<< nl << "Creating merge patch pairs" << nl << endl;

        List<Pair<word> > mergePatchPairs
        (
            meshDescription.lookup("mergePatchPairs")
        );

        if (mergePatchPairs.size() > 0)
        {
            // Create and add point and face zones and mesh modifiers
            List<pointZone*> pz(mergePatchPairs.size());
            List<faceZone*> fz(3*mergePatchPairs.size());
            List<cellZone*> cz(0);

            forAll (mergePatchPairs, pairI)
            {
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                pz[pairI] =
                    new pointZone
                    (
                        mergeName + "CutPointZone",
                        labelList(0),
                        0,
                        pMesh.pointZones()
                    );

                // Master patch
                const word masterPatchName(mergePatchPairs[pairI].first());
                const polyPatch& masterPatch =
                    pMesh.boundaryMesh()
                    [
                        pMesh.boundaryMesh().findPatchID(masterPatchName)
                    ];

                labelList isf(masterPatch.size());

                forAll (isf, i)
                {
                    isf[i] = masterPatch.start() + i;
                }

                fz[3*pairI] =
                    new faceZone
                    (
                        mergeName + "MasterZone",
                        isf,
                        boolList(masterPatch.size(), false),
                        0,
                        pMesh.faceZones()
                    );

                // Slave patch
                const word slavePatchName(mergePatchPairs[pairI].second());
                const polyPatch& slavePatch =
                    pMesh.boundaryMesh()
                    [
                        pMesh.boundaryMesh().findPatchID(slavePatchName)
                    ];

                labelList osf(slavePatch.size());

                forAll (osf, i)
                {
                    osf[i] = slavePatch.start() + i;
                }

                fz[3*pairI + 1] =
                    new faceZone
                    (
                        mergeName + "SlaveZone",
                        osf,
                        boolList(slavePatch.size(), false),
                        1,
                        pMesh.faceZones()
                    );

                // Add empty zone for cut faces
                fz[3*pairI + 2] =
                    new faceZone
                    (
                        mergeName + "CutFaceZone",
                        labelList(0),
                        boolList(0, false),
                        2,
                        pMesh.faceZones()
                    );
            }  // end of all merge pairs

            Info << "Adding point and face zones" << endl;
            pMesh.addZones(pz, fz, cz);

            List<polyMeshModifier*> tm(mergePatchPairs.size());

            forAll (mergePatchPairs, pairI)
            {
                const word mergeName
                (
                    mergePatchPairs[pairI].first()
                  + mergePatchPairs[pairI].second()
                  + name(pairI)
                );

                // Add the sliding interface mesh modifier
                tm[pairI] =
                    new slidingInterface
                    (
                        "couple" + name(pairI),
                        pairI,
                        pMesh,
                        mergeName + "MasterZone",
                        mergeName + "SlaveZone",
                        mergeName + "CutPointZone",
                        mergeName + "CutFaceZone",
                        mergePatchPairs[pairI].first(),
                        mergePatchPairs[pairI].second(),
                        slidingInterface::INTEGRAL, // always integral
                        intersection::VISIBLE
                    );
            }

            Info << "Adding topology modifiers" << endl;
            pMesh.addTopologyModifiers(tm);

            pMesh.attach();
        }
    }
    else
    {
        Info<< nl << "There are no merge patch pairs edges" << endl;
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << nl << "Writing polyMesh" << endl;
    if (!pMesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    Info<< nl << "end" << endl;
    return 0;
}


// ************************************************************************* //
