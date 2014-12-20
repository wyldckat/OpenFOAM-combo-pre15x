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
    'Stitches' a mesh.

    Takes a mesh and two patches and merges the faces on the two patches
    (if geometrically possible) so the faces become internal.

    Can do
    - 'perfect' match: faces and points on patches align exactly. Order might
    be different though.
    - 'integral' match: where the surfaces on both patches exactly
    match but the individual faces not
    - 'partial' match: where the non-overlapping part of the surface remains
    in the respective patch.

    Note : Is just a front-end to perfectInterface/slidingInterface.

    Comparable to running a meshModifier of the form
    (if masterPatch is called "M" and slavePatch "S"):

    couple
    {
        type                    slidingInterface;
        masterFaceZoneName      MSMasterZone
        slaveFaceZoneName       MSSlaveZone
        cutPointZoneName        MSCutPointZone
        cutFaceZoneName         MSCutFaceZone
        masterPatchName         M;
        slavePatchName          S;
        typeOfMatch             partial or integral
    }


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "morphMesh.H"
#include "mapPolyMesh.H"
#include "Map.H"
#include "ListSearch.H"
#include "morphMesh.H"
#include "IndirectList.H"
#include "slidingInterface.H"
#include "perfectInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Checks whether patch present
void checkPatch(const polyBoundaryMesh& bMesh, const word& name)
{
    label patchI = bMesh.findPatchID(name);

    if (patchI == -1)
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Cannot find patch " << name << endl
            << "It should be present and of non-zero size" << endl
            << "Valid patches are " << bMesh.names()
            << exit(FatalError);
    }

    if (bMesh[patchI].size() == 0)
    {
        FatalErrorIn("checkPatch(const polyBoundaryMesh&, const word&)")
            << "Patch " << name << " is present but zero size"
            << exit(FatalError);
    }
}


// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::noParallel();

    Foam::argList::validArgs.append("masterPatch");
    Foam::argList::validArgs.append("slavePatch");

    Foam::argList::validOptions.insert("partial", "");
    Foam::argList::validOptions.insert("perfect", "");

#   include "setRootCase.H"
#   include "createTime.H"

    word masterPatchName(args.args()[3]);
    word slavePatchName(args.args()[4]);

    bool partialCover = args.options().found("partial");
    bool perfectCover = args.options().found("perfect");

    if (partialCover && perfectCover)
    {
        FatalErrorIn(args.executable())
            << "Cannot both supply partial and perfect." << endl
            << "Use perfect match option if the patches perfectly align"
            << " (both vertex positions and face centres)" << endl
            << exit(FatalError);
    }


    const word mergePatchName(masterPatchName + slavePatchName);
    const word cutZoneName(mergePatchName + "CutFaceZone");

    slidingInterface::typeOfMatch tom = slidingInterface::INTEGRAL;

    if (partialCover)
    {
        Info<< "Coupling partially overlapping patches "
            << masterPatchName << " and " << slavePatchName << nl
            << "Resulting internal faces will be in faceZone " << cutZoneName
            << nl
            << "Any uncovered faces will remain in their patch"
            << endl;

        tom = slidingInterface::PARTIAL;
    }
    else if (perfectCover)
    {
        Info<< "Coupling perfectly aligned patches "
            << masterPatchName << " and " << slavePatchName << nl
            << "Resulting (internal) faces will be in faceZone " << cutZoneName
            << nl << nl
            << "Note: both patches need to align perfectly." << nl
            << "Both the vertex"
            << " positions and the face centres need to align to within" << nl
            << "a tolerance given by the minimum edge length on the patch"
            << endl;
    }
    else
    {
        Info<< "Coupling patches " << masterPatchName << " and "
            << slavePatchName << nl
            << "Resulting (internal) faces will be in faceZone " << cutZoneName
            << nl << nl
            << "Note: the overall area covered by both patches should be"
            << " identical (\"integral\" interface)." << endl
            << "If this is not the case use the -partial option" << nl << endl;
    }


    morphMesh pMesh
    (
        IOobject
        (
            morphMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );


    // Check for non-empty master and slave patches
    checkPatch(pMesh.boundaryMesh(), masterPatchName);
    checkPatch(pMesh.boundaryMesh(), slavePatchName);

    // Create and add face zones and mesh modifiers

    // Master patch
    const polyPatch& masterPatch =
        pMesh.boundaryMesh()
        [
            pMesh.boundaryMesh().findPatchID(masterPatchName)
        ];

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

    if (perfectCover)
    {
        // Add empty zone for resulting internal faces
        fz.append
        (
            new faceZone
            (
                cutZoneName,
                isf,
                boolList(masterPatch.size(), false),
                0,
                pMesh.faceZones()
            )
        );

        // Note: make sure to add the zones BEFORE constructing polyMeshModifier
        // (since looks up various zones at construction time)
        Info << "Adding point and face zones" << endl;
        pMesh.addZones(pz.shrink(), fz.shrink(), cz.shrink());

        // Add the perfect interface mesh modifier
        tm[0] =
            new perfectInterface
            (
                "couple",
                0,
                pMesh.morphEngine(),
                cutZoneName,
                masterPatchName,
                slavePatchName
            );
    }
    else
    {
        pz.append
        (
            new pointZone
            (
                mergePatchName + "CutPointZone",
                labelList(0),
                0,
                pMesh.pointZones()
            )
        );

        fz.append
        (
            new faceZone
            (
                mergePatchName + "MasterZone",
                isf,
                boolList(masterPatch.size(), false),
                0,
                pMesh.faceZones()
            )
        );

        // Slave patch
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

        fz.append
        (
            new faceZone
            (
                mergePatchName + "SlaveZone",
                osf,
                boolList(slavePatch.size(), false),
                1,
                pMesh.faceZones()
            )
        );

        // Add empty zone for cut faces
        fz.append
        (
            new faceZone
            (
                cutZoneName,
                labelList(0),
                boolList(0, false),
                2,
                pMesh.faceZones()
            )
        );


        // Note: make sure to add the zones BEFORE constructing polyMeshModifier
        // (since looks up various zones at construction time)
        Info << "Adding point and face zones" << endl;
        pMesh.addZones(pz.shrink(), fz.shrink(), cz.shrink());

        // Add the sliding interface mesh modifier
        tm[0] =
            new slidingInterface
            (
                "couple",
                0,
                pMesh.morphEngine(),
                mergePatchName + "MasterZone",
                mergePatchName + "SlaveZone",
                mergePatchName + "CutPointZone",
                cutZoneName,
                masterPatchName,
                slavePatchName,
                tom                   // integral or partial
            );
    }


    Info << "Adding topology modifiers" << endl;
    pMesh.addTopologyModifiers(tm);

    pMesh.morphEngine().write();

    runTime++;

    // Execute all polyMeshModifiers
    pMesh.polyMesh::updateTopology();

    pMesh.movePoints(pMesh.morphMap().preMotionPoints());

    Info << nl << "Writing polyMesh to time " << runTime.timeName() << endl;

    IOstream::defaultPrecision(10);
    if (!pMesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    pMesh.faceZones().write();
    pMesh.pointZones().write();
    pMesh.cellZones().write();

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
