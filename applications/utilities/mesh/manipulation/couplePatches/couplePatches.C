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
    Utility to reorder cyclic and processor patches.

    Uses dummy morph to sort things out. Note that only geometricMatch will
    work since topological match would require an unchanged face which
    we have no information about.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "morphMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "coupledPolyPatch.H"
#include "processorPolyPatch.H"
#include "PstreamReduceOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Using geometry to calculate face correspondence across"
        << " coupled boundaries (processor, cyclic)" << nl
        << "This will only work for cyclics if they are parallel or"
        << " their rotation is defined across the origin" << nl
        << endl;

    coupledPolyPatch::setGeometricMatch(true);


    Info<< "Create morphMesh for time = " << runTime.value() << nl << endl;

    morphMesh mesh
    (
        IOobject
        (
            morphMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );


    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    bool hasCoupled = false;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            hasCoupled = true;

            break;
        }
    }

    reduce(hasCoupled, orOp<bool>());

    if (hasCoupled)
    {
        Info<< "Mesh has coupled patches ..." << nl << endl;

        Info<< "Testing for correct face ordering ..." << endl;
        bool wrongOrder = mesh.boundaryMesh().checkFaceOrder(true);

        if (wrongOrder)
        {
            // Dummy topo changes container
            polyTopoChange meshMod(mesh);

            // Do all changes
            Info<< "Doing dummy mesh morph to correct face ordering ..."
                << endl;

            runTime++;

            mesh.setMorphTimeIndex(runTime.timeIndex());

            mesh.updateTopology(meshMod);

            // Move mesh (since morphing does not do this)
            mesh.movePoints(mesh.morphMap().preMotionPoints());

            // Set the precision of the points data to 10
            IOstream::defaultPrecision(10);

            // Write resulting mesh
            Info << "Writing morphed mesh to time " << runTime.value() << endl;

            mesh.write();
        }
        else
        {
            Info<< "Coupled patch face ordering ok. Nothing changed ..." << nl
                << endl;
        }
    }
    else
    {
        Info<< "Mesh has no coupled patches. Nothing changed ..." << nl << endl;
    }
        

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
