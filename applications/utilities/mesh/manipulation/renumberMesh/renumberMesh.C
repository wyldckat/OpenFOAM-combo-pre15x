/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anispulation  |
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

Application
    renumberMesh

Description
    Renumbers the cell list in order to reduce the bandwidth, reading and
    renumbering all fields from all the time directories.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "fvMeshBandCompression.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "setRootCase.H"
#   include "createTime.H"

    fvMeshBandCompression mesh
    (
        IOobject
        (
            fvMeshBandCompression::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    const unallocLabelList& oldFownerC = mesh.faceOwner();
    const unallocLabelList& oldFneighbourC = mesh.faceNeighbour();

    label band = 0;

    forAll(oldFneighbourC, faceI)
    {
        label diff = oldFneighbourC[faceI] - oldFownerC[faceI];

        if (diff > band)
        {
            band = diff;
        }
    }

    Info << "Mesh size: " << mesh.nCells()
        << "    band before renumbering: " << band << endl;

    // create a new mesh
    fvMesh* newMeshPtr = mesh.renumberedMesh();
    fvMesh& newMesh = *newMeshPtr;

    newMesh.write();

    // check the new bandwidth
    band = 0;
    const unallocLabelList& newFownerC = newMesh.owner();
    const unallocLabelList& newFneighbourC = newMesh.neighbour();

    forAll(newFownerC, faceI)
    {
        label diff = newFneighbourC[faceI] - newFownerC[faceI];

        if (diff > band)
        {
            band = diff;
        }
    }

    Info << "        Band after renumbering: " << band << endl;

    // Get times list
    instantList Times = runTime.times();

    forAll(Times, i)
    {
        runTime.setTime(Times[i], i);

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        // Renumbering volScalarField
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search list of objects for volScalarFields
        IOobjectList scalarFields(objects.lookupClass("volScalarField"));

        for
        (
            IOobjectList::iterator scalarFieldIter = scalarFields.begin();
            scalarFieldIter != scalarFields.end();
            ++scalarFieldIter
        )
        {
            Info<< "    renumbering " << scalarFieldIter()->name()
                << endl;

            // Read field
            volScalarField theta
            (
                *scalarFieldIter(),
                mesh
            );

            // Renumber field
            volScalarField thetaCpy
            (
                IOobject
                (
                    (const word&)(scalarFieldIter()->name()),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh.renumber(newMesh, theta)
            );

            // Write field
            thetaCpy.write();
        }

        // Renumbering volVectorField
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search list of objects for volVectorFields
        IOobjectList vectorFields(objects.lookupClass("volVectorField"));

        for
        (
            IOobjectList::iterator vectorFieldIter = vectorFields.begin();
            vectorFieldIter != vectorFields.end();
            ++vectorFieldIter
        )
        {
            Info<< "    renumbering " << vectorFieldIter()->name()
                << endl;

            // Read field
            volVectorField theta
            (
                *vectorFieldIter(),
                mesh
            );

            // Renumber field
            volVectorField thetaCpy
            (
                IOobject
                (
                    (const word&)(vectorFieldIter()->name()),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh.renumber(newMesh, theta)
            );

            // Write field
            thetaCpy.write();
        }

        // Renumbering volTensorField
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search list of objects for volTensorFields
        IOobjectList tensorFields(objects.lookupClass("volTensorField"));

        for
        (
            IOobjectList::iterator tensorFieldIter = tensorFields.begin();
            tensorFieldIter != tensorFields.end();
            ++tensorFieldIter
        )
        {
            Info<< "    renumbering " << tensorFieldIter()->name()
                << endl;

            // Read field
            volTensorField theta
            (
                *tensorFieldIter(),
                mesh
            );

            // Renumber field
            volTensorField thetaCpy
            (
                IOobject
                (
                    (const word&)(tensorFieldIter()->name()),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh.renumber(newMesh, theta)
            );

            // Write field
            thetaCpy.write();
        }
    }

    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
