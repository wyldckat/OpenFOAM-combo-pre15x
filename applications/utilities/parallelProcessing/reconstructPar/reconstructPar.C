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
    reconstructPar

Description
    Reconstructs a mesh and fields of a case that is decomposed for parallel
    execution of FOAM.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tetPointFields.H"
#include "IOobjectList.H"
#include "polyMeshReconstructor.H"
#include "geometricFvFieldReconstructor.H"
#include "geometricTetPointFieldReconstructor.H"
#include "reconstructLagrangian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Create mesh\n" << endl;
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    IOdictionary decompositionDict
    (
        IOobject
        (
            "decomposeParDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    int nProcs(readInt(decompositionDict.lookup("numberOfSubdomains")));


    // Create the processor databases, meshes, cell addressing
    // and boundary adressing
    PtrList<Time> databases(nProcs);
    PtrList<fvMesh> fvMeshes(nProcs);

    PtrList<labelIOList> pointProcAddressing(nProcs);
    PtrList<labelIOList> faceProcAddressing(nProcs);
    PtrList<labelIOList> cellProcAddressing(nProcs);
    PtrList<labelIOList> boundaryProcAddressing(nProcs);

    forAll (databases, procI)
    {
        databases.hook
        (
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );

        fvMeshes.hook
        (
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    databases[procI].timeName(),
                    databases[procI]
                )
            )
        );

        pointProcAddressing.hook
        (
            new labelIOList
            (
                IOobject
                (
                    "pointProcAddressing",
                    "constant",
                    fvMeshes[procI].meshSubDir,
                    fvMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        faceProcAddressing.hook
        (
            new labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    "constant",
                    fvMeshes[procI].meshSubDir,
                    fvMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        cellProcAddressing.hook
        (
            new labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    "constant",
                    fvMeshes[procI].meshSubDir,
                    fvMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        boundaryProcAddressing.hook
        (
            new labelIOList
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    "constant",
                    fvMeshes[procI].meshSubDir,
                    fvMeshes[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    // Foam version 2.1 changes the addressing of faces in faceProcAddressing
    // The following code checks and modifies the addressing for cases where
    // the decomposition has been done with the foam2.0 and earlier tools, but
    // the reconstruction is attempted with version 2.1 or later
    // 
#   include "checkFaceAddressingComp.H"

    // Get times list from first database
    const instantList Times = databases[0].times();

#   include "checkTimeOptions.H"

    polyMeshReconstructor meshReconstructor
    (
        mesh,
        reinterpret_cast<PtrList<polyMesh>&>(fvMeshes),
        pointProcAddressing,
        faceProcAddressing,
        cellProcAddressing,
        boundaryProcAddressing
    );

    // Loop over all times
    for (label timei = startTime; timei < endTime; timei++)
    {
        // Set time for global database
        runTime.setTime(Times[timei], timei);

        Info << "Time = " << runTime.timeName() << endl << endl;

        // Set time for all databases
        forAll (databases, procI)
        {
            databases[procI].setTime(Times[timei], timei);
        }

        // Reconstruct the points for moving mesh cases and write them out
        meshReconstructor.reconstructPoints();


        // Get list of objects from processor0 database
        IOobjectList objects(fvMeshes[0], databases[0].timeName());

        // If there are any FV fields, reconstruct them
        if
        (
            objects.lookupClass(volScalarField::typeName).size()
         || objects.lookupClass(volVectorField::typeName).size()
         || objects.lookupClass(volTensorField::typeName).size()
         || objects.lookupClass(surfaceScalarField::typeName).size()
        )
        {
            Info << "Reconstructing FV fields" << nl << endl;

            geometricFvFieldReconstructor fvReconstructor
            (
                mesh,
                fvMeshes,
                pointProcAddressing,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            fvReconstructor.reconstructFvVolumeFields<scalar>(objects);
            fvReconstructor.reconstructFvVolumeFields<vector>(objects);
            fvReconstructor.reconstructFvVolumeFields<tensor>(objects);

            fvReconstructor.reconstructFvSurfaceFields<scalar>(objects);
        }
        else
        {
            Info << "No FV fields" << nl << endl;
        }

        // If there are any tetFem fields, reconstruct them
        if
        (
            objects.lookupClass(tetPointScalarField::typeName).size()
         || objects.lookupClass(tetPointVectorField::typeName).size()
        )
        {
            Info << "Reconstructing tet point fields" << nl << endl;

            tetPolyMesh tetMesh(mesh);
            PtrList<tetPolyMesh> tetMeshes(fvMeshes.size());

            forAll (tetMeshes, procI)
            {
                tetMeshes.hook(new tetPolyMesh(fvMeshes[procI]));
            }

            geometricTetPointFieldReconstructor tetPointReconstructor
            (
                tetMesh,
                tetMeshes,
                pointProcAddressing,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            // Reconstruct tet point fields
            tetPointReconstructor.reconstructTetPointFields<scalar>(objects);
            tetPointReconstructor.reconstructTetPointFields<vector>(objects);
        }
        else
        {
            Info << "No tetFem fields" << nl << endl;
        }

        
        label lagrangianProcI = -1;

        forAll (databases, procI)
        {
            if 
            (
                IOobject
                (
                    "positions",
                    databases[procI].timeName(),
                    "lagrangian",
                    databases[procI],
                    IOobject::NO_READ
                ).headerOk()
            )
            {
                lagrangianProcI = procI;

                if 
                (
                    objects.lookupClass(IOField<scalar>::typeName).size()
                 || objects.lookupClass(IOField<vector>::typeName).size()
                 || objects.lookupClass(IOField<tensor>::typeName).size()
                )
                {
                    break;
                }
            }
        }

        if (lagrangianProcI != -1)
        {
            Info << "Reconstructing lagrangian fields" << nl << endl;

            // Get list of objects from processor0 database
            IOobjectList objects
            (
                databases[lagrangianProcI],
                databases[lagrangianProcI].timeName(),
                "lagrangian"
            );

            reconstructLagrangianPositions
            (
                mesh,
                reinterpret_cast<PtrList<polyMesh>&>(fvMeshes),
                faceProcAddressing,
                cellProcAddressing
            );
            reconstructLagrangianFields<scalar>(runTime, databases, objects);
            reconstructLagrangianFields<vector>(runTime, databases, objects);
            reconstructLagrangianFields<tensor>(runTime, databases, objects);
        }
        else
        {
            Info << "No lagrangian fields" << nl << endl;
        }
    }

    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
