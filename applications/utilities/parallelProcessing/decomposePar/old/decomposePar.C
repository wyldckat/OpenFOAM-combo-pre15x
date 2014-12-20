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
    decomposePar

Description
    Automatically decomposes a mesh and fields of a case for parallel execution
    of FOAM.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorFvPatchFields.H"
#include "domainDecomposition.H"
#include "labelIOList.H"

#include "geometricFvFieldDecomposer.H"

#include "tetPointFields.H"
#include "geometricTetPointFieldDecomposer.H"

#include "tensorIOField.H"
#include "lagrangianFieldDecomposer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("fields", "");
    argList::validOptions.insert("cellDist", "");

#   include "setRootCase.H"

    bool decomposeFieldsOnly(args.options().found("fields"));

    bool writecellDist(args.options().found("cellDist"));

#   include "createTime.H"

    Info<< "Time = " << runTime.timeName() << endl;

    Info<< "Create mesh\n" << endl;
    domainDecomposition mesh
    (
        IOobject
        (
            domainDecomposition::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    // Decompose the mesh
    if (!decomposeFieldsOnly)
    {
        if (dir(runTime.path()/"processor1"))
        {
            FatalErrorIn(args.executable())
                << "Case is already decomposed." << endl
                << "    Please remove processor directories before "
                   "decomposing e.g. using:" << nl
                << "    rm -rf " << runTime.path().c_str() << "/processor*"
                << exit(FatalError);
        }

        mesh.decomposeMesh();

        mesh.writeDecomposition();

        if (writecellDist)
        {
            // Write the decomposition as labelList for use with 'manual'
            // decomposition method.
            OFstream str
            (
                runTime.path()
              / mesh.facesInstance()
              / polyMesh::defaultRegion
              / "cellDecomposition"
            );

            str << mesh.cellToProc();

            Info<< nl << "Written decomposition to "
                <<  mesh.facesInstance()
                  / polyMesh::defaultRegion
                  / "cellDecomposition"
                << nl << endl;

            // Write as volScalarField for postprocessing.
            volScalarField cellDist
            (
                IOobject
                (
                    "cellDist",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("cellDist", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            );

            const labelList& procIds = mesh.cellToProc();
            forAll(procIds, celli)
            {
               cellDist[celli] = procIds[celli]; 
            }

            cellDist.write();
        }
    }

    Cloud<indexedParticle> lagrangianPositions(mesh);

    if (lagrangianPositions.size())
    {
        Info<< nl << "Creating lagrangian particle to processor addressing"
            << endl;
    }

    List<SLList<indexedParticle*>*> cellParticles
    (
        mesh.nCells(),
        reinterpret_cast<SLList<indexedParticle*>*>(NULL)
    );

    label i = 0;

    for
    (
        Cloud<indexedParticle>::iterator iter = lagrangianPositions.begin();
        iter != lagrangianPositions.end();
        ++iter
    )
    {
        iter().index() = i++;

        label celli = iter().cell();

        if (!cellParticles[celli])
        {
            cellParticles[celli] = new SLList<indexedParticle*>();
        }

        cellParticles[celli]->append(&iter());
    }


    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    // Construct the vol fields

    PtrList<volScalarField> volScalarFields;
    geometricFvFieldDecomposer::readFields(mesh, objects, volScalarFields);

    PtrList<volVectorField> volVectorFields;
    geometricFvFieldDecomposer::readFields(mesh, objects, volVectorFields);

    PtrList<volTensorField> volTensorFields;
    geometricFvFieldDecomposer::readFields(mesh, objects, volTensorFields);

    PtrList<surfaceScalarField> surfaceScalarFields;
    geometricFvFieldDecomposer::readFields(mesh, objects, surfaceScalarFields);


    // Search list of objects for tetPointTypeFields
    IOobjectList tetPScalFields(objects.lookupClass("tetPointScalarField"));
    IOobjectList tetPVecFields(objects.lookupClass("tetPointVectorField"));
    IOobjectList tetPTensFields(objects.lookupClass("tetPointTensorField"));

    PtrList<tetPointScalarField> tetPointScalarFields(tetPScalFields.size());
    PtrList<tetPointVectorField> tetPointVectorFields(tetPVecFields.size());
    PtrList<tetPointTensorField> tetPointTensorFields(tetPTensFields.size());

    tetPolyMesh* tetMeshPtr = NULL;

    if
    (
        tetPScalFields.size()
     || tetPVecFields.size()
     || tetPTensFields.size()
    )
    {
        tetMeshPtr = new tetPolyMesh(mesh);
        tetPolyMesh& tetMesh = *tetMeshPtr;

        // Construct the tet point scalar fields
        for
        (
            IOobjectList::iterator tetPScalFieldIter = tetPScalFields.begin();
            tetPScalFieldIter != tetPScalFields.end();
            ++tetPScalFieldIter
        )
        {
            tetPointScalarFields.hook
            (
                new tetPointScalarField
                (
                    *tetPScalFieldIter(),
                    tetMesh
                )
            );
        }

        // Construct the tet point vector fields
        for
        (
            IOobjectList::iterator tetPVecFieldIter = tetPVecFields.begin();
            tetPVecFieldIter != tetPVecFields.end();
            ++tetPVecFieldIter
        )
        {
            tetPointVectorFields.hook
            (
                new tetPointVectorField
                (
                    *tetPVecFieldIter(),
                    tetMesh
                )
            );
        }

        // Construct the tet point tensor fields
        for
        (
            IOobjectList::iterator tetPTensFieldIter = tetPTensFields.begin();
            tetPTensFieldIter != tetPTensFields.end();
            ++tetPTensFieldIter
        )
        {
            tetPointTensorFields.hook
            (
                new tetPointTensorField
                (
                    *tetPTensFieldIter(),
                    tetMesh
                )
            );
        }
    }


    IOobjectList lagrangianObjects(mesh, runTime.timeName(), "lagrangian");

    PtrList<scalarIOField> lagrangianScalarFields;
    lagrangianFieldDecomposer::readFields
    (
        lagrangianObjects,
        lagrangianScalarFields
    );

    PtrList<vectorIOField> lagrangianVectorFields;
    lagrangianFieldDecomposer::readFields
    (
        lagrangianObjects,
        lagrangianVectorFields
    );

    PtrList<tensorIOField> lagrangianTensorFields;
    lagrangianFieldDecomposer::readFields
    (
        lagrangianObjects,
        lagrangianTensorFields
    );


    Info<< endl;


    // split the fields over processors
    for (label procI = 0; procI < mesh.nProcs(); procI++)
    {
        Info<< "Processor " << procI << ": field transfer" << endl;

        // open the database
        Time processorDb
        (
            Time::controlDictName,
            args.rootPath(),
            args.caseName()/fileName(word("processor") + name(procI))
        );

        processorDb.setTime(runTime);

        // read the mesh
        fvMesh procMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                processorDb.timeName(),
                processorDb
            )
        );

        // read the addressing information
        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // read the addressing information
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                "constant",
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // FV fields
        if
        (
            volScalarFields.size()
         || volVectorFields.size()
         || volTensorFields.size()
         || surfaceScalarFields.size()
        )
        {
            geometricFvFieldDecomposer fieldDecomposer
            (
                procMesh,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            fieldDecomposer.decomposeFields(volScalarFields);
            fieldDecomposer.decomposeFields(volVectorFields);
            fieldDecomposer.decomposeFields(volTensorFields);
            fieldDecomposer.decomposeFields(surfaceScalarFields);
        }


        // tetPoint fields
        if
        (
            tetPointScalarFields.size()
         || tetPointVectorFields.size()
         || tetPointTensorFields.size()
        )
        {
            tetPolyMesh& tetMesh = *tetMeshPtr;

            tetPolyMesh procTetMesh(procMesh);

            // read the point addressing information
            labelIOList pointProcAddressing
            (
                IOobject
                (
                    "pointProcAddressing",
                    "constant",
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            geometricTetPointFieldDecomposer fieldDecomposer
            (
                tetMesh,
                procTetMesh,
                pointProcAddressing,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            fieldDecomposer.decomposeFields(tetPointScalarFields);
            fieldDecomposer.decomposeFields(tetPointVectorFields);
            fieldDecomposer.decomposeFields(tetPointTensorFields);
        }


        // If there is lagrangian data write it out
        if (lagrangianPositions.size())
        {
            lagrangianFieldDecomposer fieldDecomposer
            (
                mesh,
                procMesh,
                cellProcAddressing,
                lagrangianPositions,
                cellParticles
            );

            // Lagrangian fields
            if
            (
                lagrangianScalarFields.size()
             || lagrangianVectorFields.size()
             || lagrangianTensorFields.size()
            )
            {
                fieldDecomposer.decomposeFields(lagrangianScalarFields);
                fieldDecomposer.decomposeFields(lagrangianVectorFields);
                fieldDecomposer.decomposeFields(lagrangianTensorFields);
            }
        }
    }


    // Clean up the tetMesh pointer
    if (tetMeshPtr)
    {
        delete tetMeshPtr;
    }

    Info<< "\nEnd.\n" << endl;

    return(0);
}


// ************************************************************************* //
