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
    legacy VTK file format writer.
    - handles volScalar, volVector, pointScalar, pointVector, surfaceScalar
      fields.
    - mesh topo changes.
    - both ascii and binary.
    - single time step writing.
    - write subset only.
    - automatic decomposition of cells; polygons on boundary undecomposed since
      handled by vtk.


    Current options:
    -mesh           : different mesh name (instead of -region)
    -fields         : selected fields only
    -surfaceFields  : write surfaceScalarFields (e.g. phi)
    -cellSet, -faceSet : set display
    -ascii
    -nearCellValue  : output cell value on patches instead of patch value itself
    -noInternal     : do not generate file for mesh, only for patches
    -noPointValues  : no pointFields
    -allPatches     : combine all patches into one big file

    Note: mesh subset is handled by vtkMesh. Slight inconsistency in
    interpolation: on the internal field it interpolates the whole volfield
    to the whole-mesh pointField and then selects only those values it needs
    for the subMesh (using the fvMeshSubset cellMap(), pointMap() functions).
    For the patches however it uses the fvMeshSubset.interpolate function
    to directly interpolate the whole-mesh values onto the subset patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "PrimitivePatchInterpolation.H"
#include "emptyPolyPatch.H"
#include "Cloud.H"
#include "passiveParticle.H"
#include "IOobject.H"
#include "IOField.H"
#include "fvMeshSubset.H"
#include "interpolateToCell.H"
#include "faceZoneMesh.H"

#include "vtkTopo.H"
#include "vtkMesh.H"
#include "readFields.H"
#include "writeFuns.H"
#include "writeInternal.H"
#include "writeLagrangian.H"
#include "writePatch.H"
#include "writeFaceSet.H"
#include "writePatchGeom.H"
#include "writeSurfFields.H"
#include "writePatches.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


static const label VTK_TETRA      = 10;
static const label VTK_PYRAMID    = 14;
static const label VTK_WEDGE      = 13;
static const label VTK_HEXAHEDRON = 12;


template<class GeoField>
void print(Ostream& os, const PtrList<GeoField>& flds)
{
    forAll(flds, i)
    {
        os<< ' ' << flds[i].name();
    }
    os << endl;
}

// Main program:

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
    argList::validOptions.insert("mesh", "mesh name");
    argList::validOptions.insert("fields", "fields");
    argList::validOptions.insert("cellSet", "cellSet name");
    argList::validOptions.insert("faceSet", "faceSet name");
    argList::validOptions.insert("ascii","");
    argList::validOptions.insert("surfaceFields","");
    argList::validOptions.insert("nearCellValue","");
    argList::validOptions.insert("noInternal","");
    argList::validOptions.insert("noPointValues","");
    argList::validOptions.insert("allPatches","");

#   include "setRootCase.H"
#   include "createTime.H"


    bool doWriteInternal = !args.options().found("noInternal");

    bool binary = !args.options().found("ascii");

    if (binary && (sizeof(floatScalar) != 4 || sizeof(label) != 4))
    {
        FatalErrorIn(args.executable())
            << "floatScalar and/or label are not 4 bytes in size" << nl
            << "Hence cannot use binary VTK format. Please use -ascii"
            << exit(FatalError);
    }

    bool nearCellValue = args.options().found("nearCellValue");

    if (nearCellValue)
    {
        WarningIn(args.executable())
            << "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    bool noPointValues = args.options().found("noPointValues");

    if (noPointValues)
    {
        WarningIn(args.executable())
            << "Outputting cell values only" << nl << endl;
    }

    bool allPatches = args.options().found("allPatches");


    word cellSetName;
    string vtkName;

    if (args.options().found("cellSet"))
    {
        cellSetName = args.options()["cellSet"];
        vtkName = cellSetName;
    }
    else if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        vtkName = runTime.caseName();

        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }
    else
    {
        vtkName = runTime.caseName();
    }

    word meshName = fvMeshSubset::defaultRegion;

    // make a directory called VTK in the case
    fileName fvPath(runTime.path()/"VTK");

    if (args.options().found("mesh"))
    {
        meshName = args.options()["mesh"];
        fvPath = fvPath/meshName;
    }

    if (dir(fvPath))
    {
        if
        (
            args.options().found("time")
         || args.options().found("latestTime")
         || cellSetName.size() > 0
        )
        {
            Pout<< "Keeping old VTK files in " << fvPath << nl << endl;
        }
        else
        {
            Pout<< "Deleting old VTK files in " << fvPath << nl << endl;

            rmDir(fvPath);
        }
    }

    mkDir(fvPath);


    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptionsNoConstant.H"

    runTime.setTime(Times[startTime], startTime);

    // Current mesh.
    vtkMesh vMesh
    (
        IOobject
        (
            meshName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        ),
        cellSetName
    );


    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Pout<< "Time " << Times[i].name() << endl;

        // Check for new polyMesh/ and update mesh, fvMeshSubset and cell
        // decomposition.
        polyMesh::readUpdateState meshState = vMesh.readUpdate();

        const fvMesh& mesh = vMesh.mesh();

        if
        (
            meshState == polyMesh::TOPO_CHANGE
         || meshState == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Pout<< "    Read new mesh" << nl << endl;
        }


        // If faceSet: write faceSet only (as polydata)
        if (args.options().found("faceSet"))
        {
            // Load the faceSet
            faceSet set(mesh, args.options()["faceSet"]);

            // Filename as if patch with same name.
            mkDir(fvPath/set.name());

            fileName patchFileName
            (
                fvPath/set.name()/set.name()
              + "_"
              + name(runTime.timeIndex())
              + ".vtk"
            );

            Pout<< "    FaceSet   : " << patchFileName << endl;

            writeFaceSet(binary, vMesh, set, patchFileName);

            continue;
        }


        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        HashSet<word> selectedFields;
        if (args.options().found("fields"))
        {
            IStringStream(args.options()["fields"])() >> selectedFields;
        }

        // Construct the vol fields (on the original mesh if subsetted)

        PtrList<volScalarField> volScalarFields;
        readFields(vMesh, objects, selectedFields, volScalarFields);
        Pout<< "    volScalarFields   :";
        print(Pout, volScalarFields);

        PtrList<volVectorField> volVectorFields;
        readFields(vMesh, objects, selectedFields, volVectorFields);
        Pout<< "    volVectorFields   :";
        print(Pout, volVectorFields);

        PtrList<surfaceScalarField> surfScalarFields;
        PtrList<surfaceVectorField> surfVectorFields;
        if (args.options().found("surfaceFields"))
        {
            readFields
            (
                vMesh,
                objects,
                selectedFields,
                surfScalarFields
            );
            Pout<< "    surfScalarFields  :";
            print(Pout, surfScalarFields);

            readFields
            (
                vMesh,
                objects,
                selectedFields,
                surfVectorFields
            );
            Pout<< "    surfVectorFields  :";
            print(Pout, surfVectorFields);
        }

        // Construct pointMesh only if nessecary since constructs edge
        // addressing (expensive on polyhedral meshes)
        autoPtr<pointMesh> pMeshPtr(NULL);
        if (noPointValues)
        {
            Pout<< "    pointScalarFields : switched off"
                << " (\"-noPointValues\" option)\n";
            Pout<< "    pointVectorFields : switched off"
                << " (\"-noPointValues\" option)\n";
        }
        else
        {
            pMeshPtr.reset(new pointMesh(vMesh));
        }

        PtrList<pointScalarField> pointScalarFields;
        PtrList<pointVectorField> pointVectorFields;
        if (pMeshPtr.valid() && !vMesh.useSubMesh())
        {
            readFields(pMeshPtr(), objects, selectedFields, pointScalarFields);
            Pout<< "    pointScalarFields :";
            print(Pout, pointScalarFields);

            readFields(pMeshPtr(), objects, selectedFields, pointVectorFields);
            Pout<< "    pointVectorFields :";
            print(Pout, pointVectorFields);
        }

        Pout<< endl;

        if (doWriteInternal)
        {
            //
            // Create file and write header
            //
            fileName vtkFileName
            (
                fvPath/vtkName + "_" + Foam::name(runTime.timeIndex()) + ".vtk"
            );

            Pout<< "    Internal  : " << vtkFileName << endl;

            writeInternal
            (
                binary,              // write binary
                vMesh,
                pMeshPtr,
                vtkFileName,
                volScalarFields,
                volVectorFields,
                pointScalarFields,
                pointVectorFields
            );
        }

        //---------------------------------------------------------------------
        //
        // Write surface fields
        // 
        //---------------------------------------------------------------------

        if (surfScalarFields.size() + surfVectorFields.size() > 0)
        {
            mkDir(fvPath / "surfaceFields");

            fileName surfFileName
            (
                fvPath
               /"surfaceFields"
               /"surfaceFields_"
               + name(runTime.timeIndex())
               + ".vtk"
            );

            writeSurfFields
            (
                binary,
                vMesh,
                surfFileName,
                surfScalarFields,
                surfVectorFields
            );
        }


        //---------------------------------------------------------------------
        //
        // Write patches (POLYDATA file, one for each patch)
        // 
        //---------------------------------------------------------------------

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        if (allPatches)
        {
            DynamicList<label> selectedPatches(patches.size());

            Pout<< "Combining patches:" << endl;

            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if
                (
                    isType<emptyPolyPatch>(pp)
                    || (Pstream::parRun() && isType<processorPolyPatch>(pp))
                )
                {
                    Pout<< "    discarding empty/processor patch " << patchI
                        << " " << pp.name() << endl;
                }
                else
                {
                    selectedPatches.append(patchI);
                    Pout<< "    patch " << patchI << " " << pp.name() << endl;
                }
            }
            selectedPatches.shrink();


            mkDir(fvPath/"allPatches");

            fileName patchFileName;

            if (vMesh.useSubMesh())
            {
                patchFileName =
                    fvPath/"allPatches"/cellSetName
                  + "_"
                  + name(runTime.timeIndex())
                  + ".vtk";
            }
            else
            {
                patchFileName =
                    fvPath/"allPatches"/"allPatches"
                  + "_"
                  + name(runTime.timeIndex())
                  + ".vtk";
            }

            Pout<< "    Combined patches     : " << patchFileName << endl;

            writePatches
            (
                binary,
                nearCellValue,
                vMesh,
                selectedPatches,
                patchFileName,
                volScalarFields,
                volVectorFields
            );
        }
        else
        {
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                mkDir(fvPath/pp.name());

                fileName patchFileName;

                if (vMesh.useSubMesh())
                {
                    patchFileName =
                        fvPath/pp.name()/cellSetName
                      + "_"
                      + name(runTime.timeIndex())
                      + ".vtk";
                }
                else
                {
                    patchFileName =
                        fvPath/pp.name()/pp.name()
                      + "_"
                      + name(runTime.timeIndex())
                      + ".vtk";
                }

                Pout<< "    Patch     : " << patchFileName << endl;

                writePatch
                (
                    binary,
                    nearCellValue,
                    vMesh,
                    patchI,
                    patchFileName,
                    volScalarFields,
                    volVectorFields,
                    pointScalarFields,
                    pointVectorFields
                );
            }
        }

        //---------------------------------------------------------------------
        //
        // Write faceZones (POLYDATA file, one for each zone)
        // 
        //---------------------------------------------------------------------

        const faceZoneMesh& zones = mesh.faceZones();

        forAll(zones, zoneI)
        {
            const faceZone& pp = zones[zoneI];

            mkDir(fvPath/pp.name());

            fileName patchFileName;

            if (vMesh.useSubMesh())
            {
                patchFileName =
                    fvPath/pp.name()/cellSetName
                  + "_"
                  + name(runTime.timeIndex())
                  + ".vtk";
            }
            else
            {
                patchFileName =
                    fvPath/pp.name()/pp.name()
                  + "_"
                  + name(runTime.timeIndex())
                  + ".vtk";
            }

            Pout<< "    FaceZone  : " << patchFileName << endl;

            std::ofstream pStream(patchFileName.c_str());

            pStream
                << "# vtk DataFile Version 2.0" << std::endl
                << pp.name() << std::endl;
            if (binary)
            {
                pStream << "BINARY" << std::endl;
            }
            else
            {
                pStream << "ASCII" << std::endl;
            }
            pStream << "DATASET POLYDATA" << std::endl;

            const primitiveFacePatch& fp = pp();

            writePatchGeom(binary, fp.localFaces(), fp.localPoints(), pStream);
        }



        //---------------------------------------------------------------------
        //
        // Write lagrangian data
        // 
        //---------------------------------------------------------------------

        IOobjectList sprayObjects(mesh, runTime.timeName(), "lagrangian");

        if (!vMesh.useSubMesh() && sprayObjects.size())
        {
            mkDir(fvPath/"lagrangian");

            fileName lagrFileName
            (
                fvPath/"lagrangian"/vtkName
              + "_" + name(runTime.timeIndex()) + ".vtk"
            );

            wordList scalarNames(sprayObjects.names(scalarIOField::typeName));
            wordList vectorNames(sprayObjects.names(vectorIOField::typeName));

            Pout<< "    Lagrangian: " << lagrFileName << endl;
            
            writeLagrangian
            (
                binary,
                vMesh,
                lagrFileName,
                scalarNames,
                vectorNames
            );
        }
    }

    Pout << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
