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
    legacy VTK file format writer.
    - handles volScalar, volVector, pointScalar, pointVector fields.
    - mesh topo changes.
    - both ascii and binary.
    - single time step writing.
    - write subset only.
    - automatic decomposition of cells; polygons on boundary undecomposed since
      handled by vtk.

    Note: mesh subset is handled by vtkMesh. Slight inconsistency in
    interpolation: on the internal field it interpolates the whole volfield
    to the whole-mesh pointField and then selects only those values it needs
    for the subMesh (using the meshSubset cellMap(), pointMap() functions).
    For the patches however it uses the meshSubset.interpolate function
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
#include "meshSubset.H"
#include "interpolateToCell.H"
#include "faceZoneMesh.H"

#include "vtkTopo.H"
#include "vtkMesh.H"
#include "readFields.H"
#include "writeFuns.H"
#include "writeLagrangian.H"
#include "writePatch.H"
#include "writeFaceSet.H"
#include "writePatchGeom.H"
#include "writeSurfFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


static const label VTK_TETRA      = 10;
static const label VTK_PYRAMID    = 14;
static const label VTK_WEDGE      = 13;
static const label VTK_HEXAHEDRON = 12;


template<class GeoField>
void print(Ostream& os, const ptrList<GeoField>& flds)
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
    argList::noParallel();
#   include "addTimeOptions.H"
    argList::validOptions.insert("mesh", "mesh name");
    argList::validOptions.insert("fields", "fields");
    argList::validOptions.insert("cellSet", "cellSet name");
    argList::validOptions.insert("faceSet", "faceSet name");
    argList::validOptions.insert("ascii","");
    argList::validOptions.insert("surfaceFields","");

#   include "setRootCase.H"
#   include "createTime.H"


    bool binary = !args.options().found("ascii");

    if (binary && (sizeof(floatScalar) != 4 || sizeof(label) != 4))
    {
        FatalErrorIn(args.executable())
            << "floatScalar and/or label are not 4 bytes in size" << nl
            << "Hence cannot use binary VTK format. Please use -ascii"
            << exit(FatalError);
    }

    word cellSetName;
    word vtkName;

    if (args.options().found("cellSet"))
    {
        cellSetName = args.options()["cellSet"];
        vtkName = cellSetName;
    }
    else
    {
        vtkName = runTime.caseName();
    }

    word meshName = meshSubset::defaultRegion;

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
            Info<< "Keeping old VTK files in " << fvPath << nl << endl;
        }
        else
        {
            Info<< "Deleting old VTK files in " << fvPath << nl << endl;

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

        Info<< "Time " << Times[i].name() << endl;

        // Check for new polyMesh/ and update mesh, meshSubset and cell
        // decomposition.
        polyMesh::readUpdateState meshState = vMesh.readUpdate();

        const fvMesh& mesh = vMesh.mesh();
        const vtkTopo& topo = vMesh.topo();

        if
        (
            meshState == polyMesh::TOPO_CHANGE
         || meshState == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "    Original cells:" << mesh.nCells()
                << " points:" << mesh.nPoints()
                << "   Additional cells:" << topo.superCells().size()
                << "  additional points:" << topo.addPointCellLabels().size()
                << nl << endl;
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

            Info<< "    FaceSet   : " << patchFileName << endl;

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

        ptrList<volScalarField> volScalarFields;
        readFields(vMesh, objects, selectedFields, volScalarFields);
        Info<< "    volScalarFields   :";
        print(Info, volScalarFields);

        ptrList<volVectorField> volVectorFields;
        readFields(vMesh, objects, selectedFields, volVectorFields);
        Info<< "    volVectorFields   :";
        print(Info, volVectorFields);

        ptrList<surfaceScalarField> surfScalarFields;
        ptrList<surfaceVectorField> surfVectorFields;
        if (args.options().found("surfaceFields"))
        {
            readFields
            (
                vMesh,
                objects,
                selectedFields,
                surfScalarFields
            );
            Info<< "    surfScalarFields  :";
            print(Info, surfScalarFields);

            readFields
            (
                vMesh,
                objects,
                selectedFields,
                surfVectorFields
            );
            Info<< "    surfVectorFields  :";
            print(Info, surfVectorFields);
        }

        pointMesh pMesh(vMesh);

        ptrList<pointScalarField> pointScalarFields;
        ptrList<pointVectorField> pointVectorFields;
        if (!vMesh.useSubMesh())
        {
            readFields(pMesh, objects, selectedFields, pointScalarFields);
            Info<< "    pointScalarFields :";
            print(Info, pointScalarFields);

            readFields(pMesh, objects, selectedFields, pointVectorFields);
            Info<< "    pointVectorFields :";
            print(Info, pointVectorFields);
        }
        else
        {
            Info<< "    pointScalarFields : Not supported on subMesh\n";
            Info<< "    pointVectorFields : Not supported on subMesh\n";
        }

        Info<< endl;


        //
        // Create file and write header
        //
        fileName vtkFileName
        (
            fvPath/vtkName + "_" + Foam::name(runTime.timeIndex()) + ".vtk"
        );

        Info<< "    Internal  : " << vtkFileName << endl;

        std::ofstream vtkStream(vtkFileName.c_str());


        //
        // Write header
        //
        vtkStream
            << "# vtk DataFile Version 2.0" << std::endl
            << runTime.caseName() << std::endl;
        if (binary)
        {
            vtkStream << "BINARY" << std::endl;
        }
        else
        {
            vtkStream << "ASCII" << std::endl;
        }
        vtkStream << "DATASET UNSTRUCTURED_GRID" << std::endl;


        //---------------------------------------------------------------------
        //
        // Write topology
        // 
        //---------------------------------------------------------------------

        const labelList& addPointCellLabels = topo.addPointCellLabels();
        const label nTotPoints = mesh.nPoints() + addPointCellLabels.size();

        vtkStream
            << "POINTS " << nTotPoints
            << " float" << std::endl;

        DynamicList<floatScalar> ptField(3*nTotPoints);

        writeFuns::insert(mesh.points(), ptField);

        const pointField& ctrs = mesh.cellCentres();
        forAll(addPointCellLabels, api)
        {
            writeFuns::insert(ctrs[addPointCellLabels[api]], ptField);
        }
        writeFuns::write(vtkStream, binary, ptField);


        //
        // Write cells
        //

        const labelListList& vtkVertLabels = topo.vertLabels();

        // Count total number of vertices referenced.
        label nFaceVerts = 0;

        forAll(vtkVertLabels, cellI)
        {
            nFaceVerts += vtkVertLabels[cellI].size() + 1;
        }

        vtkStream << "CELLS " << vtkVertLabels.size() << ' ' << nFaceVerts
            << std::endl;


        DynamicList<label> vertLabels(nFaceVerts);

        forAll(vtkVertLabels, cellI)
        {
            const labelList& vtkVerts = vtkVertLabels[cellI];

            vertLabels.append(vtkVerts.size());

            writeFuns::insert(vtkVerts, vertLabels);
        }
        writeFuns::write(vtkStream, binary, vertLabels);



        const labelList& vtkCellTypes = topo.cellTypes();

        vtkStream << "CELL_TYPES " << vtkCellTypes.size() << std::endl;

        // Make copy since writing might swap stuff.
        DynamicList<label> cellTypes(vtkCellTypes.size());

        writeFuns::insert(vtkCellTypes, cellTypes);

        writeFuns::write(vtkStream, binary, cellTypes);



        //---------------------------------------------------------------------
        //
        // Write data
        // 
        //---------------------------------------------------------------------

        // Construct cell-mapping (= identity if no subset)
        const labelList& cMap = vMesh.cMap();


        // Cell data

        const labelList& superCells = topo.superCells();

        vtkStream
            << "CELL_DATA " << vtkCellTypes.size() << std::endl
            << "FIELD attributes "
            << 1 + volScalarFields.size() + volVectorFields.size() << std::endl;


        // Cell ids first
        vtkStream << "cellID 1 " << vtkCellTypes.size() << " int" << std::endl;

        labelList cellId(vtkCellTypes.size());
        label labelI = 0;

        forAll(mesh.cells(), cellI)
        {
            cellId[labelI++] = cMap[cellI];
        }
        forAll(superCells, superCellI)
        {
            label origCellI = cMap[superCells[superCellI]];

            cellId[labelI++] = origCellI;
        }
        writeFuns::write(vtkStream, binary, cellId);


        // VolScalarFields
        forAll(volScalarFields, fieldI)
        {
            const volScalarField& vsf = volScalarFields[fieldI];

            writeFuns::writeVSF(vtkStream, binary, vsf, vMesh);
        }

        // VolVectorFields
        forAll(volVectorFields, fieldI)
        {
            const volVectorField& vvf = volVectorFields[fieldI];

            writeFuns::writeVVF(vtkStream, binary, vvf, vMesh);
        }


        // Point data

        // Construct interpolation on the raw mesh
        volPointInterpolation pInterp(vMesh, pMesh);


        vtkStream
            << "POINT_DATA " << nTotPoints
            << std::endl
            << "FIELD attributes "
            << volScalarFields.size()
             + volVectorFields.size()
             + pointScalarFields.size()
             + pointVectorFields.size()
            << std::endl;

        // PointScalarFields
        forAll(pointScalarFields, fieldI)
        {
            const pointScalarField& psf = pointScalarFields[fieldI];

            writeFuns::writePSF(vtkStream, binary, psf, vMesh);
        }

        // PointVectorFields
        forAll(pointVectorFields, fieldI)
        {
            const pointVectorField& pvf = pointVectorFields[fieldI];

            writeFuns::writePVF(vtkStream, binary, pvf, vMesh);
        }

        // VolScalarFields
        forAll(volScalarFields, fieldI)
        {
            const volScalarField& vsf = volScalarFields[fieldI];

            // Interpolate to points
            pointScalarField psf(pInterp.interpolate(vsf));

            writeFuns::writeIntVSF(vtkStream, binary, vsf, psf, vMesh);
        }

        // VolVectorFields
        forAll(volVectorFields, fieldI)
        {
            const volVectorField& vvf = volVectorFields[fieldI];

            // Interpolate to points
            pointVectorField pvf(pInterp.interpolate(vvf));

            writeFuns::writeIntVVF(vtkStream, binary, vvf, pvf, vMesh);
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
               /"internalFaces"
               + "_"
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

            Info<< "    Patch     : " << patchFileName << endl;

            writePatch
            (
                binary,
                vMesh,
                patchI,
                patchFileName,
                volScalarFields,
                volVectorFields
            );
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

            Info<< "    FaceZone  : " << patchFileName << endl;

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

        wordList scalarNames(sprayObjects.names(scalarIOField::typeName));
        wordList vectorNames(sprayObjects.names(vectorIOField::typeName));

        if
        (
            !vMesh.useSubMesh()
         && scalarNames.size() > 0
         && vectorNames.size() > 0
        )
        {
            mkDir(fvPath / "lagrangian");

            fileName lagrFileName
            (
                fvPath/"lagrangian"/vtkName
              + "_" + name(runTime.timeIndex()) + ".vtk"
            );

            Info<< "    Lagrangian: " << lagrFileName << endl;
            
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

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
