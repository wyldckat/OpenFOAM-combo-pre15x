/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "writeInternal.H"
#include "OFstream.H"
#include "floatScalar.H"
#include "writeFuns.H"
#include "vtkTopo.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writeInternal
(
    const bool binary,
    const vtkMesh& vMesh,
    const autoPtr<pointMesh>& pMeshPtr,
    const fileName& fileName,
    const PtrList<volScalarField>& volScalarFields,
    const PtrList<volVectorField>& volVectorFields,
    const PtrList<pointScalarField>& pointScalarFields,
    const PtrList<pointVectorField>& pointVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();
    const vtkTopo& topo = vMesh.topo();

    std::ofstream vtkStream(fileName.c_str());

    //
    // Write header
    //
    vtkStream
        << "# vtk DataFile Version 2.0" << std::endl
        << mesh.time().caseName() << std::endl;
    if (binary)
    {
        vtkStream << "BINARY" << std::endl;
    }
    else
    {
        vtkStream << "ASCII" << std::endl;
    }
    vtkStream << "DATASET UNSTRUCTURED_GRID" << std::endl;


    //------------------------------------------------------------------
    //
    // Write topology
    // 
    //------------------------------------------------------------------

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



    //------------------------------------------------------------------
    //
    // Write data
    // 
    //------------------------------------------------------------------

    // Construct cell-mapping (= identity if no subset)
    const labelList& cMap = vMesh.cMap();


    // Cell data

    const labelList& superCells = topo.superCells();

    vtkStream
        << "CELL_DATA " << vtkCellTypes.size() << std::endl
        << "FIELD attributes "
        << 1 + volScalarFields.size() + volVectorFields.size()
        << std::endl;


    // Cell ids first
    vtkStream << "cellID 1 " << vtkCellTypes.size() << " int"
        << std::endl;

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


    if (pMeshPtr.valid())
    {
        // Construct interpolation on the raw mesh
        volPointInterpolation pInterp(vMesh.baseMesh(), pMeshPtr());


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
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
