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

\*---------------------------------------------------------------------------*/

#include "writePatch.H"
#include "OFstream.H"
#include "PrimitivePatchInterpolation.H"
#include "floatScalar.H"
#include "writeFuns.H"
#include "writePatchGeom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writePatch
(
    const bool binary,
    const bool nearCellValue,
    const vtkMesh& vMesh,
    const label patchI,
    const fileName& fileName,
    const PtrList<volScalarField>& volScalarFields,
    const PtrList<volVectorField>& volVectorFields,
    const PtrList<pointScalarField>& pointScalarFields,
    const PtrList<pointVectorField>& pointVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();

    const polyPatch& pp = mesh.boundaryMesh()[patchI];

    std::ofstream pStream(fileName.c_str());

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


    //------------------------------------------------------------------
    //
    // Write topology
    // 
    //------------------------------------------------------------------

    writePatchGeom(binary, pp.localFaces(), pp.localPoints(), pStream);


    //-----------------------------------------------------------------
    //
    // Write data
    // 
    //-----------------------------------------------------------------

    if (!isType<emptyPolyPatch>(pp))
    {
        // Face data
        pStream
            << "CELL_DATA " << pp.size() << std::endl
            << "FIELD attributes "
            << volScalarFields.size() + volVectorFields.size()
            << std::endl;

        // VolScalarFields
        forAll(volScalarFields, fieldI)
        {
            const volScalarField& vsf = volScalarFields[fieldI];

            const fvPatchScalarField& pField = vsf.boundaryField()[patchI];

            pStream
                << vsf.name() << " 1 "
                << pField.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(pField.size());

            if (nearCellValue)
            {
                // Write value of neighbouring cell instead of patch value
                writeFuns::insert
                (
                    pField.patch().patchInternalField
                    (
                        vsf.internalField()
                    )(),
                    fField
                );
            }
            else
            {
                writeFuns::insert(pField, fField);
            }

            writeFuns::write(pStream, binary, fField);
        }

        // VolVectorFields
        forAll(volVectorFields, fieldI)
        {
            const volVectorField& vvf = volVectorFields[fieldI];

            const fvPatchVectorField& pField = vvf.boundaryField()[patchI];

            pStream
                << vvf.name() << " 3 "
                << pField.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(3*pField.size());

            if (nearCellValue)
            {
                writeFuns::insert
                (
                    pField.patch().patchInternalField
                    (
                        vvf.internalField()
                    )(),
                    fField
                );
            }
            else
            {
                writeFuns::insert(pField, fField);
            }

            writeFuns::write(pStream, binary, fField);
        }

        // Vertex data
        pStream
            << "POINT_DATA " << pp.nPoints() << std::endl
            << "FIELD attributes "
            << volScalarFields.size()
             + volVectorFields.size()
             + pointScalarFields.size()
             + pointVectorFields.size()
            << std::endl;

        PrimitivePatchInterpolation<primitivePatch> pInter(pp);

        // VolScalarFields
        forAll(volScalarFields, fieldI)
        {
            tmp<volScalarField> tvsf =
                vMesh.interpolate(volScalarFields[fieldI]);

            const volScalarField& vsf = tvsf();

            const fvPatchScalarField& bField = vsf.boundaryField()[patchI];

            DynamicList<floatScalar> fField;

            if (nearCellValue)
            {
                tmp<scalarField> tpField =
                    pInter.faceToPointInterpolate
                    (
                        bField.patch().patchInternalField
                        (
                            vsf.internalField()
                        )()
                    );
                const scalarField& pField = tpField();

                pStream
                    << vsf.name() << " 1 "
                    << pField.size() << " float" << std::endl;

                fField.setSize(pField.size());

                writeFuns::insert(pField, fField);
            }
            else
            {
                tmp<scalarField> tpField = pInter.faceToPointInterpolate(bField);
                const scalarField& pField = tpField();

                pStream
                    << vsf.name() << " 1 "
                    << pField.size() << " float" << std::endl;

                fField.setSize(pField.size());

                writeFuns::insert(pField, fField);
            }
    
            writeFuns::write(pStream, binary, fField);
        }

        // VolVectorFields
        forAll(volVectorFields, fieldI)
        {
            tmp<volVectorField> tvvf =
                vMesh.interpolate(volVectorFields[fieldI]);
            const volVectorField& vvf = tvvf();

            const fvPatchVectorField& bField = vvf.boundaryField()[patchI];

            // Storage for components of vector
            DynamicList<floatScalar> fField;

            if (nearCellValue)
            {
                tmp<vectorField> tpField =
                    pInter.faceToPointInterpolate
                    (
                        bField.patch().patchInternalField
                        (
                            vvf.internalField()
                        )()
                    );
                const vectorField& pField = tpField();

                pStream
                    << vvf.name() << " 3 "
                    << pField.size() << " float" << std::endl;

                fField.setSize(3*pField.size());

                writeFuns::insert(pField, fField);
            }
            else
            {
                tmp<vectorField> tpField = pInter.faceToPointInterpolate(bField);
                const vectorField& pField = tpField();

                pStream
                    << vvf.name() << " 3 "
                    << pField.size() << " float" << std::endl;

                fField.setSize(3*pField.size());

                writeFuns::insert(pField, fField);
            }

            writeFuns::write(pStream, binary, fField);
        }

        // PointScalarFields
        forAll(pointScalarFields, fieldI)
        {
            const pointScalarField& psf = pointScalarFields[fieldI];

            scalarField pf = psf.boundaryField()[patchI].patchInternalField();

            pStream
                << psf.name() << " 1 "
                << pf.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(pf.size());

            writeFuns::insert(pf, fField);

            writeFuns::write(pStream, binary, fField);
        }

        // PointVectorFields
        forAll(pointVectorFields, fieldI)
        {
            const pointVectorField& psf = pointVectorFields[fieldI];

            vectorField pf = psf.boundaryField()[patchI].patchInternalField();

            pStream
                << psf.name() << " 3 "
                << pf.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(3*pf.size());

            writeFuns::insert(pf, fField);

            writeFuns::write(pStream, binary, fField);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
