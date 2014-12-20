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
    const vtkMesh& vMesh,
    const label patchI,
    const fileName& fileName,
    const ptrList<volScalarField>& volScalarFields,
    const ptrList<volVectorField>& volVectorFields
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

    if (typeid(pp) != typeid(emptyPolyPatch))
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

            const scalarField& pField = vsf.boundaryField()[patchI];

            pStream
                << vsf.name() << " 1 "
                << pField.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(pField.size());

            writeFuns::insert(pField, fField);

            writeFuns::write(pStream, binary, fField);
        }

        // VolVectorFields
        forAll(volVectorFields, fieldI)
        {
            const volVectorField& vvf = volVectorFields[fieldI];

            const vectorField& pField = vvf.boundaryField()[patchI];

            pStream
                << vvf.name() << " 3 "
                << pField.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(3*pField.size());

            writeFuns::insert(pField, fField);

            writeFuns::write(pStream, binary, fField);
        }

        // Vertex data
        pStream
            << "POINT_DATA " << pp.nPoints() << std::endl
            << "FIELD attributes "
            << volScalarFields.size() + volVectorFields.size()
            << std::endl;

        PrimitivePatchInterpolation<primitivePatch> pInter(pp);

        // VolScalarFields
        forAll(volScalarFields, fieldI)
        {
            tmp<volScalarField> tvsf =
                vMesh.interpolate(volScalarFields[fieldI]);

            const volScalarField& vsf = tvsf();

            tmp<scalarField> tpField =
                pInter.faceToPointInterpolate
                (
                    vsf.boundaryField()[patchI]
                );

            const scalarField& pField = tpField();

            pStream
                << vsf.name() << " 1 "
                << pField.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(pField.size());

            writeFuns::insert(pField, fField);

            writeFuns::write(pStream, binary, fField);
        }

        // VolVectorFields
        forAll(volVectorFields, fieldI)
        {
            tmp<volVectorField> tvvf =
                vMesh.interpolate(volVectorFields[fieldI]);
            const volVectorField& vvf = tvvf();

            tmp<vectorField> tpField =
                pInter.faceToPointInterpolate
                (
                    vvf.boundaryField()[patchI]
                );

            const vectorField& pField = tpField();

            pStream
                << vvf.name() << " 3 "
                << pField.size() << " float" << std::endl;

            DynamicList<floatScalar> fField(3*pField.size());

            writeFuns::insert(pField, fField);

            writeFuns::write(pStream, binary, fField);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
