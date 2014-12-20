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

\*---------------------------------------------------------------------------*/

#include "writeSurfFields.H"
#include "OFstream.H"
#include "floatScalar.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writeSurfFields
(
    const bool binary,
    const vtkMesh& vMesh,
    const fileName& fileName,
    const PtrList<surfaceScalarField>& surfScalarFields,
    const PtrList<surfaceVectorField>& surfVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();

    std::ofstream str(fileName.c_str());

    str << "# vtk DataFile Version 2.0" << std::endl
        << "internalFaces" << std::endl;

    if (binary)
    {
        str << "BINARY" << std::endl;
    }
    else
    {
        str << "ASCII" << std::endl;
    }
    str << "DATASET POLYDATA" << std::endl;

    const pointField& fc = mesh.faceCentres();

    str << "POINTS " << mesh.nInternalFaces() << " float" << std::endl;

    DynamicList<floatScalar> pField(3*mesh.nInternalFaces());

    for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
    {
        writeFuns::insert(fc[faceI], pField);   
    }

    writeFuns::write(str, binary, pField);

    str << "POINT_DATA " << mesh.nInternalFaces() << std::endl
        << "FIELD attributes "
        << surfScalarFields.size() + surfVectorFields.size()
        << std::endl;

    // surfScalarFields
    forAll(surfScalarFields, fieldI)
    {
        const surfaceScalarField& ssf = surfScalarFields[fieldI];

        str << ssf.name() << " 3 "
            << mesh.nInternalFaces() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*mesh.nInternalFaces());

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            vector v(ssf[faceI]*mesh.Sf()[faceI]/mesh.magSf()[faceI]);

            writeFuns::insert(v, fField);
        }

        writeFuns::write(str, binary, fField);
    }
    // surfVectorFields
    forAll(surfVectorFields, fieldI)
    {
        const surfaceVectorField& svf = surfVectorFields[fieldI];

        str << svf.name() << " 3 "
            << mesh.nInternalFaces() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*mesh.nInternalFaces());

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            writeFuns::insert(svf[faceI], fField);
        }

        writeFuns::write(str, binary, fField);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
