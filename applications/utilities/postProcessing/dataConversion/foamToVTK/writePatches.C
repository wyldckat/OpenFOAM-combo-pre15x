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

#include "writePatches.H"
#include "OFstream.H"
#include "floatScalar.H"
#include "writeFuns.H"
#include "writePatchGeom.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writePatches
(
    const bool binary,
    const bool nearCellValue,
    const vtkMesh& vMesh,
    const labelList& patches,
    const fileName& fileName,
    const PtrList<volScalarField>& volScalarFields,
    const PtrList<volVectorField>& volVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();

    std::ofstream pStream(fileName.c_str());

    pStream
        << "# vtk DataFile Version 2.0" << std::endl
        << "patches" << std::endl;
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

    // Create PrimitivePatch for all selected patches

    labelList faceLabels(mesh.nFaces() - mesh.nInternalFaces());
    label bFaceI = 0;

    forAll(patches, i)
    {
        label patchI = patches[i];

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if (!isType<emptyPolyPatch>(pp))
        {
            forAll(pp, patchFaceI)
            {
                faceLabels[bFaceI++] = pp.start() + patchFaceI;
            }
        }
    }
    faceLabels.setSize(bFaceI);

    indirectPrimitivePatch allPp
    (
        IndirectList<face>
        (
            mesh.faces(),
            faceLabels
        ),
        mesh.points()
    );

    writePatchGeom(binary, allPp.localFaces(), allPp.localPoints(), pStream);


    //-----------------------------------------------------------------
    //
    // Write data
    // 
    //-----------------------------------------------------------------

    pStream
        << "CELL_DATA " << allPp.size() << std::endl
        << "FIELD attributes "
        << 1 + volScalarFields.size() + volVectorFields.size()
        << std::endl;

    // PatchID
    {
        DynamicList<floatScalar> fField(allPp.size());

        pStream
            << "patchID 1 "
            << allPp.size() << " float" << std::endl;

        forAll(patches, i)
        {
            label patchI = patches[i];

            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            if (!isType<emptyPolyPatch>(pp))
            {
                writeFuns::insert(scalarField(pp.size(), patchI), fField);
            }
        }
        writeFuns::write(pStream, binary, fField);
    }

    // VolScalarFields
    forAll(volScalarFields, fieldI)
    {
        const volScalarField& vsf = volScalarFields[fieldI];

        pStream
            << vsf.name() << " 1 "
            << allPp.size() << " float" << std::endl;

        DynamicList<floatScalar> fField(allPp.size());

        forAll(patches, i)
        {
            label patchI = patches[i];

            const fvPatchScalarField& pField = vsf.boundaryField()[patchI];

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
        }
        writeFuns::write(pStream, binary, fField);
    }

    // VolVectorFields
    forAll(volVectorFields, fieldI)
    {
        const volVectorField& vvf = volVectorFields[fieldI];

        pStream
            << vvf.name() << " 3 "
            << allPp.size() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*allPp.size());

        forAll(patches, i)
        {
            label patchI = patches[i];

            const fvPatchVectorField& pField = vvf.boundaryField()[patchI];

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
        }

        writeFuns::write(pStream, binary, fField);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
