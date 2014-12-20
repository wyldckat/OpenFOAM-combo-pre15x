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

Class
    writePatches

Description
    Write set of patches with their data to a single file.
    (does not do interpolated fields)

SourceFiles
    writePatches.C

\*---------------------------------------------------------------------------*/

#ifndef writePatches_H
#define writePatches_H

#include "vtkMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Write patch fields.
// nearCellValue writes the patchInternalField which is useful for e.g.
// walls where the wall velocity is usually not terribly interesting but
// the neighbouring cell value is.
void writePatches
(
    const bool binary,              // write binary
    const bool nearCellValue,       // write near cell value
    const vtkMesh& vMesh,
    const labelList& patches,
    const fileName& fileName,
    const PtrList<volScalarField>& volScalarFields,
    const PtrList<volVectorField>& volVectorFields
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
