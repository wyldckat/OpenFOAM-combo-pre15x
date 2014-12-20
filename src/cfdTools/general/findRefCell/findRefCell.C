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

\*---------------------------------------------------------------------------*/

#include "findRefCell.H"
#include "polyMesh.H"
#include "cyclicPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::onConstraintPatch(const polyMesh& mesh, const label refCelli)
{
    bool ocp = false;

    const cell& refCell = mesh.cells()[refCelli];

    forAll(refCell, i)
    {
        label facei = refCell[i];

        if (!mesh.isInternalFace(facei))
        {
            const std::type_info& patchType = typeid
                (mesh.boundaryMesh()[mesh.boundaryMesh().whichPatch(facei)]);

            if
            (
                patchType == typeid(cyclicPolyPatch)
             || patchType == typeid(symmetryPolyPatch)
             || patchType == typeid(processorPolyPatch)
            )
            {
                ocp = true;
            }
        }
    }

    return ocp;
}


Foam::label Foam::findRefCell(const polyMesh& mesh, const label refCelli)
{
    if (Pstream::master())
    {
        label celli = refCelli;

        while (celli < mesh.nCells() && onConstraintPatch(mesh, celli))
        {
            celli++;
        }

        if (celli == mesh.nCells())
        {
            FatalErrorIn("findRefCell(const polyMesh&, const label refCelli)")
                << "Cannot find a cell not on a constraint boundary "
                   "starting from cell " << refCelli
                << exit(FatalError);
        }

        if (celli != refCelli)
        {
            WarningIn
            (
                "Foam::findRefCell(const polyMesh& mesh, const label refCelli)"
            )   << "Requested reference cell " << refCelli
                << " is on a constraint boundary, "
                   "selecting reference cell " << celli
                << endl;
        }

        return celli;
    }
    else
    {
        return refCelli;
    }
}


void Foam::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue
)
{
    if (field.needReference())
    {
        word refCellName = field.name() + "RefCell";
        word refValueName = field.name() + "RefValue";

        refCelli = readLabel(dict.lookup(refCellName));

        label refCellNew = findRefCell(field.mesh(), refCelli);

        if (refCellNew != refCelli)
        {
            refCelli = refCellNew;

            const_cast<dictionary&>(dict).remove(refCellName);
            const_cast<dictionary&>(dict).add(refCellName, refCelli);
        }

        refValue = readScalar(dict.lookup(refValueName));
    }
}


// ************************************************************************* //
