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

\*----------------------------------------------------------------------------*/

#include "syncTools.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Define dummy transform for labels.
Foam::label Foam::transform(const tensor&, const label val)
{
    return val;
}


// Does anyone have couples? Since meshes might have 0 cells and 0 proc
// boundaries need to reduce this info.
bool Foam::syncTools::hasCouples(const polyBoundaryMesh& patches)
{
    bool hasAnyCouples = false;

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            hasAnyCouples = true;
            break;
        }
    }
    return returnReduce(hasAnyCouples, orOp<bool>());
}


void Foam::syncTools::checkTransform
(
    const coupledPolyPatch& pp,
    const bool applySeparation
)
{
    if (!pp.parallel() && pp.forwardT().size() > 1)
    {
        FatalErrorIn("syncTools::checkTransform(const coupledPolyPatch&)")
            << "Non-uniform transformation not supported for point or edge"
            << " fields." << endl
            << "Patch:" << pp.name()
            << abort(FatalError);
    }
    if (applySeparation && pp.separated() && pp.separation().size() > 1)
    {
        FatalErrorIn("syncTools::checkTransform(const coupledPolyPatch&)")
            << "Non-uniform separation vector not supported for point or edge"
            << " fields." << endl
            << "Patch:" << pp.name()
            << abort(FatalError);
    }
}


template <>
void Foam::syncTools::separateList
(
    const vectorField& separation,
    UList<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.

        forAll(field, i)
        {
            field[i] += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] += separation[i];
        }
    }
    else
    {
        FatalErrorIn
        (
            "syncTools::separateList(const vectorField&, UList<T>&)"
        )   << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


// ************************************************************************* //
