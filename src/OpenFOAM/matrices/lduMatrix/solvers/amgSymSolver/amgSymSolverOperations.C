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
    Restriction, prolongation and scaling for the amg solver.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "amgSymSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// restrict the fine mesh field to coarse mesh
scalarField amgSymSolver::restrictField
(
    const scalarField& fineField,
    const label fineLevelIndex
) const
{
    // Restriction operator is zero-order (summation)

    // Get addressing
    const labelField& fineToCoarse = restrictAddressing_[fineLevelIndex];

    // Debugging
    if (fineField.size() != fineToCoarse.size())
    {
        FatalErrorIn
        (
            "scalarField amgSymSolver::restrictField(const scalarField& f, "
            "const label fineLevelIndex) const"
        )   << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << fineField.size() << " level = "
            << fineToCoarse.size()
            << abort(FatalError);
    }

    // Dimension the result; the addressing and matrices indices are shifted
    // by one because of separate fine matrix. 
    scalarField result(addrLevels_[fineLevelIndex].size(), 0.0);

    forAll (fineField, cellI)
    {
        result[fineToCoarse[cellI]] += fineField[cellI];
    }

    return result;
}


scalarField amgSymSolver::prolongField
(
    const scalarField& f,
    const label coarseLevelIndex
) const
{
    // Prolongation operator is zero order

    // Get addressing
    const labelField& fineToCoarse = restrictAddressing_[coarseLevelIndex];

    // Dimension the result; the addressing and matrices indices are shifted
    // by one because of separate fine matrix. 
    // Note. Setting to -GREAT is for debigging only; the initial value
    // can be removed
    scalarField result(fineToCoarse.size(), -GREAT);

    forAll (fineToCoarse, cellI)
    {
        result[cellI] = f[fineToCoarse[cellI]];
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
