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

#include "gamgSolver.H"
#include "vector2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::gamgSolver::scalingFactor
(
    scalarField& field,
    const scalarField& source,
    const scalarField& A,
    const scalarField& D
) const
{
    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    forAll(field, i)
    {
        scalingFactorNum += source[i]*field[i];
        scalingFactorDenom += A[i]*field[i];
        field[i] += (source[i] - A[i])/D[i];
    }

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    reduce(scalingVector, sumOp<vector2D>());
    return scalingVector.x()/stabilise(scalingVector.y(), VSMALL);
}


// ************************************************************************* //
