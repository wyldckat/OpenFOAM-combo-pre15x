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
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void lduMatrix::Amul
(
    scalarField& Apsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt
) const
{
    scalar* restrict ApsiPtr = Apsi.begin();

    const scalarField& psi = tpsi();
    const scalar* restrict psiPtr = psi.begin();

    const scalar* restrict diagPtr = diag().begin();

    const label* restrict uPtr = lduAddr_.upperAddr().begin();
    const label* restrict lPtr = lduAddr_.lowerAddr().begin();

    const scalar* restrict lowerPtr = lower().begin();
    const scalar* restrict upperPtr = upper().begin();

    register const label nCells = diag().size();
    register const label nFaces = upper().size();

    for (register label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        coupleBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    for (register label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }

    // Update couple interfaces
    updateMatrixInterfaces
    (
        coupleBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    tpsi.clear();
}


void lduMatrix::Tmul
(
    scalarField& Tpsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt
) const
{
    scalar* restrict TpsiPtr = Tpsi.begin();

    const scalarField& psi = tpsi();
    const scalar* restrict psiPtr = psi.begin();

    const scalar* restrict diagPtr = diag().begin();

    const label* restrict uPtr = lduAddr_.upperAddr().begin();
    const label* restrict lPtr = lduAddr_.lowerAddr().begin();

    const scalar* restrict lowerPtr = lower().begin();
    const scalar* restrict upperPtr = upper().begin();

    register const label nCells = diag().size();
    register const label nFaces = upper().size();

    for (register label cell=0; cell<nCells; cell++)
    {
        TpsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        coupleIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    for (register label face=0; face<nFaces; face++)
    {
        TpsiPtr[uPtr[face]] += upperPtr[face]*psiPtr[lPtr[face]];
        TpsiPtr[lPtr[face]] += lowerPtr[face]*psiPtr[uPtr[face]];
    }

    // Update couple interfaces
    updateMatrixInterfaces
    (
        coupleIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    tpsi.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
