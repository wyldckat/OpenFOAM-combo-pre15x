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
    const scalarField& psi = tpsi();
    const scalarField& Diag = diag();

    for (register label cell=0; cell<psi.size(); cell++)
    {
        Apsi[cell] = Diag[cell]*psi[cell];
    }


#   ifndef scalarATmul

    const label* restrict U = lduAddr_.upperAddr().begin();
    const label* restrict L = lduAddr_.lowerAddr().begin();

    const scalar* restrict Lower = lower().begin();
    const scalar* restrict Upper = upper().begin();
    const scalar* restrict Psi = psi.begin();
    scalar* restrict APsi = Apsi.begin();

    register label nFaces = upper().size();

    for (register label face=0; face<nFaces; face++)
    {
        APsi[U[face]] += Lower[face]*Psi[L[face]];
        APsi[L[face]] += Upper[face]*Psi[U[face]];
    }

#   else

    const scalarField& Lower = lower();
    const scalarField& Upper = upper();
    scalarField& Apsi = tApsi();

    const unallocLabelList& u = lduAddr_.upperAddr();
    const unallocLabelList& l = lduAddr_.lowerAddr();

    register label nFaces = upper().size();

    for (register label face=0; face<nFaces; face++)
    {
        Apsi[u[face]] += Lower[face]*psi[l[face]];
        Apsi[l[face]] += Upper[face]*psi[u[face]];
    }

#   endif

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
    const scalarField& psi = tpsi();
    const scalarField& Diag = diag();

    for (register label cell=0; cell<psi.size(); cell++)
    {
        Tpsi[cell] = Diag[cell]*psi[cell];
    }


#   ifndef scalarATmul

    const label* restrict U = lduAddr_.upperAddr().begin();
    const label* restrict L = lduAddr_.lowerAddr().begin();

    const scalar* restrict Lower = lower().begin();
    const scalar* restrict Upper = upper().begin();
    const scalar* restrict Psi = psi.begin();
    scalar* restrict TPsi = Tpsi.begin();

    register label nFaces = upper().size();

    for (register label face=0; face<nFaces; face++)
    {
        TPsi[U[face]] += Upper[face]*Psi[L[face]];
        TPsi[L[face]] += Lower[face]*Psi[U[face]];
    }

#   else

    const scalarField& Lower = lower();
    const scalarField& Upper = upper();
    scalarField& Tpsi = tTpsi();

    const unallocLabelList& u = lduAddr_.upperAddr();
    const unallocLabelList& l = lduAddr_.lowerAddr();

    register label nFaces = upper().size();

    for (register label face=0; face<l_.size(); face++)
    {
        Tpsi[u[face]] += Upper[face]*psi[l[face]];
        Tpsi[l[face]] += Lower[face]*psi[u[face]];
    }

#   endif

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
