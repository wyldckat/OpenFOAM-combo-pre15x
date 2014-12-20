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
    lduMatrix member operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void lduMatrix::sumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();

    const unallocLabelList& l = lduAddr_.lowerAddr();
    const unallocLabelList& u = lduAddr_.upperAddr();

    for (register label face=0; face<l.size(); face++)
    {
        Diag[l[face]] += Lower[face];
        Diag[u[face]] += Upper[face];
    }
}


void lduMatrix::negSumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();

    const unallocLabelList& l = lduAddr_.lowerAddr();
    const unallocLabelList& u = lduAddr_.upperAddr();

    for (register label face=0; face<l.size(); face++)
    {
        Diag[l[face]] -= Lower[face];
        Diag[u[face]] -= Upper[face];
    }
}


void lduMatrix::relax
(
    const FieldField<Field, scalar>& intCoeffsCmptAvg,
    const FieldField<Field, scalar>& magCoupleBouCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const scalar alpha
)
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

    scalarField& Diag = diag();
    scalarField sumOff(Diag.size(), 0.0);

    const unallocLabelList& l = lduAddr_.lowerAddr();
    const unallocLabelList& u = lduAddr_.upperAddr();

    for (register label face = 0; face < l.size(); face++)
    {
        sumOff[u[face]] += mag(Lower[face]);
        sumOff[l[face]] += mag(Upper[face]);
    }

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the summ-off-diagonal
    forAll (interfaces, patchI)
    {
        if (interfaces[patchI]->coupled())
        {
            const unallocLabelList& pa = lduAddr_.patchAddr(patchI);
            const scalarField& iCoeffs = intCoeffsCmptAvg[patchI];
            const scalarField& pCoeffs = magCoupleBouCoeffs[patchI];

            for (register label face = 0; face < pa.size(); face++)
            {
                Diag[pa[face]] += iCoeffs[face];
                sumOff[pa[face]] += pCoeffs[face];
            }
        }
    }

    Diag = max(Diag, sumOff);
    Diag /= alpha;

    // Remove the interface internal coefficients from the diagonal
    forAll (interfaces, patchI)
    {
        if (interfaces[patchI]->coupled())
        {
 
            const unallocLabelList& pa = lduAddr_.patchAddr(patchI);
            const scalarField& iCoeffs = intCoeffsCmptAvg[patchI];

            for (register label face = 0; face < pa.size(); face++)
            {
                Diag[pa[face]] -= iCoeffs[face];
            }
        }
    }
}


tmp<scalarField> lduMatrix::residual
(
    const scalarField& psi,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt
) const
{
    scalarField Apsi(psi.size());
    Amul(Apsi, psi, coupleBouCoeffs, interfaces, cmpt);
    return -Apsi;
}


tmp<scalarField> lduMatrix::residual
(
    const scalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt
) const
{
    return source + residual(psi, coupleBouCoeffs, interfaces, cmpt);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// member operators

void lduMatrix::operator=(const lduMatrix& A)
{
    if (this == &A)
    {
        FatalError
            << "lduMatrix::operator=(const lduMatrix&) : "
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (A.lowerPtr_)
    {
        lower() = A.lower();
    }
    else if (lowerPtr_)
    {
        delete lowerPtr_;
        lowerPtr_ = NULL;
    }

    if (A.upperPtr_)
    {
        upper() = A.upper();
    }
    else if (upperPtr_)
    {
        delete upperPtr_;
        upperPtr_ = NULL;
    }

    if (A.diagPtr_)
    {
        diag() = A.diag();
    }
}


void lduMatrix::negate()
{
    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (diagPtr_)
    {
        diagPtr_->negate();
    }
}


void lduMatrix::operator+=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() += A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        if (upperPtr_)
        {
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
        }
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() += A.upper();
        lower() += A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.upperPtr_)
        {
            lower() += A.upper();
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
            upper() += A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() += A.lower();
        upper() += A.upper();
    }
    else if (diagonal())
    {
        if (A.upperPtr_)
        {
            upper() = A.upper();
        }

        if (A.lowerPtr_)
        {
            lower() = A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        FatalErrorIn("lduMatrix::operator+=(const lduMatrix& A)")
            << "Unknown matrix type combination"
            << abort(FatalError);
    }
}


void lduMatrix::operator-=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() -= A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        if (upperPtr_)
        {
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
        }
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() -= A.upper();
        lower() -= A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.upperPtr_)
        {
            lower() -= A.upper();
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
            upper() -= A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() -= A.lower();
        upper() -= A.upper();
    }
    else if (diagonal())
    {
        if (A.upperPtr_)
        {
            upper() = -A.upper();
        }

        if (A.lowerPtr_)
        {
            lower() = -A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        FatalErrorIn("lduMatrix::operator-=(const lduMatrix& A)")
            << "Unknown matrix type combination"
            << abort(FatalError);
    }
}


void lduMatrix::operator*=(const scalarField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (upperPtr_)
    {
        scalarField& upper = *upperPtr_;

        const unallocLabelList& l = lduAddr_.lowerAddr();

        for (register label face=0; face<upper.size(); face++)
        {
            upper[face] *= sf[l[face]];
        }
    }

    if (lowerPtr_)
    {
        scalarField& lower = *lowerPtr_;

        const unallocLabelList& u = lduAddr_.upperAddr();

        for (register label face=0; face<lower.size(); face++)
        {
            lower[face] *= sf[u[face]];
        }
    }
}


void lduMatrix::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
