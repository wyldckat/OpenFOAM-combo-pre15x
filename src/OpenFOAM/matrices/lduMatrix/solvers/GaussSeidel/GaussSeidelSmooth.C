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
    Gauss-Seidel smoother with fixed number of sweeps.

\*---------------------------------------------------------------------------*/

#include "GaussSeidel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void GaussSeidel::smooth
(
    const word& fieldName,
    scalarField& psi,
    const lduMatrix& m,
    const scalarField& Source,
    const FieldField<Field, scalar>& bouCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt,
    const label nSweeps
)
{
    if (m.diagonal())
    {
        psi = Source/m.diag();
    }
    else if (m.symmetric() || m.asymmetric())
    {
        scalar* restrict psiPtr = psi.begin();

        const scalar* restrict diagPtr = m.diag().begin();

        scalarField bPrime(psi.size());
        scalar* restrict bPrimePtr = bPrime.begin();

        const label* restrict uPtr = m.lduAddr().upperAddr().begin();
        const label* restrict ownStartPtr =
            m.lduAddr().ownerStartAddr().begin();

        const scalar* restrict lowerPtr = m.lower().begin();
        const scalar* restrict upperPtr = m.upper().begin();

        const label nCells = psi.size();

        for (label sweep=0; sweep<nSweeps; sweep++)
        {
            bPrime = Source;

            // Parallel boundary update.  The parallel boundary is treated
            // as an effective jacobi interface in the boundary.
            // Note: there is a change of sign in the coupled
            // interface update.  The reason for this is that the
            // internal coefficients are all located at the l.h.s. of
            // the matrix whereas the "implicit" coefficients on the
            // coupled boundaries are all created as if the
            // coefficient contribution is of a source-kind (i.e. they
            // have a sign as if they are on the r.h.s. of the matrix.
            // To compensate for this, it is necessary to turn the
            // sign of the contribution.
            {
                FieldField<Field, scalar> mBouCoeffs = -bouCoeffs;

                m.initMatrixInterfaces
                (
                    mBouCoeffs,
                    interfaces,
                    psi,
                    bPrime,
                    cmpt
                );

                m.updateMatrixInterfaces
                (
                    mBouCoeffs,
                    interfaces,
                    psi,
                    bPrime,
                    cmpt
                );
            }

            for (register label cellI = 0; cellI < nCells; cellI++)
            {
                // lCell is equal to cellI
                scalar& curPsi = psiPtr[cellI];

                // Grab the accumulated neighbour side
                curPsi = bPrimePtr[cellI];

                // Start and end of this row
                label fStart = ownStartPtr[cellI];
                label fEnd = ownStartPtr[cellI + 1];

                // Accumulate the owner product side
                for (register label curFace = fStart; curFace < fEnd; curFace++)
                {
                    curPsi -= upperPtr[curFace]*psiPtr[uPtr[curFace]];
                }

                // Finish current psi
                curPsi /= diagPtr[cellI];

                // Distribute the neighbour side using current psi
                for (register label curFace = fStart; curFace < fEnd; curFace++)
                {
                    bPrimePtr[uPtr[curFace]] -= lowerPtr[curFace]*curPsi;
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "lduMatrix::smooth(scalarField&, const scalarField&, const label)"
        )   << "cannot solve incomplete matrix, no diagonal"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
