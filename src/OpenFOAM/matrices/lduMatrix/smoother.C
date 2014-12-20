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
    Gauss-Seidel smoother with fixed number of sweeps.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void lduMatrix::smooth
(
    scalarField& psi,
    const scalarField& Source,
    const FieldField<Field, scalar>& bouCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt,
    const label nSweeps
) const
{
    if (diagonal())
    {
        psi = Source/diag();
    }
    else if (symmetric() || asymmetric())
    {
        scalarField dD = 1.0/diag();

        scalarField bPrime(psi.size());

        const unallocLabelList& u = lduAddr_.upperAddr();
        const unallocLabelList& ownStart = lduAddr_.ownerStartAddr();

        const scalarField& Lower = lower();
        const scalarField& Upper = upper();
        const label nCells = ownStart.size() - 1;

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
            updateMatrixInterfaces
            (
                bouCoeffs,
                interfaces,
                -psi,
                bPrime,
                cmpt
            );

            for (register label cellI = 0; cellI < nCells; cellI++)
            {
                register label curFace = ownStart[cellI];
                label fEnd = ownStart[cellI + 1];

                // lCell is equal to cellI

                scalar& curPsi = psi[cellI];

                // Grab the accumulated neighbour side
                curPsi = bPrime[cellI];

                // Accumulate the owner product side
                for (; curFace < fEnd; curFace++)
                {
                    register label uCell = u[curFace];
                    curPsi -= Upper[curFace]*psi[uCell];
                }

                // Finish current psi
                curPsi *= dD[cellI];

                // Distribute the neighbour side using current psi
                curFace = ownStart[cellI];

                for (; curFace < fEnd; curFace++)
                {
                    register label uCell = u[curFace];
                    bPrime[uCell] -= Lower[curFace]*curPsi;
                }
            }
        }
    }
    else
    {
        FatalError
            << "lduMatrix::smooth"
            << "(scalarField&, const scalarField&, const label) : "
            << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
