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
    This is an optimised BICCG solver for asymmetric lduMatrices using
    LU-decomposion preconditioning.

\*---------------------------------------------------------------------------*/

#include "BICCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(BICCG, 0);

label BICCG::maxIter_ = 1000;

lduMatrix::solver::addasymMatrixConstructorToTable<BICCG>
    addBICCGAsymMatrixConstructorToTable_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
BICCG::BICCG
(
    const word& fieldName,
    scalarField& psi,
    const lduMatrix& matrix,
    const scalarField& source,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt,
    Istream& solverData
)
:
    ICCG
    (
        fieldName,
        psi,
        matrix,
        source,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        cmpt,
        solverData
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduMatrix::solverPerformance BICCG::solve()
{
    // --- Setup class containing solver performance data
    lduMatrix::solverPerformance solverPerf(typeName, fieldName_);


    register label cell, nCells = psi_.size();

#   define forField for (cell=0; cell<nCells; cell++)

#   ifndef DONT_IC_PRECONDITION

    register label face, sface;
    register label nFaces = matrix_.upper().size(), nFacesM1 = nFaces - 1;

#   define forFaces for (face=0; face<nFaces; face++)
#   define forFacesReverse for (face=nFacesM1; face>=0; face--)

    const label* const restrict uPtr = matrix_.lduAddr().upperAddr().begin();
    const label* const restrict lPtr = matrix_.lduAddr().lowerAddr().begin();
    const label* const restrict losortPtr =
        matrix_.lduAddr().losortAddr().begin();

    // Matrix is symmetric with the elements stored only in the upper triangle
    const scalar* const restrict upperPtr = matrix_.upper().begin();
    const scalar* const restrict lowerPtr = matrix_.lower().begin();

#   endif

    scalar* restrict psiPtr = psi_.begin();

    scalarField pA(psi_.size(), 0.0);
    scalar* restrict pAPtr = pA.begin();

    scalarField pT(psi_.size(), 0.0);
    scalar* restrict pTPtr = pT.begin();

    scalarField wA(nCells);
    scalar* restrict wAPtr = wA.begin();

    scalarField wT(nCells);
    scalar* restrict wTPtr = wT.begin();

    scalar wArT = matrix_.great_;
    scalar wArTold = wArT;

    scalar alpha, beta, wApT;


    // --- make a copy of the diagonal for preconditioning
    scalarField dD(matrix_.diag());
    scalar* restrict dDPtr = dD.begin();


    // --- Calculate reference value of psi_
    scalar psiRef = gAverage(psi_);

    // --- Calculate A.psi_, T.psi_ and A.psiRef
    matrix_.Amul(wA, psi_, coupleBouCoeffs_, interfaces_, cmpt_);
    matrix_.Tmul(wT, psi_, coupleIntCoeffs_, interfaces_, cmpt_);
    matrix_.Amul
    (
        pA,
        scalarField(psi_.size(), psiRef),
        coupleBouCoeffs_,
        interfaces_,
        cmpt_
    );

    // --- Calculate initial residual field
    scalarField rA(source_ - wA);
    scalarField rT(source_ - wT);
    scalar* restrict rAPtr = rA.begin();
    scalar* restrict rTPtr = rT.begin();

    // --- Calculate normalisation factor
    scalar normFactor = gSum(mag(wA - pA) + mag(source_ - pA)) + matrix_.small_;

    // --- reset pA after it's temporary as A.psiRef
    pA = 0.0;

    // --- Calculate residual magnitude
    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }


    // --- Check convergence, solve if not converged

    if (!solverPerf.checkConvergence(tolerance_, relTol_))
    {
        // --- Precondition diagonal as part of IC preconditioning
#       ifndef DONT_IC_PRECONDITION
        forFaces
        {
            dDPtr[uPtr[face]] -=
                upperPtr[face]*lowerPtr[face]/dDPtr[lPtr[face]];
        }
#       endif

        // --- Generate reciprocal diagonal for future use
        // --- This is only useful on systems where * is cheaper than /
        forField
        {
            dDPtr[cell] = 1.0/dDPtr[cell];
        }


        // --- Iteration loop

        do
        {
            wArTold = wArT;

            // --- Initialize diagonal-preconditioned residuals
            forField
            {
                wAPtr[cell] = dDPtr[cell]*rAPtr[cell];
                wTPtr[cell] = dDPtr[cell]*rTPtr[cell];
            }

            // --- IC preconditioning
#           ifndef DONT_IC_PRECONDITION
            forFaces
            {
                sface = losortPtr[face];
                wAPtr[uPtr[sface]] -=
                    dDPtr[uPtr[sface]]*lowerPtr[sface]*wAPtr[lPtr[sface]];

                wTPtr[uPtr[face]] -=
                    dDPtr[uPtr[face]]*upperPtr[face]*wTPtr[lPtr[face]];
            }

            forFacesReverse
            {
                wAPtr[lPtr[face]] -=
                    dDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];

                sface = losortPtr[face];
                wTPtr[lPtr[sface]] -=
                    dDPtr[lPtr[sface]]*lowerPtr[sface]*wTPtr[uPtr[sface]];
            }
#           endif

            // --- Update search directions:

            wArT = gSumProd(wA, rT);

            beta = wArT/wArTold;

            forField
            {
                pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                pTPtr[cell] = wTPtr[cell] + beta*pTPtr[cell];
            }


            // --- Update preconditioned residuals

            matrix_.Amul(wA, pA, coupleBouCoeffs_, interfaces_, cmpt_);
            matrix_.Tmul(wT, pT, coupleIntCoeffs_, interfaces_, cmpt_);

            wApT = gSumProd(wA, pT);


            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApT)/normFactor)) break;


            // --- Update solution and raw residuals:

            alpha = wArT/wApT;

            forField
            {
                psiPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
                rTPtr[cell] -= alpha*wTPtr[cell];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
          //solverPerf.finalResidual() =
          //    gSumMag
          //    (
          //          source
          //        - Amul(psi_, ?, coupleBouCoeffs, interfaces, cmpt)
          //      )/normFactor;

        } while
        (
            solverPerf.nIterations()++ < maxIter_
        && !(solverPerf.checkConvergence(tolerance_, relTol_))
        );

        solverPerf.print(matrix_.quietOperation());
    }

    return solverPerf;
}

#undef forField
#undef forFaces


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
