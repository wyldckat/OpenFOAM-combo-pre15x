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
    This is an optimised ICCG solver for symmetric lduMatrices using
    LU-decomposion preconditioning.

\*---------------------------------------------------------------------------*/

#include "ICCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ICCG, 0);

label ICCG::maxIter_ = 5000;

lduMatrix::solver::addsymMatrixConstructorToTable<ICCG>
    addICCGSymMatrixConstructorToTable_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
ICCG::ICCG
(
    const word& fieldName,
    scalarField& psi,
    const lduMatrix& matrix,
    const scalarField& source,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt,
    const scalar tolerance,
    const scalar relTol
)
:
    lduMatrix::solver
    (
        fieldName,
        psi,
        matrix,
        source,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        cmpt
    ),
    tolerance_(tolerance),
    relTol_(relTol)
{}


//- Construct from matrix and solver data stream
ICCG::ICCG
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
    lduMatrix::solver
    (
        fieldName,
        psi,
        matrix,
        source,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        cmpt
    ),
    tolerance_(readScalar(solverData)),
    relTol_(readScalar(solverData))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduMatrix::solverPerformance ICCG::solve()
{
    // --- Setup class containing solver performance data
    lduMatrix::solverPerformance solverPerf(typeName, fieldName_);

    register label cell, nCells = psi_.size();

#   define forField for (cell=0; cell<nCells; cell++)

#   ifndef DONT_IC_PRECONDITION

    register label face, nFaces = matrix_.upper().size(), nFacesM1 = nFaces - 1;

    const label* const restrict uPtr = matrix_.lduAddr().upperAddr().begin();
    const label* const restrict lPtr = matrix_.lduAddr().lowerAddr().begin();

    // Matrix is symmetric with the elements stored only in the upper triangle
    const scalar* const restrict upperPtr = matrix_.upper().begin();

#   define forFaces for (face=0; face<nFaces; face++)
#   define forFacesReverse for (face=nFacesM1; face>=0; face--)

#   endif

    scalar* restrict psiPtr = psi_.begin();

    scalarField pA(psi_.size());
    scalar* restrict pAPtr = pA.begin();

    scalarField wA(nCells);
    scalar* restrict wAPtr = wA.begin();

    scalar wArA = matrix_.great_;
    scalar wArAold = wArA;

    scalar alpha, beta, wApA;


    // --- make a copy of the diagonal for preconditioning
    scalarField dD(matrix_.diag());
    scalar* restrict dDPtr = dD.begin();


    // --- Calculate reference value of psi
    scalar psiRef = gAverage(psi_);

    // --- Calculate A.psi and A.psiRef
    matrix_.Amul(wA, psi_, coupleBouCoeffs_, interfaces_, cmpt_);
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
    scalar* restrict rAPtr = rA.begin();

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
            dDPtr[uPtr[face]] -= sqr(upperPtr[face])/dDPtr[lPtr[face]];
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
            wArAold = wArA;

            // --- Initialize diagonal-preconditioned residuals
            forField
            {
                wAPtr[cell] = dDPtr[cell]*rAPtr[cell];
            }

            // --- IC preconditioning
#           ifndef DONT_IC_PRECONDITION
            forFaces
            {
                wAPtr[uPtr[face]] -=
                    dDPtr[uPtr[face]]*upperPtr[face]*wAPtr[lPtr[face]];
            }

            forFacesReverse
            {
                wAPtr[lPtr[face]] -=
                    dDPtr[lPtr[face]]*upperPtr[face]*wAPtr[uPtr[face]];
            }
#           endif

            // --- Update search directions:

            wArA = gSumProd(wA, rA);

            beta = wArA/wArAold;

            forField
            {
                pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
            }


            // --- Update preconditioned residuals

            matrix_.Amul(wA, pA, coupleBouCoeffs_, interfaces_, cmpt_);

            wApA = gSumProd(wA, pA);


            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


            // --- Update solution and raw residuals:

            alpha = wArA/wApA;

            forField
            {
                psiPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            //scalarField hmm(psi_.size());
            //Amul(hmm, psi_, coupleBouCoeffs, interfaces, cmpt);
            //solverPerf.finalResidual() = gSumMag(source - hmm)/normFactor;
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
