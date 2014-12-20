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
    Gauss-Seidel solver for symmetric and assymetric matrices.  In
    order to improve efficiency, the residual is evaluated after every
    nSweeps sweeps.

\*---------------------------------------------------------------------------*/

#include "GaussSeidel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(GaussSeidel, 0);

label GaussSeidel::maxIter_ = 5000;

lduMatrix::solver::addsymMatrixConstructorToTable<GaussSeidel>
    addGaussSeidelSymMatrixConstructorToTable_;

lduMatrix::solver::addasymMatrixConstructorToTable<GaussSeidel>
    addGaussSeidelAsymMatrixConstructorToTable_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
GaussSeidel::GaussSeidel
(
    const word& fieldName,
    scalarField& psi,
    const lduMatrix& matrix,
    const scalarField& source,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt,
    const label nSweeps,
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
    relTol_(relTol),
    nSweeps_(nSweeps)
{}


//- Construct from matrix and solver data stream
GaussSeidel::GaussSeidel
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
    relTol_(readScalar(solverData)),
    nSweeps_(readLabel(solverData))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduMatrix::solverPerformance GaussSeidel::solve()
{
    // --- Setup class containing solver performance data
    lduMatrix::solverPerformance solverPerf(typeName, fieldName_);

    scalar normFactor = 0;

    {
        scalarField pA(psi_.size());
        scalarField wA(psi_.size());

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

        normFactor = gSum(mag(wA - pA) + mag(source_ - pA)) + matrix_.small_;

        // --- Calculate residual magnitude
        solverPerf.initialResidual() = gSumMag(source_ - wA)/normFactor;
        solverPerf.finalResidual() = solverPerf.initialResidual();
    }

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }


    // --- Check convergence, solve if not converged

    if (!solverPerf.checkConvergence(tolerance_, relTol_))
    {
        // --- Iteration loop

        do
        {
            matrix_.smooth
            (
                psi_,
                source_,
                coupleBouCoeffs_,
                interfaces_,
                cmpt_,
                nSweeps_
            );

            solverPerf.finalResidual() = gSumMag
            (
                matrix_.residual
                (
                    psi_,
                    source_,
                    coupleBouCoeffs_,
                    interfaces_,
                    cmpt_
                )
            )/normFactor;

        } while
        (
            (solverPerf.nIterations() += nSweeps_) < maxIter_
        && !(solverPerf.checkConvergence(tolerance_, relTol_))
        );

        solverPerf.print(matrix_.quietOperation());
    }

    return solverPerf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
