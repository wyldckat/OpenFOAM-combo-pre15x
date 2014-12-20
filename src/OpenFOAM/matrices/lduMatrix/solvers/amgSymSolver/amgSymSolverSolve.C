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
    Agglomerated algebraic multigrid solver tuned for the FV elliptic matrices.

\*---------------------------------------------------------------------------*/

#include "amgSymSolver.H"
#include "ICCG.H"
#include "GaussSeidel.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduMatrix::solverPerformance amgSymSolver::solve()
{
    // If no levels are created, bail out by calling ICCG
    if (matrixLevels_.size() < 1)
    {
        if (debug >= 2)
        {
            Pout<< "lduMatrix::solverPerformance amgSymSolver::solve : "
                << "No coarse levels. Matrix too small for AMG; "
                << "Reverting to CG."
                << endl;
        }

        return ICCG
        (
            fieldName_,
            psi_,
            matrix_,
            source_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaces_,
            cmpt_,
            tolerance_,
            relTol_
        ).solve();
    }


    // Setup class containing solver performance data
    lduMatrix::solverPerformance solverPerf(typeName, fieldName_);

    // Calculate reference value of psi
    scalar psiRef = gAverage(psi_);

    // Calculate A.psi and A.psiRef
    scalarField Apsi(psi_.size());
    matrix_.Amul(Apsi, psi_, coupleBouCoeffs_, interfaces_, cmpt_);

    scalarField ApsiRef(psi_.size());
    matrix_.Amul
    (
        ApsiRef,
        scalarField(psi_.size(), psiRef),
        coupleBouCoeffs_,
        interfaces_,
        cmpt_
    );
        
    // Calculate normalisation factor
    scalar normFactor =
        gSum(mag(Apsi - ApsiRef) + mag(source_ - ApsiRef))
      + lduMatrix::small_;

    // Calculate initial residual field
    scalarField fineResidual(source_ - Apsi);

    // Calculate residual magnitude
    solverPerf.initialResidual() = gSumMag(fineResidual)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (debug >= 2)
    {
        Pout<< "   Normalisation factor = " << normFactor << endl;
    }


    // Check convergence, solve if not converged

    if (!solverPerf.checkConvergence(tolerance_, relTol_))
    {
        // Create coarse fields
        PtrList<scalarField> coarseFields(matrixLevels_.size());

        // Create coarse Sources
        PtrList<scalarField> coarseSources(matrixLevels_.size());

        forAll (matrixLevels_, levelI)
        {
            coarseFields.hook
            (
                new scalarField(addrLevels_[levelI].size(), 0.0)
            );

            coarseSources.hook
            (
                new scalarField(addrLevels_[levelI].size(), 0.0)
            );
        }

        const label topLevel = matrixLevels_.size() - 1;

        if (debug >= 3)
        {
            Pout<< "preparatory operations " << cpu_.cpuTimeIncrement() << endl;
        }

        do
        {
            if (matrixLevels_.size() > 0)
            {
                // Restrict residual for the next level up
                restrictField(coarseSources[0], fineResidual, 0);

                // Solver V-cycle

                // Residual restriction (going to coarser levels)
                for (label levelI = 0; levelI < topLevel; levelI++)
                {
                    // Residual is equal to source
                    restrictField
                    (
                        coarseSources[levelI + 1],
                        coarseSources[levelI],
                        levelI + 1
                    );
                }

                // Top level

                if (debug >= 3)
                {
                    Pout<< "restriction " << cpu_.cpuTimeIncrement() << endl;
                }

                // clean up the field
                coarseFields[topLevel] = 0;

                lduMatrix::solverPerformance coarseSolverPerf = ICCG
                (
                    "topLevelCorr",
                    coarseFields[topLevel],
                    matrixLevels_[topLevel],
                    coarseSources[topLevel],
                    interfaceCoeffs_[topLevel],
                    interfaceCoeffs_[topLevel],
                    interfaceLevels_[topLevel],
                    cmpt_,
                    tolerance_,
                    relTol_
                ).solve();

                if (debug >= 2)
                {
                    coarseSolverPerf.print();

                    Pout<< "Scaling factors: ";
                }

                // Post-smoothing (going to finer levels)
                for (label levelI = topLevel - 1; levelI >= 0; levelI--)
                {
                    // Note:
                    // This form of scaling assumes there is no pre-smoothing
                    // For full form see scaling for the finest mesh

                    prolongField
                    (
                        coarseFields[levelI],
                        coarseFields[levelI + 1],
                        levelI + 1
                    );

                    if (debug >= 3)
                    {
                        Pout<< "prolongation " << cpu_.cpuTimeIncrement()
                            << endl;
                    }

                    // Calculate scaling factor
                    scalarField::subField ACf
                    (
                        Apsi,
                        coarseFields[levelI].size()
                    );

                    matrixLevels_[levelI].Amul
                    (
                        reinterpret_cast<scalarField&>(ACf),
                        coarseFields[levelI],
                        interfaceCoeffs_[levelI],
                        interfaceLevels_[levelI],
                        cmpt_
                    );

                    scalar sf = scalingFactor
                    (
                        coarseFields[levelI],
                        coarseSources[levelI],
                        ACf
                    );

                    if (debug >= 2)
                    {
                        Pout<< sf << " ";
                    }

                    coarseFields[levelI] *= sf;

                    if (debug >= 3)
                    {
                        Pout<< "scaling " << cpu_.cpuTimeIncrement() << endl;
                    }

                    // Smooth the solution
                    GaussSeidel::smooth
                    (
                        fieldName_,
                        coarseFields[levelI],
                        matrixLevels_[levelI],
                        coarseSources[levelI],
                        interfaceCoeffs_[levelI],
                        interfaceLevels_[levelI],
                        cmpt_,
                        nPostSweeps_
                    );

                    if (debug >= 3)
                    {
                        Pout<< "smoothing " << cpu_.cpuTimeIncrement() << endl;
                    }
                }

                // Prolong the finest level correction
                scalarField& fineCorrection = ApsiRef;
                prolongField(fineCorrection, coarseFields[0], 0);

                if (debug >= 3)
                {
                    Pout<< "finest prolongation " << cpu_.cpuTimeIncrement()
                        << endl;
                }

                // Calculate fine scaling factor
                matrix_.Amul
                (
                    Apsi,
                    fineCorrection,
                    coupleBouCoeffs_,
                    interfaces_,
                    cmpt_
                );

                scalar fsf = scalingFactor
                (
                    fineCorrection,
                    fineResidual,
                    Apsi
                );

                if (debug >= 2)
                {
                    Pout<< fsf << " ";
                }

                forAll(psi_, i)
                {
                    psi_[i] += fsf*fineCorrection[i];
                }
            }

            if (debug >= 3)
            {
                Pout<< "finest scaling " << cpu_.cpuTimeIncrement() << endl;
            }

            // smooth fine matrix
            GaussSeidel::smooth
            (
                fieldName_,
                psi_,
                matrix_,
                source_,
                coupleBouCoeffs_,
                interfaces_,
                cmpt_,
                nBottomSweeps_
            );

            if (debug >= 3)
            {
                Pout<< "finest smoothing " << cpu_.cpuTimeIncrement() << endl;
            }

            // Calculate fine residual
            matrix_.Amul(Apsi, psi_, coupleBouCoeffs_, interfaces_, cmpt_);
            fineResidual = source_;
            fineResidual -= Apsi;

            solverPerf.finalResidual() = gSumMag(fineResidual)/normFactor;

            if (debug >= 3)
            {
                Pout<< "finest residual " << cpu_.cpuTimeIncrement() << endl;
            }
        } while
        (
            solverPerf.nIterations()++ < maxCycles_
         && !(solverPerf.checkConvergence(tolerance_, relTol_))
        );
    }

    return solverPerf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
