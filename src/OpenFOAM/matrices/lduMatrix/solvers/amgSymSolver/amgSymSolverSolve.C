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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduMatrix::solverPerformance amgSymSolver::solve()
{
    // If no levels are created, bail out by calling ICCG
    if (matrixLevels_.size() < 1)
    {
#       ifdef FULLDEBUG
        if (lduMatrix::debug >= 2)
        {
        Info
            << "lduMatrix::solverPerformance amgSymSolver::solve : "
            << "No coarse levels. Matrix too small for AMG; "
            << "Reverting to CG."
            << endl;
        }
#       endif

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
    scalar psiRef = average(psi_);

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

    // Calculate initial residual field
    scalarField fineResidual(source_ - Apsi);

    // Calculate normalisation factor
    scalar normFactor =
        gSum(mag(Apsi - ApsiRef) + mag(source_ - ApsiRef)) + lduMatrix::small_;

    // Calculate residual magnitude
    solverPerf.initialResidual() = gSum(mag(fineResidual))/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
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

#       ifdef FULLDEBUG
        if (lduMatrix::debug >= 3)
        {
            Info<< "preparatory operations " << cpu_.cpuTimeIncrement() << endl;
        }
#       endif

        do
        {
            if (matrixLevels_.size() > 0)
            {
                // Restrict residual for the next level up
                coarseSources[0] = restrictField(fineResidual, 0);

                // Solver V-cycle

                // Residual restriction (going to coarser levels)
                for (label levelI = 0; levelI < topLevel; levelI++)
                {
                    // Clean up the field
                    coarseFields[levelI] = 0;

                    // Residual is equal to source
                    coarseSources[levelI + 1] = restrictField
                    (
                        coarseSources[levelI],
                        levelI + 1
                    );
                }

                // Top level
#               ifdef FULLDEBUG
                if (lduMatrix::debug >= 3)
                {
                    Info << "restriction " << cpu_.cpuTimeIncrement() << endl;
                }
#               endif

                // clean up the field
                coarseFields[topLevel] = 0;

                matrixLevels_[topLevel].quietOperation() = true;

                ICCG
                (
                    "topLevelCorr",
                    coarseFields[topLevel],
                    matrixLevels_[topLevel],
                    coarseSources[topLevel],
                    interfaceCoeffs_[topLevel],
                    interfaceCoeffs_[topLevel],
                    interfaceLevels_[topLevel],
                    cmpt_,
                    1e-12
                ).solve();

#               ifdef FULLDEBUG
                if (lduMatrix::debug >= 2)
                {
                    Info << "Scaling factors: ";
                }
#               endif

                // Post-smoothing (going to finer levels)
                for (label levelI = topLevel - 1; levelI >= 0; levelI--)
                {
                    // Note:
                    // This form of scaling assumes there is no pre-smoothing
                    // For full form see scaling for the finest mesh
                    // 

                    coarseFields[levelI] = prolongField
                    (
                        coarseFields[levelI + 1],
                        levelI + 1
                    );
                    
#                   ifdef FULLDEBUG
                    if (lduMatrix::debug >= 3)
                    {
                        Info << "prolongation " << cpu_.cpuTimeIncrement()
                            << endl;
                    }
#                   endif

                    // Calculate scaling factor
                    scalarField ACf(coarseFields[levelI].size());
                    matrixLevels_[levelI].Amul
                    (
                        ACf,
                        coarseFields[levelI],
                        interfaceCoeffs_[levelI],
                        interfaceLevels_[levelI],
                        cmpt_
                    );

                    scalar scalingFactor =
                        gSum(coarseSources[levelI]*coarseFields[levelI])
                       /stabilise(gSum(coarseFields[levelI]*ACf), VSMALL);

#                   ifdef FULLDEBUG
                    if (lduMatrix::debug >= 2)
                    {
                        Info << scalingFactor << " ";
                    }
#                   endif

                    coarseFields[levelI] *= scalingFactor;

#                   ifdef FULLDEBUG
                    if (lduMatrix::debug >= 3)
                    {
                        Info << "scaling " << cpu_.cpuTimeIncrement() << endl;
                    }
#                   endif

                    // Smooth the solution
                    matrixLevels_[levelI].smooth
                    (
                        coarseFields[levelI],
                        coarseSources[levelI],
                        interfaceCoeffs_[levelI],
                        interfaceLevels_[levelI],
                        cmpt_,
                        nPostSweeps_
                    );

#                   ifdef FULLDEBUG
                    if (lduMatrix::debug >= 3)
                    {
                        Info << "smoothing " << cpu_.cpuTimeIncrement() << endl;
                    }
#                   endif
                }

                // Prolong the finest level correction
                scalarField fineCorrection = prolongField(coarseFields[0], 0);

#               ifdef FULLDEBUG
                if (lduMatrix::debug >= 3)
                {
                    Info << "finest prolongation " << cpu_.cpuTimeIncrement()
                        << endl;
                }
#               endif

                // Calculate fine scaling factor
                scalarField AFine(psi_.size());
                matrix_.Amul
                (
                    AFine,
                    fineCorrection,
                    coupleBouCoeffs_,
                    interfaces_,
                    cmpt_
                );

                scalar fineScalingFactor =
                    gSum(fineResidual*fineCorrection)
                   /stabilise(gSum(fineCorrection*AFine), VSMALL);

#               ifdef FULLDEBUG
                if (lduMatrix::debug >= 2)
                {
                    Info << fineScalingFactor << " ";
                }
#               endif

                psi_ += fineScalingFactor*fineCorrection;
            }

#           ifdef FULLDEBUG
            if (lduMatrix::debug >= 3)
            {
                Info << "finest scaling " << cpu_.cpuTimeIncrement() << endl;
            }
#           endif

            // smooth fine matrix
            matrix_.smooth
            (
                psi_,
                source_,
                coupleBouCoeffs_,
                interfaces_,
                cmpt_,
                nBottomSweeps_
            );

#           ifdef FULLDEBUG
            if (lduMatrix::debug >= 3)
            {
                Info << "finest smoothing " << cpu_.cpuTimeIncrement() << endl;
            }
#           endif

            // Calculate fine residual
            fineResidual = matrix_.residual
            (
                psi_,
                source_,
                coupleBouCoeffs_,
                interfaces_,
                cmpt_
            );

            solverPerf.finalResidual() = gSum(mag(fineResidual))/normFactor;

#           ifdef FULLDEBUG
            if (lduMatrix::debug >= 3)
            {
                Info << "finest residual " << cpu_.cpuTimeIncrement() << endl;
            }
#           endif
        } while
        (
            solverPerf.nIterations()++ < maxCycles_
         && !(solverPerf.checkConvergence(tolerance_, relTol_))
        );

        solverPerf.print();
    }

    return solverPerf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
