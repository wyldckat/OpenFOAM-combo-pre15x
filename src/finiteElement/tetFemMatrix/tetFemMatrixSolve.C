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
    Tetrahedral Finite Element matrix basic solvers.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

template<class Type>
lduMatrix::solverPerformance tetFemMatrix<Type>::solve()
{
#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve : solving tetFemMatrix<Type>"
            << endl;
    }
#   endif

    lduMatrix::solverPerformance solverPerfVec
    (
        "tetFemMatrix<Type>::solve",
        psi_.name()
    );

    // Add boundary source for gradient-type conditions
    addBoundarySourceDiag();

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    for
    (
        solvingComponent = 0;
        solvingComponent < Type::nComponents;
        solvingComponent++
    )
    {
        scalarField psiCmpt = psi_.internalField().component(solvingComponent);
        scalarField sourceCmpt = source_.component(solvingComponent);

        // Set component boundary conditions
        setComponentBoundaryConditions(solvingComponent, psiCmpt, sourceCmpt);

        // Add the coupling coefficients
        addCouplingCoeffs();
        addCouplingSource(sourceCmpt);

        // Prepare for coupled interface update
        FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());
        FieldField<Field, scalar> coupledIntCoeffs(psi_.boundaryField().size());
        lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

        forAll(psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            coupledBouCoeffs.hook
            (
                new scalarField(ptf.cutBouCoeffs(*this))
            );

            coupledIntCoeffs.hook
            (
                new scalarField(ptf.cutIntCoeffs(*this))
            );

            interfaces[patchI] = &psi_.boundaryField()[patchI];
        }

        eliminateCouplingCoeffs();

        lduMatrix::solverPerformance solverPerf = lduMatrix::solver::New
        (
            psi_.name() + pTraits<Type>::componentNames[solvingComponent],
            psiCmpt,
            *this,
            sourceCmpt,
            coupledBouCoeffs,
            coupledIntCoeffs,
            interfaces,
            solvingComponent,
            psi_.mesh().solver(psi_.name())
        )->solve();
        
        if
        (
            solverPerf.initialResidual() > solverPerfVec.initialResidual()
         && !solverPerf.singular()
        )
        {
            solverPerfVec = solverPerf;
        }

        psi_.internalField().replace(solvingComponent, psiCmpt);

        reconstructMatrix();
    }

#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve : correcting boundary conditions"
            << endl;
    }
#   endif

    psi_.correctBoundaryConditions();

    return solverPerfVec;
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > tetFemMatrix<Type>::residual()
{
    tmp<Field<Type> > tres(psi_.size());

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        scalarField PsiInternalCmpt = psi_.internalField().component(cmpt);
        scalarField sourceCmpt = source_.component(cmpt);

        setComponentBoundaryConditions(cmpt, PsiInternalCmpt, sourceCmpt);

        // Add the coupling coefficients
        addCouplingCoeffs();
        addCouplingSource(sourceCmpt);

        // Prepare for coupled interface update
        FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());
        lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

        forAll(psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];
            coupledBouCoeffs.hook(new scalarField(ptf.cutBouCoeffs(*this)));

            interfaces[patchI] = &psi_.boundaryField()[patchI];
        }

        eliminateCouplingCoeffs();

        tres().replace
        (
            cmpt,
            lduMatrix::residual
            (
                psi_.internalField().component(cmpt),
                sourceCmpt,
                coupledBouCoeffs,
                interfaces,
                cmpt
            )
        );

        reconstructMatrix();
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
