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
    Finite-Area matrix basic solvers.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<class Type>
void faMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    internalCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

    boundaryCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]*value;
}


template<class Type>
lduMatrix::solverPerformance faMatrix<Type>::solve(Istream& solverControls)
{
#   ifdef DEBUGfaMatrix
    if (debug)
    {
        Info<< "faMatrix<Type>::solve(Istream& solverControls) : "
               "solving faMatrix<Type>"
            << endl;
    }
#   endif

    lduMatrix::solverPerformance solverPerfaec
    (
        "faMatrix<Type>::solve",
        psi_.name()
    );

    scalarField saveDiag = diag();

    Field<Type> source = source_;
    addBoundarySource(source);

    for
    (
        solvingComponent=0;
        solvingComponent<Type::nComponents;
        solvingComponent++
    )
    {
        // copy field and source

        scalarField psiCmpt = psi_.internalField().component(solvingComponent);
        addBoundaryDiag(diag(), solvingComponent);

        scalarField sourceCmpt = source.component(solvingComponent);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(solvingComponent)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(solvingComponent)
        );

        lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

        forAll (bouCoeffsCmpt, patchI)
        {
            interfaces[patchI] = &psi_.boundaryField()[patchI];
        }

        updateMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            solvingComponent
        );

        lduMatrix::solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi_.name() + pTraits<Type>::componentNames[solvingComponent],
            psiCmpt,
            *this,
            sourceCmpt,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solvingComponent,
            solverControls.rewind()
        )->solve();

        if
        (
            solverPerf.initialResidual() > this->solverPerfVec.initialResidual()
         && !solverPerf.singular()
        )
        {
            this->solverPerfVec = solverPerf;
        }

        psi_.internalField().replace(solvingComponent, psiCmpt);
        diag() = saveDiag;
    }

    psi_.correctBoundaryConditions();

    return this->solverPerfVec;
}


template<class Type>
lduMatrix::solverPerformance faMatrix<Type>::solve()
{
    return solve(psi_.mesh().solver(psi_.name()));
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > faMatrix<Type>::residual() const
{
    tmp<Field<Type> > tres(source_);
    Field<Type>& res = tres();

    addBoundarySource(res);

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.internalField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(solvingComponent)
        );

        lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

        forAll (bouCoeffsCmpt, patchI)
        {
            interfaces[patchI] = &psi_.boundaryField()[patchI];
        }

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                interfaces,
                cmpt
            )
        );
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
