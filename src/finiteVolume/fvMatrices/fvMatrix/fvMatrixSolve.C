/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<class Type>
void fvMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            internalCoeffs_[patchi][facei].component(cmpt) +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

            boundaryCoeffs_[patchi][facei].component(cmpt) +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]*value;
        }
    }
}


template<class Type>
lduMatrix::solverPerformance fvMatrix<Type>::solve(Istream& solverControls)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::solve(Istream& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    lduMatrix::solverPerformance solverPerfVec
    (
        "fvMatrix<Type>::solve",
        psi_.name()
    );

    scalarField saveDiag = diag();

    Field<Type> source = source_;
    addBoundarySource(source, false);

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh().directions(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        // copy field and source

        scalarField psiCmpt = psi_.internalField().component(cmpt);
        addBoundaryDiag(diag(), cmpt);

        scalarField sourceCmpt = source.component(cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        lduInterfaceFieldPtrsList interfaces = 
            psi_.boundaryField().interfaces();

        lduMatrix::solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls.rewind()
        )->solve(psiCmpt, sourceCmpt, cmpt);

        solverPerf.print();

        if
        (
            solverPerf.initialResidual() > solverPerfVec.initialResidual()
         && !solverPerf.singular()
        )
        {
            solverPerfVec = solverPerf;
        }

        psi_.internalField().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi_.correctBoundaryConditions();

    return solverPerfVec;
}


template<class Type>
autoPtr<typename fvMatrix<Type>::fvSolver> fvMatrix<Type>::solver()
{
    return solver(psi_.mesh().solver(psi_.name()));
}

template<class Type>
lduMatrix::solverPerformance fvMatrix<Type>::fvSolver::solve()
{
    return solve(psi_.mesh().solver(psi_.name()));
}


template<class Type>
lduMatrix::solverPerformance fvMatrix<Type>::solve()
{
    return solve(psi_.mesh().solver(psi_.name()));
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > fvMatrix<Type>::residual() const
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
            boundaryCoeffs_.component(cmpt)
        );

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                psi_.boundaryField().interfaces(),
                cmpt
            )
        );
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
