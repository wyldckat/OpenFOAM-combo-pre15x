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
    Top level solver class which selects the solver relevant to the particular
    matrix structure.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"
#include "diagonalSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define the constructor function hash tables

defineRunTimeSelectionTable
(
    lduMatrix::solver,
    symMatrix
);

defineRunTimeSelectionTable
(
    lduMatrix::solver,
    asymMatrix
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<lduMatrix::solver> lduMatrix::solver::New
(
    const word& fieldName,
    scalarField& psi,
    lduMatrix& matrix,
    const scalarField& source,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const direction cmpt,
    Istream& solverData
)
{
    word solverName(solverData);

    if (matrix.diagonal())
    {
        return autoPtr<lduMatrix::solver>
        (
            new diagonalSolver
            (
                fieldName,
                psi,
                matrix,
                source,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces,
                cmpt
            )
        );
    }
    else if (matrix.symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(solverName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::solver::New(const fvMesh&, Istream&)",
                solverData
            )   << "Unknown symmetric matrix solver " << solverName
                << endl << endl
                << "Valid symmetric matrix solvers are :" << endl
                << symMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::solver>
        (
            constructorIter()
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
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(solverName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::solver::New(const fvMesh&, Istream&)",
                solverData
            )   << "Unknown asymmetric matrix solver " << solverName
                << endl << endl
                << "Valid asymmetric matrix solvers are :" << endl
                << asymMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::solver>
        (
            constructorIter()
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
        );
    }
    else
    {
        FatalErrorIn
        (
            "lduMatrix::solver::New"
            "(lduMatrix& matrix, const word& solverName)"
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalError);

        return autoPtr<lduMatrix::solver>(NULL);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
