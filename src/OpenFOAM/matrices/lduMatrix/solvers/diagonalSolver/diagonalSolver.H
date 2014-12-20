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

Class
    diagonalSolver

Description

SourceFiles
    diagonalSolver.C

\*---------------------------------------------------------------------------*/

#ifndef diagonalSolver_H
#define diagonalSolver_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diagonalSolver Declaration
\*---------------------------------------------------------------------------*/

class diagonalSolver
:
    public lduMatrix::solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        diagonalSolver(const diagonalSolver&);

        //- Disallow default bitwise assignment
        void operator=(const diagonalSolver&);


public:

    //- Runtime type information
    TypeName("diagonal");


    // Constructors

        //- Construct from matrix
        diagonalSolver
        (
            const word& fieldName,
            scalarField& psi,
            const lduMatrix& matrix,
            const scalarField& source,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduCoupledInterfacePtrsList& interfaces,
            const direction cmpt
        );


    // Member Functions

        //- Solve the matrix with this solver
        lduMatrix::solverPerformance solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
