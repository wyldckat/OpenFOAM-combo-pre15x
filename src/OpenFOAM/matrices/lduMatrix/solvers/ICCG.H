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

Class
    ICCG

Description

SourceFiles
    ICCG.C

\*---------------------------------------------------------------------------*/

#ifndef ICCG_H
#define ICCG_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ICCG Declaration
\*---------------------------------------------------------------------------*/

class ICCG
:
    public lduMatrix::solver
{
private:

    // Static data

        //- Maximum number of iterations in the solver
        static label maxIter_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ICCG(const ICCG&);

        //- Disallow default bitwise assignment
        void operator=(const ICCG&);


protected:

    // Protected data

        scalar tolerance_;
        scalar relTol_;


public:

    //- Runtime type information
#   ifndef ICCGname
#   define ICCGname "ICCG"
#   endif
    TypeName(ICCGname);
#   undef ICCGname


    // Constructors

        //- Construct from components
        ICCG
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
            const scalar relTol = 0.0
        );

        //- Construct from matrix components and solver data stream
        ICCG
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
        );


    // Member Functions

        //- Solve the matrix with this solver
        lduMatrix::solverPerformance solve();

        //- Set maximum number of iterations
        static label setMaxIter(const label mIter)
        {
            label oldMaxIter = maxIter_;
            maxIter_ = mIter;

            return oldMaxIter;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
