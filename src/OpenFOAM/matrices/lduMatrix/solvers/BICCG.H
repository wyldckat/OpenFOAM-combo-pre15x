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
    BICCG

Description

SourceFiles
    BICCG.C

\*---------------------------------------------------------------------------*/

#ifndef BICCG_H
#define BICCG_H

#include "ICCG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BICCG Declaration
\*---------------------------------------------------------------------------*/

class BICCG
:
    public ICCG
{
private:

    // Static data

        //- Maximum number of iterations in the solver
        static label maxIter_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BICCG(const BICCG&);

        //- Disallow default bitwise assignment
        void operator=(const BICCG&);


public:

    //- Runtime type information
#   ifndef BICCGname
#   define BICCGname "BICCG"
#   endif
    TypeName(BICCGname);
#   undef BICCGname


    // Constructors

        //- Construct from matrix and solver data stream
        BICCG
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
