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
    simpleMatrix

Description

SourceFiles
    simpleMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMatrix_H
#define simpleMatrix_H

#include "Matrix.H"
#include "Field.H"
#include "scalar.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T>
class simpleMatrix;

template<class T>
simpleMatrix<T> operator+
(
    const simpleMatrix<T>&,
    const simpleMatrix<T>&
);

template<class T>
simpleMatrix<T> operator-
(
    const simpleMatrix<T>&,
    const simpleMatrix<T>&
);

template<class T>
simpleMatrix<T> operator*
(
    const scalar,
    const simpleMatrix<T>&
);

template<class T>
Ostream& operator<<
(
    Ostream&,
    const simpleMatrix<T>&
);


/*---------------------------------------------------------------------------*\
                           Class simpleMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class simpleMatrix
{
    // Private data

        Matrix<scalar> matrix_;
        Field<T> source_;


public:

    // Constructors

        //- Construct given size
        simpleMatrix(const label);

        //- Construct from components
        simpleMatrix(const Matrix<scalar>&, const Field<T>&);

        //- Construct from Istream
        simpleMatrix(Istream&);

        //- Construct as copy
        simpleMatrix(const simpleMatrix<T>&);


    // Member Functions

        // Access

            Matrix<scalar>& matrix()
            {
                return matrix_;
            }

            Field<T>& source()
            {
                return source_;
            }

            const Matrix<scalar>& matrix() const
            {
                return matrix_;
            }

            const Field<T>& source() const
            {
                return source_;
            }


        //- Solve the matrix using Gaussian elimination with pivoting
        static void solve(Matrix<scalar>& matrix, Field<T>& source);

        //- Solve the matrix using Gaussian elimination with pivoting
        //  and return the solution
        Field<T> solve() const;


        //- LU decompose the matrix with pivoting
        static void LUDecompose
        (
            Matrix<scalar>& matrix,
            labelList& pivotIndices
        );

        //- LU back-substitution with given source, returning the solution
        //  in the source
        static void LUBacksubstitute
        (
            const Matrix<scalar>& luMmatrix,
            const labelList& pivotIndices,
            Field<T>& source
        );

        //- Solve the matrix using LU decomposition with pivoting
        //  returning the LU form of the matrix and the solution in the source
        static void LUsolve(Matrix<scalar>& matrix, Field<T>& source);

        //- Solve the matrix using LU decomposition with pivoting
        //  and return the solution
        Field<T> LUsolve() const;


    // Member Operators

        void operator=(const simpleMatrix<T>&);


    // Friend Operators

        friend simpleMatrix<T> operator+ <T>
        (
            const simpleMatrix<T>&,
            const simpleMatrix<T>&
        );

        friend simpleMatrix<T> operator- <T>
        (
            const simpleMatrix<T>&,
            const simpleMatrix<T>&
        );

        friend simpleMatrix<T> operator* <T>
        (
            const scalar,
            const simpleMatrix<T>&
        );


    // Ostream Operator

        friend Ostream& operator<< <T>
        (
            Ostream&,
            const simpleMatrix<T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "simpleMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
