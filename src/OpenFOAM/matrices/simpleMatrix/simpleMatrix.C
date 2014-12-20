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


\*---------------------------------------------------------------------------*/

#include "simpleMatrix.H"
#include "scalarField.H"
#include "swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given size
template<class T>
simpleMatrix<T>::simpleMatrix(const label mSize)
:
    matrix_(mSize, mSize, 0.0),
    source_(mSize, pTraits<T>::zero)
{}


// Construct from components
template<class T>
simpleMatrix<T>::simpleMatrix
(
    const Matrix<scalar>& matrix,
    const Field<T>& source
)
:
    matrix_(matrix),
    source_(source)
{}


// Construct from Istream
template<class T>
simpleMatrix<T>::simpleMatrix(Istream& is)
:
    matrix_(is),
    source_(is)
{}


// Construct as copy
template<class T>
simpleMatrix<T>::simpleMatrix(const simpleMatrix& M)
:
    matrix_(M.matrix_),
    source_(M.source_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void simpleMatrix<T>::solve
(
    Matrix<scalar>& tmpMatrix,
    Field<T>& sourceSol
)
{
    label n = tmpMatrix.n();

    // Elimination
    for (register label i=0; i<n; i++)
    {
        label iMax = i;
        scalar largestCoeff = mag(tmpMatrix[iMax][i]);

        // Swap entries around to find a good pivot
        for (register label j=i+1; j<n; j++)
        {
            if (mag(tmpMatrix[j][i]) > largestCoeff)
            {
                iMax = j;
                largestCoeff = mag(tmpMatrix[iMax][i]);
            }
        }

        if (i != iMax)
        {
            //Info<< "Pivoted on " << i << " " << iMax << endl;

            for (register label k=i; k<n; k++)
            {
                swap(tmpMatrix[i][k], tmpMatrix[iMax][k]);
            }
            swap(sourceSol[i], sourceSol[iMax]);
        }

        // Check that the system of equations isn't singular
        if (mag(tmpMatrix[i][i]) < 1e-20)
        {
            FatalErrorIn("simpleMatrix<T>::solve()")
                << "Singular Matrix"
                << exit(FatalError);
        }

        // Reduce to upper triangular form
        for (register label j=i+1; j<n; j++)
        {
            sourceSol[j] -= sourceSol[i]*(tmpMatrix[j][i]/tmpMatrix[i][i]);

            for (register label k=n-1; k>=i; k--)
            {
                tmpMatrix[j][k] -=
                    tmpMatrix[i][k]*tmpMatrix[j][i]/tmpMatrix[i][i];
            }
        }
    }

    // Back-substitution
    for (register label j=n-1; j>=0; j--)
    {
        T ntempvec = pTraits<T>::zero;

        for (register label k=j+1; k<n; k++)
        {
            ntempvec += tmpMatrix[j][k]*sourceSol[k];
        }

        sourceSol[j] = (sourceSol[j] - ntempvec)/tmpMatrix[j][j];
    }
}


template<class T>
Field<T> simpleMatrix<T>::solve() const
{
    Matrix<scalar> tmpMatrix = matrix_;
    Field<T> sourceSol = source_;

    solve(tmpMatrix, sourceSol);

    return sourceSol;
}


template<class T>
void simpleMatrix<T>::LUDecompose
(
    Matrix<scalar>& matrix,
    labelList& pivotIndices
)
{
    label n = matrix.n();
    scalarField vv(n);

    for (register label i=0; i<n; i++)
    {
        scalar largestCoeff = 0.0;
        scalar temp;
        for (register label j=0; j<n; j++)
        {
            if ((temp = mag(matrix[i][j])) > largestCoeff)
            {
                largestCoeff = temp;
            }
        }

        if (largestCoeff == 0.0)
        {
            FatalErrorIn
            (
                "simpleMatrix<T>::LUdecompose"
                "(Matrix<scalar>& matrix, labelList& rowIndices)"
            )   << "Singular matrix" << exit(FatalError);
        }

        vv[i] = 1.0/largestCoeff;
    }

    for (register label j=0; j<n; j++)
    {
        for (register label i=0; i<j; i++)
        {
            scalar sum = matrix[i][j];
            for (register label k=0; k<i; k++)
            {
                sum -= matrix[i][k]*matrix[k][j];
            }
            matrix[i][j] = sum;
        }

        label iMax = 0;

        scalar largestCoeff = 0.0;
        for (register label i=j; i<n; i++)
        {
            scalar sum = matrix[i][j];

            for (register label k=0; k<j; k++)
            {
                sum -= matrix[i][k]*matrix[k][j];
            }

            matrix[i][j] = sum;

            scalar temp;
            if ((temp = vv[i]*mag(sum)) >= largestCoeff)
            {
                largestCoeff = temp;
                iMax = i;
            }
        }

        pivotIndices[j] = iMax;

        if (j != iMax)
        {
            //Info<< "Pivoted on " << j << " " << iMax << endl;

            for (register label k=0; k<n; k++)
            {
                swap(matrix[j][k], matrix[iMax][k]);
            }
            
            vv[iMax] = vv[j];
        }

        if (matrix[j][j] == 0.0)
        {
            matrix[j][j] = SMALL;
        }

        if (j != n-1)
        {
            scalar rDiag = 1.0/(matrix[j][j]);
            for (register label i=j+1; i<n; i++)
            {
                matrix[i][j] *= rDiag;
            }
        }
    }
}


template<class T>
void simpleMatrix<T>::LUBacksubstitute
(
    const Matrix<scalar>& luMatrix,
    const labelList& pivotIndices,
    Field<T>& sourceSol
)
{
    label n = luMatrix.n();

    label ii = 0;

    for (register label i=0; i<n; i++)
    {
        label ip = pivotIndices[i];
        T sum = sourceSol[ip];
        sourceSol[ip] = sourceSol[i];

        if (ii != 0)
        {
            for (label j=ii-1; j<i; j++)
            {
                sum -= luMatrix[i][j]*sourceSol[j];
            }
        }
        else if (sum != pTraits<T>::zero)
        {
            ii = i+1;
        }

        sourceSol[i] = sum;
    }

    for (register label i=n-1; i>=0; i--)
    {
        T sum = sourceSol[i];

        for (register label j=i+1; j<n; j++)
        { 
            sum -= luMatrix[i][j]*sourceSol[j];
        }

        sourceSol[i] = sum/luMatrix[i][i];
    }
}


template<class T>
void simpleMatrix<T>::LUsolve
(
    Matrix<scalar>& matrix,
    Field<T>& sourceSol
)
{
    labelList pivotIndices(matrix.n());
    LUDecompose(matrix, pivotIndices);
    LUBacksubstitute(matrix, pivotIndices, sourceSol);
}


template<class T>
Field<T> simpleMatrix<T>::LUsolve() const
{
    Matrix<scalar> luMatrix = matrix_;
    Field<T> sourceSol = source_;

    LUsolve(luMatrix, sourceSol);

    return sourceSol;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void simpleMatrix<T>::operator=(const simpleMatrix<T>& m)
{
    if (this == &m)
    {
        FatalErrorIn("simpleMatrix<T>::operator=(const simpleMatrix<T>&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (matrix_.n() != m.matrix_.n())
    {
        FatalErrorIn("simpleMatrix<T>::operator=(const simpleMatrix<T>&)")
            << "Different size matrices"
            << abort(FatalError);
    }

    if (source_.size() != m.source_.size())
    {
        FatalErrorIn("simpleMatrix<T>::operator=(const simpleMatrix<T>&)")
            << "Different size source vectors"
            << abort(FatalError);
    }

    matrix_ = m.matrix_;
    source_ = m.source_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T>
simpleMatrix<T> operator+(const simpleMatrix<T>& m1, const simpleMatrix<T>& m2)
{
    return simpleMatrix<T>(m1.matrix_ + m2.matrix, m1.source_ + m2.source_);
}


template<class T>
simpleMatrix<T> operator-(const simpleMatrix<T>& m1, const simpleMatrix<T>& m2)
{
    return simpleMatrix<T>(m1.matrix_ - m2.matrix, m1.source_ - m2.source_);
}


template<class T>
simpleMatrix<T> operator*(const scalar s, const simpleMatrix<T>& m)
{
    return simpleMatrix<T>(s*m.matrix_, s*m.source_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Ostream& operator<<(Ostream& os, const simpleMatrix<T>& m)
{
    os << m.matrix_ << nl << m.source_;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
