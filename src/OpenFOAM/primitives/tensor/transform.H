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
    transform

Description
    3D tensor transformation operations.

\*---------------------------------------------------------------------------*/

#ifndef transform_H
#define transform_H

#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline tensor transformationTensor
(
    const vector& n1,
    const vector& n2
)
{
    return
        (n1 & n2)*I
      + (1 - (n1 & n2))*sqr(n1 ^ n2)/(magSqr(n1 ^ n2) + VSMALL)
      + (n2*n1 - n1*n2);
}


inline scalar transform(const tensor&, const scalar s)
{
    return s;
}


template<class Cmpt>
inline Vector<Cmpt> transform(const tensor& tt, const Vector<Cmpt>& v)
{
    return tt & v;
}


template<class Cmpt>
inline Tensor<Cmpt> transform(const tensor& tt, const Tensor<Cmpt>& t)
{
    return tt & t & tt.T();
}


template<class Cmpt>
inline SphericalTensor<Cmpt> transform
(
    const tensor& tt,
    const SphericalTensor<Cmpt>& st
)
{
    return st;
}


template<class Type1, class Type2>
inline Type1 transformMask(const Type2& t)
{
    return t;
}


template<>
inline SphericalTensor<scalar> transformMask(const tensor& t)
{
    return SphericalTensor<scalar>::one;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
