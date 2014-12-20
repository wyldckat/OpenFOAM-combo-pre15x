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
    Specialisation of Field<T> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarField.H"
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<scalarField> scalarField::component(const direction) const
{
    return *this;
}

void component(scalarField& sf, const UList<scalar>& f, const direction)
{
    sf = f;
}

template<>
void scalarField::replace(const direction, const UList<scalar>& sf)
{
    *this = sf;
}


void stabilise(scalarField& Res, const UList<scalar>& sf, const scalar s)
{
    TFOR_ALL_F_OP_FUNC_S_F
        (scalar, Res, =, ::Foam::stabilise, scalar, s, scalar, sf)
}

tmp<scalarField> stabilise(const UList<scalar>& sf, const scalar s)
{
    tmp<scalarField> Res(new scalarField(sf.size()));
    TFOR_ALL_F_OP_FUNC_S_F
        (scalar, Res(), =, ::Foam::stabilise, scalar, s, scalar, sf)
    return Res;
}

tmp<scalarField> stabilise(const tmp<scalarField>& sf, const scalar s)
{
    tmp<scalarField> Res(sf.ptr());
    TFOR_ALL_F_OP_FUNC_S_F
        (scalar, Res(), =, ::Foam::stabilise, scalar, s, scalar, Res())
    return Res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR_RR(scalar, scalar, scalar, +, add)
BINARY_TYPE_OPERATOR_RR(scalar, scalar, scalar, -, subtract)

BINARY_OPERATOR_RR(scalar, scalar, scalar, *, multiply)
BINARY_OPERATOR_RR(scalar, scalar, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, divide)
BINARY_TYPE_OPERATOR_SR(scalar, scalar, scalar, /, divide)

BINARY_FUNCTION_RR(scalar, scalar, scalar, pow)
BINARY_TYPE_FUNCTION_RR(scalar, scalar, scalar, pow)

BINARY_FUNCTION_RR(scalar, scalar, scalar, atan2)
BINARY_TYPE_FUNCTION_RR(scalar, scalar, scalar, atan2)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION_R(scalar, scalar, pow3)
UNARY_FUNCTION_R(scalar, scalar, pow4)
UNARY_FUNCTION_R(scalar, scalar, sqrt)
UNARY_FUNCTION_R(scalar, scalar, sign)
UNARY_FUNCTION_R(scalar, scalar, pos)
UNARY_FUNCTION_R(scalar, scalar, neg)
UNARY_FUNCTION_R(scalar, scalar, exp)
UNARY_FUNCTION_R(scalar, scalar, log)
UNARY_FUNCTION_R(scalar, scalar, log10)
UNARY_FUNCTION_R(scalar, scalar, sin)
UNARY_FUNCTION_R(scalar, scalar, cos)
UNARY_FUNCTION_R(scalar, scalar, tan)
UNARY_FUNCTION_R(scalar, scalar, asin)
UNARY_FUNCTION_R(scalar, scalar, acos)
UNARY_FUNCTION_R(scalar, scalar, atan)
UNARY_FUNCTION_R(scalar, scalar, sinh)
UNARY_FUNCTION_R(scalar, scalar, cosh)
UNARY_FUNCTION_R(scalar, scalar, tanh)
UNARY_FUNCTION_R(scalar, scalar, asinh)
UNARY_FUNCTION_R(scalar, scalar, acosh)
UNARY_FUNCTION_R(scalar, scalar, atanh)
UNARY_FUNCTION_R(scalar, scalar, erf)
UNARY_FUNCTION_R(scalar, scalar, erfc)
UNARY_FUNCTION_R(scalar, scalar, lgamma)
UNARY_FUNCTION_R(scalar, scalar, j0)
UNARY_FUNCTION_R(scalar, scalar, j1)
UNARY_FUNCTION_R(scalar, scalar, y0)
UNARY_FUNCTION_R(scalar, scalar, y1)


#define BesselFunc(func)                                                      \
void func(scalarField& Res, const int n, const UList<scalar>& sf)             \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Res, =, ::Foam::func, int, n, scalar, sf)  \
}                                                                             \
                                                                              \
tmp<scalarField> func(const int n, const UList<scalar>& sf)                   \
{                                                                             \
    tmp<scalarField> Res(new scalarField(sf.size()));                         \
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Res(), =, ::Foam::func, int, n, scalar, sf)\
    return Res;                                                               \
}                                                                             \
                                                                              \
tmp<scalarField> func(const int n, const tmp<scalarField>& sf)                \
{                                                                             \
    tmp<scalarField> Res(sf.ptr());                                           \
    TFOR_ALL_F_OP_FUNC_S_F                                                    \
        (scalar, Res(), =, ::Foam::func, int, n, scalar, Res())               \
    return Res;                                                               \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
