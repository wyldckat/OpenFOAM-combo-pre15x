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
    Specialisation of Field<T> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarField.H"
#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<scalarField> scalarField::component(const direction) const
{
    return *this;
}

template<>
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
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Res, =,
        ::Foam::stabilise, scalar, s, scalar, sf)
}

tmp<scalarField> stabilise(const UList<scalar>& sf, const scalar s)
{
    tmp<scalarField> Res(new scalarField(sf.size()));
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Res(), =,
        ::Foam::stabilise, scalar, s, scalar, sf)
    return Res;
}

tmp<scalarField> stabilise(const tmp<scalarField>& sf, const scalar s)
{
    tmp<scalarField> Res(sf.ptr());
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Res(), =,
        ::Foam::stabilise, scalar, s, scalar, Res())
    return Res;
}


void divide(scalarField& f, const scalar s, const UList<scalar>& f1)
{
    TFOR_ALL_F_OP_S_OP_F(scalar, f, =, scalar, s, /, scalar, f1)
}

tmp<scalarField> operator/(const scalar s, const UList<scalar>& f1)
{
    tmp<scalarField> tf(new scalarField(f1.size()));
    TFOR_ALL_F_OP_S_OP_F(scalar, tf(), =, scalar, s, /, scalar, f1)
    return tf;
}

tmp<scalarField> operator/(const scalar s, const tmp<scalarField >& tf1)
{
    tmp<scalarField> tf(tf1.ptr());
    TFOR_ALL_F_OP_S_OP_F(scalar, tf(), =, scalar, s, /, scalar, tf())
    return tf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global functions, result returned in first argument

void pow(scalarField& Pow, const UList<scalar>& sf1, const UList<scalar>& sf2)
{
    TFOR_ALL_F_OP_FUNC_F_F(scalar, Pow, =,::Foam::pow, scalar, sf1, scalar, sf2)
}

tmp<scalarField> pow(const UList<scalar>& sf1, const UList<scalar>& sf2)
{
    tmp<scalarField> Pow(new scalarField(sf1.size()));
    pow(Pow(), sf1, sf2);
    return Pow;
}

tmp<scalarField> pow(const UList<scalar>& sf1, const tmp<scalarField>& sf2)
{
    tmp<scalarField> Pow(sf2.ptr());
    pow(Pow(), sf1, Pow());
    return Pow;
}

tmp<scalarField> pow(const tmp<scalarField>& sf1, const UList<scalar>& sf2)
{
    tmp<scalarField> Pow(sf1.ptr());
    pow(Pow(), Pow(), sf2);
    return Pow;
}

tmp<scalarField> pow(const tmp<scalarField>& sf1, const tmp<scalarField>& sf2)
{
    tmp<scalarField> Pow(sf1.ptr());
    pow(Pow(), Pow(), sf2());
    sf2.clear();
    return Pow;
}


void pow(scalarField& Pow, const UList<scalar>& sf, const scalar& s)
{
    TFOR_ALL_F_OP_FUNC_F_S(scalar, Pow, =, ::Foam::pow, scalar, sf, scalar, s)
}

tmp<scalarField> pow(const UList<scalar>& sf, const scalar& s)
{
    tmp<scalarField> Pow(new scalarField(sf.size()));
    pow(Pow(), sf, s);
    return Pow;
}

tmp<scalarField> pow(const tmp<scalarField>& sf, const scalar& s)
{
    tmp<scalarField> Pow(sf.ptr());
    pow(Pow(), Pow(), s);
    return Pow;
}


void pow(scalarField& Pow, const scalar& s, const UList<scalar>& sf)
{
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Pow, =, ::Foam::pow, scalar, s, scalar, sf)
}

tmp<scalarField> pow(const scalar& s, const UList<scalar>& sf)
{
    tmp<scalarField> Pow(new scalarField(sf.size()));
    pow(Pow(), s, sf);
    return Pow;
}

tmp<scalarField> pow(const scalar& s, const tmp<scalarField>& sf)
{
    tmp<scalarField> Pow(sf.ptr());
    pow(Pow(), s, Pow());
    return Pow;
}


void atan2(scalarField& Atan2, const UList<scalar>& sf1, const UList<scalar>& sf2)
{
    TFOR_ALL_F_OP_FUNC_F_F(scalar, Atan2, =,::Foam::atan2, scalar, sf1, scalar, sf2)
}

tmp<scalarField> atan2(const UList<scalar>& sf1, const UList<scalar>& sf2)
{
    tmp<scalarField> Atan2(new scalarField(sf1.size()));
    atan2(Atan2(), sf1, sf2);
    return Atan2;
}

tmp<scalarField> atan2(const UList<scalar>& sf1, const tmp<scalarField>& sf2)
{
    tmp<scalarField> Atan2(sf2.ptr());
    atan2(Atan2(), sf1, Atan2());
    return Atan2;
}

tmp<scalarField> atan2(const tmp<scalarField>& sf1, const UList<scalar>& sf2)
{
    tmp<scalarField> Atan2(sf1.ptr());
    atan2(Atan2(), Atan2(), sf2);
    return Atan2;
}

tmp<scalarField> atan2(const tmp<scalarField>& sf1, const tmp<scalarField>& sf2)
{
    tmp<scalarField> Atan2(sf1.ptr());
    atan2(Atan2(), Atan2(), sf2());
    sf2.clear();
    return Atan2;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define transFunc(func)                                                       \
void func(scalarField& Res, const UList<scalar>& sf)                          \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F(scalar, Res, =, ::Foam::func, scalar, sf)            \
}                                                                             \
                                                                              \
tmp<scalarField> func(const UList<scalar>& sf)                                \
{                                                                             \
    tmp<scalarField> Res(new scalarField(sf.size()));                         \
    TFOR_ALL_F_OP_FUNC_F(scalar, Res(), =, ::Foam::func, scalar, sf)          \
    return Res;                                                               \
}                                                                             \
                                                                              \
tmp<scalarField> func(const tmp<scalarField>& sf)                             \
{                                                                             \
    tmp<scalarField> Res(sf.ptr());                                           \
    TFOR_ALL_F_OP_FUNC_F(scalar, Res(), =, ::Foam::func, scalar, Res())       \
    return Res;                                                               \
}

transFunc(pow3)
transFunc(pow4)
transFunc(sqrt)
transFunc(sign)
transFunc(pos)
transFunc(neg)
transFunc(exp)
transFunc(log)
transFunc(log10)
transFunc(sin)
transFunc(cos)
transFunc(tan)
transFunc(asin)
transFunc(acos)
transFunc(atan)
transFunc(sinh)
transFunc(cosh)
transFunc(tanh)
transFunc(asinh)
transFunc(acosh)
transFunc(atanh)
transFunc(erf)
transFunc(erfc)
transFunc(lgamma)
transFunc(j0)
transFunc(j1)
transFunc(y0)
transFunc(y1)

#undef transFunc


#define transFunc(func)                                                       \
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
    TFOR_ALL_F_OP_FUNC_S_F(scalar, Res(), =,                                  \
        ::Foam::func, int, n, scalar, Res())                                  \
    return Res;                                                               \
}

transFunc(jn)
transFunc(yn)

#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
