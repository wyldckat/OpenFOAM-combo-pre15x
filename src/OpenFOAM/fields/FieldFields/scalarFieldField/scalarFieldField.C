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
    Specialisation of FieldField<T> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<template<class> class Field>
void stabilise
(
    FieldField<Field, scalar>& f,
    const FieldField<Field, scalar>& f1,
    const scalar s
)
{
    forAll(f, i)
    {
        stabilise(f[i], f1[i], s);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > stabilise
(
    const FieldField<Field, scalar>& f1,
    const scalar s
)
{
    tmp<FieldField<Field, scalar> > tf
    (
        FieldField<Field, scalar>::NewCalculatedType(f1)
    );
    stabilise(tf(), f1, s);
    return tf;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > stabilise
(
    const tmp<FieldField<Field, scalar> >& tf1,
    const scalar s
)
{
    tmp<FieldField<Field, scalar> > tf(tf1.ptr());
    stabilise(tf(), tf(), s);
    return tf;
}


template<template<class> class Field>
void divide
(
    FieldField<Field, scalar>& f,
    const scalar s,
    const FieldField<Field, scalar>& f1
)
{
    forAll(f, i)
    {
        divide(f[i], s, f1[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > operator/
(
    const scalar s,
    const FieldField<Field, scalar>& f1
)
{
    tmp<FieldField<Field, scalar> > tf
    (
        FieldField<Field, scalar>::NewCalculatedType(f1)
    );
    divide(tf(), s, f1);
    return tf;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > operator/
(
    const scalar s,
    const tmp<FieldField<Field, scalar> >& tf1
)
{
    tmp<FieldField<Field, scalar> > tf(tf1.ptr());
    divide(tf(), s, tf());
    return tf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global functions, result returned in first argument

template<template<class> class Field>
void pow
(
    FieldField<Field, scalar>& Pow,
    const FieldField<Field, scalar>& sf1,
    const FieldField<Field, scalar>& sf2
)
{
    forAll(Pow, i)
    {
        pow(Pow[i], sf1[i], sf2[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const FieldField<Field, scalar>& sf1,
    const FieldField<Field, scalar>& sf2
)
{
    tmp<FieldField<Field, scalar> > Pow
    (
        FieldField<Field, scalar>::NewCalculatedType(sf1)
    );
    pow(Pow(), sf1, sf2);
    return Pow;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const FieldField<Field, scalar>& sf1,
    const tmp<FieldField<Field, scalar> >& sf2
)
{
    tmp<FieldField<Field, scalar> > Pow(sf2.ptr());
    pow(Pow(), sf1, Pow());
    return Pow;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const tmp<FieldField<Field, scalar> >& sf1,
    const FieldField<Field, scalar>& sf2
)
{
    tmp<FieldField<Field, scalar> > Pow(sf1.ptr());
    pow(Pow(), Pow(), sf2);
    return Pow;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const tmp<FieldField<Field, scalar> >& sf1,
    const tmp<FieldField<Field, scalar> >& sf2
)
{
    tmp<FieldField<Field, scalar> > Pow(sf1.ptr());
    pow(Pow(), Pow(), sf2());
    sf2.clear();
    return Pow;
}


template<template<class> class Field>
void pow
(
    FieldField<Field, scalar>& Pow,
    const FieldField<Field, scalar>& sf,
    const scalar& s
)
{
    forAll(Pow, i)
    {
        pow(Pow[i], sf[i], s);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const FieldField<Field, scalar>& sf,
    const scalar& s
)
{
    tmp<FieldField<Field, scalar> > Pow
    (
        FieldField<Field, scalar>::NewCalculatedType(sf)
    );
    pow(Pow(), sf, s);
    return Pow;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const tmp<FieldField<Field, scalar> >& sf,
    const scalar& s
)
{
    tmp<FieldField<Field, scalar> > Pow(sf.ptr());
    pow(Pow(), Pow(), s);
    return Pow;
}


template<template<class> class Field>
void pow
(
    FieldField<Field, scalar>& Pow,
    const scalar& s,
    const FieldField<Field, scalar>& sf
)
{
    forAll(Pow, i)
    {
        pow(Pow[i], s, sf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const scalar& s,
    const FieldField<Field, scalar>& sf
)
{
    tmp<FieldField<Field, scalar> > Pow
    (
        FieldField<Field, scalar>::NewCalculatedType(sf)
    );
    pow(Pow(), s, sf);
    return Pow;
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const scalar& s,
    const tmp<FieldField<Field, scalar> >& sf
)
{
    tmp<FieldField<Field, scalar> > Pow(sf.ptr());
    pow(Pow(), s, Pow());
    return Pow;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define transFunc(func)                                                       \
                                                                              \
template<template<class> class Field>                                         \
void func(FieldField<Field, scalar>& Res, const FieldField<Field, scalar>& sf)\
{                                                                             \
    forAll(Res, i)                                                            \
    {                                                                         \
        func(Res[i], sf[i]);                                                  \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, scalar> > func(const FieldField<Field, scalar>& sf)     \
{                                                                             \
    tmp<FieldField<Field, scalar> > Res                                       \
    (                                                                         \
        FieldField<Field, scalar>::NewCalculatedType(sf)                      \
    );                                                                        \
    forAll(Res, i)                                                            \
    {                                                                         \
        func(Res[i], sf[i]);                                                  \
    }                                                                         \
    return Res;                                                               \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, scalar> > func                                          \
(                                                                             \
    const tmp<FieldField<Field, scalar> >& sf                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, scalar> > tRes(sf.ptr());                           \
    FieldField<Field, scalar>& Res = tRes();                                  \
    forAll(Res, i)                                                            \
    {                                                                         \
        func(Res[i], Res[i]);                                                 \
    }                                                                         \
    return tRes;                                                              \
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
                                                                              \
template<template<class> class Field>                                         \
void func                                                                     \
(                                                                             \
    FieldField<Field, scalar>& Res,                                           \
    const int n,                                                              \
    const FieldField<Field, scalar>& sf                                       \
)                                                                             \
{                                                                             \
    forAll(Res, i)                                                            \
    {                                                                         \
        func(Res[i], n, sf[i]);                                               \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, scalar> > func                                          \
(                                                                             \
    const int n,                                                              \
    const FieldField<Field, scalar>& sf                                       \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, scalar> > Res                                       \
    (                                                                         \
        FieldField<Field, scalar>::NewCalculatedType(sf)                      \
    );                                                                        \
    forAll(Res, i)                                                            \
    {                                                                         \
        func(Res[i], n, sf[i]);                                               \
    }                                                                         \
    return Res;                                                               \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, scalar> > func                                          \
(                                                                             \
    const int n,                                                              \
    const tmp<FieldField<Field, scalar> >& sf                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, scalar> > tRes(sf.ptr());                           \
    FieldField<Field, scalar>& Res = tRes();                                  \
    forAll(Res, i)                                                            \
    {                                                                         \
        func(Res[i], n, Res[i]);                                              \
    }                                                                         \
    return tRes;                                                              \
}

transFunc(jn)
transFunc(yn)

#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
