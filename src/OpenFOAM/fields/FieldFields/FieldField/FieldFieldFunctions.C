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
    Generic FieldField type.

\*---------------------------------------------------------------------------*/

#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<template<class> class Field, class Type>
void component
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& sf,
    const FieldField<Field, Type>& f,
    const direction d
)
{
    forAll(sf, i)
    {
        component(sf[i], f[i], d);
    }
}


template<template<class> class Field, class Type>
void T(FieldField<Field, Type>& f1, const FieldField<Field, Type>& f2)
{
    forAll(f1, i)
    {
        T(f1[i], f2[i]);
    }
}


template<template<class> class Field, class Type, int r>
void pow
(
    FieldField<Field, typename powProduct<Type, r>::type>& f,
    const FieldField<Field, Type>& vf
)
{
    forAll(f, i)
    {
        pow(f[i], vf[i]);
    }
}

template<template<class> class Field, class Type, int r>
tmp<FieldField<Field, typename powProduct<Type, r>::type> >
pow
(
    const FieldField<Field, Type>& f, typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<FieldField<Field, powProductType> > result
    (
        FieldField<Field, powProductType>::NewCalculatedType(f)
    );
    pow<Type, r>(result(), f);
    return result;
}

template<template<class> class Field, class Type, int r>
tmp<FieldField<Field, typename powProduct<Type, r>::type> >
pow
(
    const tmp<FieldField<Field, Type> >& tf, typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<FieldField<Field, powProductType> > result
    (
        FieldField<Field, powProductType>::NewCalculatedType(tf())
    );
    pow<Type, r>(result(), tf());
    tf.clear();
    return result;
}


template<template<class> class Field, class Type>
void sqr
(
    FieldField<Field, typename outerProduct<Type, Type>::type>& f,
    const FieldField<Field, Type>& vf
)
{
    forAll(f, i)
    {
        sqr(f[i], vf[i]);
    }
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename outerProduct<Type, Type>::type> >
sqr(const FieldField<Field, Type>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<FieldField<Field, outerProductType> > result
    (
        FieldField<Field, outerProductType>::NewCalculatedType(f)
    );
    sqr(result(), f);
    return result;
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename outerProduct<Type, Type>::type> >
sqr(const tmp<FieldField<Field, Type> >& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<FieldField<Field, outerProductType> > result
    (
        FieldField<Field, outerProductType>::NewCalculatedType(tf())
    );
    sqr(result(), tf());
    tf.clear();
    return result;
}


template<template<class> class Field, class Type>
void magSqr(FieldField<Field, scalar>& sf, const FieldField<Field, Type>& f)
{
    forAll(sf, i)
    {
        magSqr(sf[i], f[i]);
    }
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > magSqr(const FieldField<Field, Type>& f)
{
    tmp<FieldField<Field, scalar> > result
    (
        FieldField<Field, scalar>::NewCalculatedType(f)
    );

    magSqr(result(), f);
    return result;
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > magSqr(const tmp<FieldField<Field, Type> >& tf)
{
    tmp<FieldField<Field, scalar> > result
    (
        FieldField<Field, scalar>::NewCalculatedType(tf())
    );

    magSqr(result(), tf());
    tf.clear();
    return result;
}


template<template<class> class Field, class Type>
void mag(FieldField<Field, scalar>& sf, const FieldField<Field, Type>& f)
{
    forAll(sf, i)
    {
        mag(sf[i], f[i]);
    }
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > mag(const FieldField<Field, Type>& f)
{
    tmp<FieldField<Field, scalar> > result
    (
        FieldField<Field, scalar>::NewCalculatedType(f)
    );

    mag(result(), f);
    return result;
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > mag(const tmp<FieldField<Field, Type> >& tf)
{
    tmp<FieldField<Field, scalar> > result
    (
        FieldField<Field, scalar>::NewCalculatedType(tf())
    );

    mag(result(), tf());
    tf.clear();
    return result;
}


template<template<class> class Field, class Type>
void cmptAv
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& cf,
    const FieldField<Field, Type>& f
)
{
    forAll(cf, i)
    {
        cmptAv(cf[i], f[i]);
    }
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptAv
(
    const FieldField<Field, Type>& f
)
{
    typedef typename FieldField<Field, Type>::cmptType cmptType;
    tmp<FieldField<Field, cmptType> > result
    (
        FieldField<Field, cmptType>::NewCalculatedType(f)
    );
    cmptAv(result(), f);
    return result;
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptAv
(
    const tmp<FieldField<Field, Type> >& tf
)
{
    typedef typename FieldField<Field, Type>::cmptType cmptType;
    tmp<FieldField<Field, cmptType> > result
    (
        FieldField<Field, cmptType>::NewCalculatedType(tf())
    );
    cmptAv(result(), tf());
    tf.clear();
    return result;
}


template<template<class> class Field, class Type>
void cmptMag
(
    FieldField<Field, Type>& cf,
    const FieldField<Field, Type>& f
)
{
    forAll(cf, i)
    {
        cmptMag(cf[i], f[i]);
    }
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type> > cmptMag
(
    const FieldField<Field, Type>& f
)
{
    tmp<FieldField<Field, Type> > result
    (
        FieldField<Field, Type>::NewCalculatedType(f)
    );
    cmptMag(result(), f);
    return result;
}

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type> > cmptMag
(
    const tmp<FieldField<Field, Type> >& tf
)
{
    tmp<FieldField<Field, Type> > result
    (
        FieldField<Field, Type>::NewCalculatedType(tf())
    );
    cmptMag(result(), tf());
    tf.clear();
    return result;
}


#define BINARY_FUNCTION(func)                                                 \
                                                                              \
template<template<class> class Field, class Type>                             \
void func                                                                     \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1,                                        \
    const FieldField<Field, Type>& f2                                         \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        func(f[i], f1[i], f2[i]);                                             \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const FieldField<Field, Type>& f2                                         \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf                                          \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f1)                        \
    );                                                                        \
    func(tf(), f1, f2);                                                       \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const tmp<FieldField<Field, Type> >& tf2                                  \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf2.ptr());                              \
    FieldField<Field, Type>& f = tf();                                        \
    func(f, f1, f);                                                           \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const FieldField<Field, Type>& f2                                         \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    FieldField<Field, Type>& f = tf();                                        \
    func(f, f, f2);                                                           \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const tmp<FieldField<Field, Type> >& tf2                                  \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf2.ptr());                              \
    FieldField<Field, Type>& f = tf();                                        \
    func(f, tf1(), f);                                                        \
    tf1.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
void func                                                                     \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1,                                        \
    const Type& s                                                             \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        func(f[i], f1[i], s);                                                 \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const Type& s                                                             \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf                                          \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f1)                        \
    );                                                                        \
    func(tf(), f1, s);                                                        \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const Type& s                                                             \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    FieldField<Field, Type>& f = tf();                                        \
    func(f, f, s);                                                            \
    return tf;                                                                \
}

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(scale)

#undef BINARY_FUNCTION


#define TMP_UNARY_FUNCTION(returnType, func)                                  \
                                                                              \
template<template<class> class Field, class Type>                             \
returnType func(const tmp<FieldField<Field, Type> >& tf1)                     \
{                                                                             \
    returnType res = func(tf1());                                             \
    tf1.clear();                                                              \
    return res;                                                               \
}

template<template<class> class Field, class Type>
Type max(const FieldField<Field, Type>& f)
{
    if (f.size())
    {
        label i = 0;
        while(!f[i].size()) i++;

        Type Max(max(f[i]));

        for (label j=i+1; j<f.size(); j++)
        {
            if (f[j].size())
            {
                Max = max(max(f[j]), Max);
            }
        }

        return Max;
    }
    else
    {
        WarningIn("max(const FieldField<Field, Type>&) const")
            << "empty fieldField, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, max)

template<template<class> class Field, class Type>
Type min(const FieldField<Field, Type>& f)
{
    if (f.size())
    {
        label i = 0;
        while(!f[i].size()) i++;

        Type Min(min(f[i]));

        for (label j=i+1; j<f.size(); j++)
        {
            if (f[j].size())
            {
                Min = min(min(f[j]), Min);
            }
        }

        return Min;
    }
    else
    {
        WarningIn("min(const FieldField<Field, Type>&) const")
            << "empty fieldField, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, min)

template<template<class> class Field, class Type>
Type sum(const FieldField<Field, Type>& f)
{
    if (f.size())
    {
        Type Sum = pTraits<Type>::zero;

        forAll(f, i)
        {
            Sum += sum(f[i]);
        }

        return Sum;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, sum)

template<template<class> class Field, class Type>
scalar sumMag(const FieldField<Field, Type>& f)
{
    if (f.size())
    {
        scalar SumMag = 0.0;

        forAll(f, i)
        {
            SumMag += sumMag(f[i]);
        }

        return SumMag;
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumMag)

template<template<class> class Field, class Type>
Type average(const FieldField<Field, Type>& f)
{
    if (f.size())
    {
        label n = 0;

        forAll(f, i)
        {
            n += f[i].size();
        }

        if (n == 0)
        {
            WarningIn("average(const FieldField<Field, Type>&) const")
                << "empty fieldField, returning zero" << endl;

            return pTraits<Type>::zero;
        }

        Type avrg = sum(f)/n;

        return avrg;
    }
    else
    {
        WarningIn("average(const FieldField<Field, Type>&) const")
            << "empty fieldField, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, average)


#include "PstreamReduceOps.H"

#define G_UNARY_FUNCTION(returnType, gFunc, func, rFunc)                      \
                                                                              \
template<template<class> class Field, class Type>                             \
returnType gFunc(const FieldField<Field, Type>& f)                            \
{                                                                             \
    returnType res = func(f);                                                 \
    reduce(res, rFunc##Op<Type>());                                           \
    return res;                                                               \
}                                                                             \
TMP_UNARY_FUNCTION(returnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION


template<template<class> class Field, class Type>
Type gAverage(const FieldField<Field, Type>& f)
{
    label n = 0;

    forAll(f, i)
    {
        n += f[i].size();
    }

    reduce(n, sumOp<label>());

    if (n > 0)
    {
        Type avrg = gSum(f)/n;

        return avrg;
    }
    else
    {
        WarningIn("gAverage(const FieldField<Field, Type>&) const")
            << "empty fieldField, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(op, opFunc)                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1                                         \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], f1[i]);                                                  \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type>& f1                                         \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf                                          \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f1)                        \
    );                                                                        \
    opFunc(tf(), f1);                                                         \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1                                  \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    opFunc(tf(), tf());                                                       \
    return tf;                                                                \
}

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<template<class> class Field, class Type>                             \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], f1[i], f2[i]);                                           \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf                                          \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f1)                        \
    );                                                                        \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf2.ptr());                              \
    opFunc(tf(), f1, tf());                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FT(Type1, Type2, op, opFunc)                          \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf = f1 op tf2();                           \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    opFunc(tf(), tf(), f2);                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TF(Type1, Type2, op, opFunc)                          \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf = tf1() op f2;                           \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    opFunc(tf(), tf(), tf2());                                                \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TTR(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf2.ptr());                              \
    opFunc(tf(), tf1(), tf());                                                \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_R(Type1, Type2, op, opFunc)                           \
    BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                              \
    BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                             \
    BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                             \
    BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)

BINARY_OPERATOR_R(Type, Type, +, add)
BINARY_OPERATOR_R(Type, Type, -, subtract)

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
BINARY_OPERATOR_FF(scalar, Type, *, multiply)
BINARY_OPERATOR_FTR(scalar, Type, *, multiply)
BINARY_OPERATOR_TF(scalar, Type, *, multiply)
BINARY_OPERATOR_TTR(scalar, Type, *, multiply)
#endif

BINARY_OPERATOR_FF(Type, scalar, /, divide)
BINARY_OPERATOR_FT(Type, scalar, /, divide)
BINARY_OPERATOR_TRF(Type, scalar, /, divide)
BINARY_OPERATOR_TRT(Type, scalar, /, divide)

#undef BINARY_OPERATOR_R
#undef BINARY_OPERATOR_FF
#undef BINARY_OPERATOR_FTR
#undef BINARY_OPERATOR_TF
#undef BINARY_OPERATOR_TTR
#undef BINARY_OPERATOR_FT
#undef BINARY_OPERATOR_TRF
#undef BINARY_OPERATOR_TRT


#define BINARY_TYPE_OPERATOR_TF(TYPE, op, opFunc)                             \
                                                                              \
template<template<class> class Field, class Type>                             \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const TYPE& s,                                                            \
    const FieldField<Field, Type>& f1                                         \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], s, f1[i]);                                               \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const TYPE& s,                                                            \
    const FieldField<Field, Type>& f1                                         \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf                                          \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f1)                        \
    );                                                                        \
    opFunc(tf(), s, f1);                                                      \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const TYPE& s,                                                            \
    const tmp<FieldField<Field, Type> >& tf1                                  \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    opFunc(tf(), s, tf());                                                    \
    return tf;                                                                \
}

#define BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)                             \
                                                                              \
template<template<class> class Field, class Type>                             \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1,                                        \
    const TYPE& s                                                             \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], f1[i], s);                                               \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const TYPE& s                                                             \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf                                          \
    (                                                                         \
        FieldField<Field, Type>::NewCalculatedType(f1)                        \
    );                                                                        \
    opFunc(tf(), f1, s);                                                      \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const TYPE& s                                                             \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, Type> > tf(tf1.ptr());                              \
    opFunc(tf(), tf(), s);                                                    \
    return tf;                                                                \
}

#define BINARY_TYPE_OPERATOR(TYPE, op, opFunc)                                \
    BINARY_TYPE_OPERATOR_TF(TYPE, op, opFunc)                                 \
    BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)

BINARY_TYPE_OPERATOR(Type, +, add)
BINARY_TYPE_OPERATOR(Type, -, subtract)

BINARY_TYPE_OPERATOR(scalar, *, multiply)
BINARY_TYPE_OPERATOR_FT(scalar, /, divide)

#undef BINARY_TYPE_OPERATOR
#undef BINARY_TYPE_OPERATOR_TF
#undef BINARY_TYPE_OPERATOR_FT


#define PRODUCT_OPERATOR(product, op, opFunc)                                 \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, typename product<Type1, Type2>::type>& f,               \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], f1[i], f2[i]);                                           \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(f1)                 \
    );                                                                        \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(f1)                 \
    );                                                                        \
    opFunc(tf(), f1, tf2());                                                  \
    tf2.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(tf1())              \
    );                                                                        \
    opFunc(tf(), tf1(), f2);                                                  \
    tf1.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(tf1())              \
    );                                                                        \
    opFunc(tf(), tf1(), tf2());                                               \
    tf1.clear();                                                              \
    tf2.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template                                                                      \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>  \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, typename product<Type, Form>::type>& f,                 \
    const FieldField<Field, Type>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], f1[i], vs);                                              \
    }                                                                         \
}                                                                             \
                                                                              \
template                                                                      \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>  \
tmp<FieldField<Field, typename product<Type, Form>::type> >                   \
operator op                                                                   \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(f1)                 \
    );                                                                        \
    opFunc(tf(), f1, static_cast<const Form&>(vs));                           \
    return tf;                                                                \
}                                                                             \
                                                                              \
template                                                                      \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>  \
tmp<FieldField<Field, typename product<Type, Form>::type> >                   \
operator op                                                                   \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(tf1())              \
    );                                                                        \
    opFunc(tf(), tf1(), static_cast<const Form&>(vs));                        \
    tf1.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template                                                                      \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>  \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, typename product<Form, Type>::type>& f,                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const FieldField<Field, Type>& f1                                         \
)                                                                             \
{                                                                             \
    forAll(f, i)                                                              \
    {                                                                         \
        opFunc(f[i], vs, f1[i]);                                              \
    }                                                                         \
}                                                                             \
                                                                              \
template                                                                      \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>  \
tmp<FieldField<Field, typename product<Form, Type>::type> >                   \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const FieldField<Field, Type>& f1                                         \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(f1)                 \
    );                                                                        \
    opFunc(tf(), static_cast<const Form&>(vs), f1);                           \
    return tf;                                                                \
}                                                                             \
                                                                              \
template                                                                      \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>  \
tmp<FieldField<Field, typename product<Form, Type>::type> >                   \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const tmp<FieldField<Field, Type> >& tf1                                  \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<FieldField<Field, productType> > tf                                   \
    (                                                                         \
        FieldField<Field, productType>::NewCalculatedType(tf1())              \
    );                                                                        \
    opFunc(tf(), static_cast<const Form&>(vs), tf1());                        \
    tf1.clear();                                                              \
    return tf;                                                                \
}

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
