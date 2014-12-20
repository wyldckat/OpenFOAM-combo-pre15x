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
    Generic field type.

\*---------------------------------------------------------------------------*/

#include "FieldM.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& sf,
    const UList<Type>& f,
    const direction d
)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_F_FUNC_S(cmptType, sf, =,
        Type, f, .component, const direction, d)
}


template<class Type>
void T(Field<Type>& f1, const UList<Type>& f2)
{
    TFOR_ALL_F_OP_F_FUNC(Type, f1, =, Type, f2, T)
}


template<class Type, int r>
void pow
(
    Field<typename powProduct<Type, r>::type>& f,
    const UList<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    TFOR_ALL_F_OP_FUNC_F_S(powProductType, f, =, pow, Type, vf, \
        powProductType, pTraits<powProductType>::zero)
}

template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
    const UList<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType> > result
    (
        new Field<powProductType>(f.size())
    );
    pow<Type, r>(result(), f);
    return result;
}

template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
    const tmp<Field<Type> >& tf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    tmp<Field<powProductType> > result
    (
        new Field<powProductType>(tf().size())
    );
    pow<Type, r>(result(), tf());
    tf.clear();
    return result;
}


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& f,
    const UList<Type>& vf
)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    TFOR_ALL_F_OP_FUNC_F(outerProductType, f, =, sqr, Type, vf)
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const UList<Type>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<Field<outerProductType> > result
    (
        new Field<outerProductType>(f.size())
    );
    sqr(result(), f);
    return result;
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const tmp<Field<Type> >& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    tmp<Field<outerProductType> > result
    (
        new Field<outerProductType>(tf().size())
    );
    sqr(result(), tf());
    tf.clear();
    return result;
}


template<class Type>
void magSqr(Field<scalar>& sf, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, sf, =, magSqr, Type, f)
}

template<class Type>
tmp<Field<scalar> > magSqr(const UList<Type>& f)
{
    tmp<Field<scalar> > result(new Field<scalar>(f.size()));
    magSqr(result(), f);
    return result;
}

template<class Type>
tmp<Field<scalar> > magSqr(const tmp<Field<Type> >& tf)
{
    tmp<Field<scalar> > result(new Field<scalar>(tf().size()));
    magSqr(result(), tf());
    tf.clear();
    return result;
}


template<class Type>
void mag(Field<scalar>& sf, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, sf, =, mag, Type, f)
}

template<class Type>
tmp<Field<scalar> > mag(const UList<Type>& f)
{
    tmp<Field<scalar> > result(new Field<scalar>(f.size()));
    mag(result(), f);
    return result;
}

template<class Type>
tmp<Field<scalar> > mag(const tmp<Field<Type> >& tf)
{
    tmp<Field<scalar> > result(new Field<scalar>(tf().size()));
    mag(result(), tf());
    tf.clear();
    return result;
}


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& cf, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, cf, =, cmptAv, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > result(new Field<cmptType>(f.size()));
    cmptAv(result(), f);
    return result;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const tmp<Field<Type> >& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    tmp<Field<cmptType> > result(new Field<cmptType>(tf().size()));
    cmptAv(result(), tf());
    tf.clear();
    return result;
}


template<class Type>
void cmptMag(Field<Type>& cf, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(Type, cf, =, cmptMag, Type, f)
}

template<class Type>
tmp<Field<Type> > cmptMag(const UList<Type>& f)
{
    tmp<Field<Type> > result(new Field<Type>(f.size()));
    cmptMag(result(), f);
    return result;
}

template<class Type>
tmp<Field<Type> > cmptMag(const tmp<Field<Type> >& tf)
{
    tmp<Field<Type> > result(new Field<Type>(tf().size()));
    cmptMag(result(), tf());
    tf.clear();
    return result;
}


#define BINARY_FUNCTION(func)                                                 \
                                                                              \
template<class Type>                                                          \
void func(Field<Type>& f, const UList<Type>& f1, const UList<Type>& f2)       \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F_F(Type, f, =, func, Type, f1, Type, f2)              \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > func(const UList<Type>& f1, const UList<Type>& f2)          \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    ::Foam::func(tf(), f1, f2);                                               \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > func(const UList<Type>& f1, const tmp<Field<Type> >& tf2)   \
{                                                                             \
    tmp<Field<Type> > tf(tf2.ptr());                                          \
    TFOR_ALL_F_OP_FUNC_F_F(Type, tf(), =, func, Type, f1, Type, tf())         \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > func(const tmp<Field<Type> >& tf1, const UList<Type>& f2)   \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    TFOR_ALL_F_OP_FUNC_F_F(Type, tf(), =, func, Type, tf(), Type, f2)         \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > func                                                        \
(                                                                             \
    const tmp<Field<Type> >& tf1,                                             \
    const tmp<Field<Type> >& tf2                                              \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf2.ptr());                                          \
    TFOR_ALL_F_OP_FUNC_F_F(Type, tf(), =, func, Type, tf1(), Type, tf())      \
    tf1.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
void func(Field<Type>& f, const UList<Type>& f1, const Type& s)               \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F_S(Type, f, =, func, Type, f1, Type, s)               \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > func(const UList<Type>& f1, const Type& s)                  \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    ::Foam::func(tf(), f1, s);                                                \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > func(const tmp<Field<Type> >& tf1, const Type& s)           \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    TFOR_ALL_F_OP_FUNC_F_S(Type, tf(), =, func, Type, tf(), Type, s)          \
    return tf;                                                                \
}

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(scale)

#undef BINARY_FUNCTION


#define TMP_UNARY_FUNCTION(returnType, func)                                  \
                                                                              \
template<class Type>                                                          \
returnType func(const tmp<Field<Type> >& tf1)                                 \
{                                                                             \
    returnType res = func(tf1());                                             \
    tf1.clear();                                                              \
    return res;                                                               \
}

template<class Type>
Type max(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Max, =, max, Type, f, Type, Max)

        return Max;
    }
    else
    {
        WarningIn("max(const UList<Type>&)")
            << "empty field, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Min, =, min, Type, f, Type, Min)

        return Min;
    }
    else
    {
        WarningIn("min(const UList<Type>&)")
            << "empty field, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f)
{
    if (f.size())
    {
        Type Sum = pTraits<Type>::zero;
        TFOR_ALL_S_OP_F(Type, Sum, +=, Type, f)

        return Sum;
    }
    else
    {
        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, sum)


template<class Type>
scalar sumProd(const UList<Type>& f1, const UList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        scalar SumProd = 0.0;
        TFOR_ALL_S_OP_F_OP_F(scalar, SumProd, +=, Type, f1, *, Type, f2)

        return SumProd;
    }
    else
    {
        return 0.0;
    }
}


template<class Type>
scalar sumSqr(const UList<Type>& f)
{
    if (f.size())
    {
        scalar SumSqr = 0.0;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumSqr, +=, sqr, Type, f)

        return SumSqr;
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumSqr)

template<class Type>
scalar sumMag(const UList<Type>& f)
{
    if (f.size())
    {
        scalar SumMag = 0.0;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, mag, Type, f)

        return SumMag;
    }
    else
    {
        return 0.0;
    }
}

TMP_UNARY_FUNCTION(scalar, sumMag)

template<class Type>
Type average(const UList<Type>& f)
{
    if (f.size())
    {
        Type avrg = sum(f)/f.size();

        return avrg;
    }
    else
    {
        WarningIn("average(const UList<Type>&)")
            << "empty field, returning zero" << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, average)


#define G_UNARY_FUNCTION(returnType, gFunc, func, rFunc)                      \
                                                                              \
template<class Type>                                                          \
returnType gFunc(const UList<Type>& f)                                        \
{                                                                             \
    returnType res = func(f);                                                 \
    reduce(res, rFunc##Op<Type>());                                           \
    return res;                                                               \
}                                                                             \
TMP_UNARY_FUNCTION(returnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION

template<class Type>
scalar gSumProd(const UList<Type>& f1, const UList<Type>& f2)
{
    Type SumProd = sumProd(f1, f2);
    reduce(SumProd, sumOp<Type>());
    return SumProd;
}


template<class Type>
Type gAverage(const UList<Type>& f)
{
    label n = f.size();
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        Type avrg = gSum(f)/n;

        return avrg;
    }
    else
    {
        WarningIn("gAverage(const UList<Type>&)")
            << "empty field, returning zero." << endl;

        return pTraits<Type>::zero;
    }
}

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(op, opFunc)                                            \
                                                                              \
template<class Type>                                                          \
void opFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const UList<Type>& f1                                                     \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_OP_F(Type, f, =, op, Type, f1)                              \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const UList<Type>& f1                                                     \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), f1);                                                         \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type> >& tf1                                              \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    opFunc(tf(), tf());                                                       \
    return tf;                                                                \
}

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<class Type>                                                          \
void opFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_F_OP_F(Type, f, =, Type1, f1, op, Type2, f2)                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf2.ptr());                                          \
    opFunc(tf(), f1, tf());                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FT(Type1, Type2, op, opFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf = f1 op tf2();                                       \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    opFunc(tf(), tf(), f2);                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TF(Type1, Type2, op, opFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf = tf1() op f2;                                       \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    opFunc(tf(), tf(), tf2());                                                \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TTR(Type1, Type2, op, opFunc)                         \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf2.ptr());                                          \
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

BINARY_OPERATOR_FF(scalar, Type, *, multiply)
BINARY_OPERATOR_FTR(scalar, Type, *, multiply)
BINARY_OPERATOR_TF(scalar, Type, *, multiply)
BINARY_OPERATOR_TTR(scalar, Type, *, multiply)

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
template<class Type>                                                          \
void opFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const TYPE& s,                                                            \
    const UList<Type>& f1                                                     \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_S_OP_F(Type, f, =, TYPE, s, op, Type, f1)                   \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const TYPE& s,                                                            \
    const UList<Type>& f1                                                     \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), s, f1);                                                      \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const TYPE& s,                                                            \
    const tmp<Field<Type> >& tf1                                              \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
    opFunc(tf(), s, tf());                                                    \
    return tf;                                                                \
}

#define BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)                             \
                                                                              \
template<class Type>                                                          \
void opFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const UList<Type>& f1,                                                    \
    const TYPE& s                                                             \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_F_OP_S(Type, f, =, Type, f1, op, TYPE, s)                   \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const UList<Type>& f1,                                                    \
    const TYPE& s                                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(new Field<Type>(f1.size()));                         \
    opFunc(tf(), f1, s);                                                      \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator op                                                 \
(                                                                             \
    const tmp<Field<Type> >& tf1,                                             \
    const TYPE& s                                                             \
)                                                                             \
{                                                                             \
    tmp<Field<Type> > tf(tf1.ptr());                                          \
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
template<class Type1, class Type2>                                            \
void opFunc                                                                   \
(                                                                             \
    Field<typename product<Type1, Type2>::type>& f,                           \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    TFOR_ALL_F_OP_F_OP_F(productType, f, =, Type1, f1, op, Type2, f2)         \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator op(const UList<Type1>& f1, const UList<Type2>& f2)                   \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tf(new Field<productType>(f1.size()));           \
    opFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator op(const UList<Type1>& f1, const tmp<Field<Type2> >& tf2)            \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tf(new Field<productType>(f1.size()));           \
    opFunc(tf(), f1, tf2());                                                  \
    tf2.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator op(const tmp<Field<Type1> >& tf1, const UList<Type2>& f2)            \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tf(new Field<productType>(tf1().size()));        \
    opFunc(tf(), tf1(), f2);                                                  \
    tf1.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator op(const tmp<Field<Type1> >& tf1, const tmp<Field<Type2> >& tf2)     \
{                                                                             \
    typedef typename product<Type1, Type2>::type productType;                 \
    tmp<Field<productType> > tf(new Field<productType>(tf1().size()));        \
    opFunc(tf(), tf1(), tf2());                                               \
    tf1.clear();                                                              \
    tf2.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
void opFunc                                                                   \
(                                                                             \
    Field<typename product<Type, Form>::type>& f,                             \
    const UList<Type>& f1,                                                    \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
)                                                                             \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    TFOR_ALL_F_OP_F_OP_S(productType, f, =,Type, f1, op, Form,static_cast<const Form&>(vs)) \
}                                                                             \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)    \
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    tmp<Field<productType> > tf(new Field<productType>(f1.size()));           \
    opFunc(tf(), f1, static_cast<const Form&>(vs));                           \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator op(const tmp<Field<Type> >&tf1,const VectorSpace<Form,Cmpt,nCmpt>&vs)\
{                                                                             \
    typedef typename product<Type, Form>::type productType;                   \
    tmp<Field<productType> > tf(new Field<productType>(tf1().size()));        \
    opFunc(tf(), tf1(), static_cast<const Form&>(vs));                        \
    tf1.clear();                                                              \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
void opFunc                                                                   \
(                                                                             \
    Field<typename product<Form, Type>::type>& f,                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const UList<Type>& f1                                                     \
)                                                                             \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    TFOR_ALL_F_OP_S_OP_F(productType, f, =,Form,static_cast<const Form&>(vs), op, Type, f1) \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1)    \
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<Field<productType> > tf(new Field<productType>(f1.size()));           \
    opFunc(tf(), static_cast<const Form&>(vs), f1);                           \
    return tf;                                                                \
}                                                                             \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator op(const VectorSpace<Form,Cmpt,nCmpt>&vs,const tmp<Field<Type> >&tf1)\
{                                                                             \
    typedef typename product<Form, Type>::type productType;                   \
    tmp<Field<productType> > tf(new Field<productType>(tf1().size()));        \
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
