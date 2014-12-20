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
);

template<template<class> class Field, class Type>
void T(FieldField<Field, Type>& f1, const FieldField<Field, Type>& f2);


template<template<class> class Field, class Type, int r>
void pow
(
    FieldField<Field, typename powProduct<Type, r>::type>& f,
    const FieldField<Field, Type>& vf
);

template<template<class> class Field, class Type, int r>
tmp<FieldField<Field, typename powProduct<Type, r>::type> >
pow
(
    const FieldField<Field, Type>& f, typename powProduct<Type, r>::type
        = pTraits<typename powProduct<Type, r>::type>::zero
);

template<template<class> class Field, class Type, int r>
tmp<FieldField<Field, typename powProduct<Type, r>::type> >
pow
(
    const tmp<FieldField<Field, Type> >& tf, typename powProduct<Type, r>::type
        = pTraits<typename powProduct<Type, r>::type>::zero
);


template<template<class> class Field, class Type>
void sqr
(
    FieldField<Field, typename outerProduct<Type, Type>::type>& f,
    const FieldField<Field, Type>& vf
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename outerProduct<Type, Type>::type> >
sqr(const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename outerProduct<Type, Type>::type> >
sqr(const tmp<FieldField<Field, Type> >& tf);

template<template<class> class Field, class Type>
void magSqr(FieldField<Field, scalar>& sf, const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > magSqr(const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > magSqr(const tmp<FieldField<Field, Type> >& tf);

template<template<class> class Field, class Type>
void mag(FieldField<Field, scalar>& sf, const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > mag(const FieldField<Field, Type>& f);

template<template<class> class Field, class Type>
tmp<FieldField<Field, scalar> > mag(const tmp<FieldField<Field, Type> >& tf);


template<template<class> class Field, class Type>
void cmptAv
(
    FieldField<Field, typename FieldField<Field, Type>::cmptType>& cf,
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptAv
(
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, typename FieldField<Field, Type>::cmptType> > cmptAv
(
    const tmp<FieldField<Field, Type> >& tf
);

template<template<class> class Field, class Type>
void cmptMag
(
    FieldField<Field, Type>& cf,
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type> > cmptMag
(
    const FieldField<Field, Type>& f
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type> > cmptMag
(
    const tmp<FieldField<Field, Type> >& tf
);

#define BINARY_FUNCTION(func)                                                 \
                                                                              \
template<template<class> class Field, class Type>                             \
void func                                                                     \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1,                                        \
    const FieldField<Field, Type>& f2                                         \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const FieldField<Field, Type>& f2                                         \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const tmp<FieldField<Field, Type> >& tf2                                  \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const FieldField<Field, Type>& f2                                         \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const tmp<FieldField<Field, Type> >& tf2                                  \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
void func                                                                     \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1,                                        \
    const Type& s                                                             \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const Type& s                                                             \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > func                                            \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const Type& s                                                             \
);

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(scale)

#undef BINARY_FUNCTION


#define TMP_UNARY_FUNCTION(returnType, func)                                  \
                                                                              \
template<template<class> class Field, class Type>                             \
returnType func(const tmp<FieldField<Field, Type> >& tf1);

template<template<class> class Field, class Type>
Type max(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, max)

template<template<class> class Field, class Type>
Type min(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, min)

template<template<class> class Field, class Type>
Type sum(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, sum)

template<template<class> class Field, class Type>
scalar sumMag(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(scalar, sumMag)

template<template<class> class Field, class Type>
Type average(const FieldField<Field, Type>& f);

TMP_UNARY_FUNCTION(Type, average)


#include "PstreamReduceOps.H"

#define G_UNARY_FUNCTION(returnType, gFunc, func, rFunc)                      \
                                                                              \
template<template<class> class Field, class Type>                             \
returnType gFunc(const FieldField<Field, Type>& f);                           \
TMP_UNARY_FUNCTION(returnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION


template<template<class> class Field, class Type>
Type gAverage(const FieldField<Field, Type>& f);

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
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type>& f1                                         \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1                                  \
);

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
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
);

#define BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);

#define BINARY_OPERATOR_FT(Type1, Type2, op, opFunc)                          \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);

#define BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
);

#define BINARY_OPERATOR_TF(Type1, Type2, op, opFunc)                          \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
);

#define BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);

#define BINARY_OPERATOR_TTR(Type1, Type2, op, opFunc)                         \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);

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
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const TYPE& s,                                                            \
    const FieldField<Field, Type>& f1                                         \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const TYPE& s,                                                            \
    const tmp<FieldField<Field, Type> >& tf1                                  \
);

#define BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)                             \
                                                                              \
template<template<class> class Field, class Type>                             \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, Type>& f,                                               \
    const FieldField<Field, Type>& f1,                                        \
    const TYPE& s                                                             \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const TYPE& s                                                             \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type>                             \
tmp<FieldField<Field, Type> > operator op                                     \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const TYPE& s                                                             \
);

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
);                                                                            \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const FieldField<Field, Type1>& f1,                                       \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const FieldField<Field, Type2>& f2                                        \
);                                                                            \
                                                                              \
template<template<class> class Field, class Type1, class Type2>               \
tmp<FieldField<Field, typename product<Type1, Type2>::type> >                 \
operator op                                                                   \
(                                                                             \
    const tmp<FieldField<Field, Type1> >& tf1,                                \
    const tmp<FieldField<Field, Type2> >& tf2                                 \
);                                                                            \
                                                                              \
template                                                                      \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>  \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, typename product<Type, Form>::type>& f,                 \
    const FieldField<Field, Type>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template                                                                      \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>  \
tmp<FieldField<Field, typename product<Type, Form>::type> >                   \
operator op                                                                   \
(                                                                             \
    const FieldField<Field, Type>& f1,                                        \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template                                                                      \
<template<class> class Field, class Type, class Form, class Cmpt, int nCmpt>  \
tmp<FieldField<Field, typename product<Type, Form>::type> >                   \
operator op                                                                   \
(                                                                             \
    const tmp<FieldField<Field, Type> >& tf1,                                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template                                                                      \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>  \
void opFunc                                                                   \
(                                                                             \
    FieldField<Field, typename product<Form, Type>::type>& f,                 \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const FieldField<Field, Type>& f1                                         \
);                                                                            \
                                                                              \
template                                                                      \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>  \
tmp<FieldField<Field, typename product<Form, Type>::type> >                   \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const FieldField<Field, Type>& f1                                         \
);                                                                            \
                                                                              \
template                                                                      \
<template<class> class Field, class Form, class Cmpt, int nCmpt, class Type>  \
tmp<FieldField<Field, typename product<Form, Type>::type> >                   \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const tmp<FieldField<Field, Type> >& tf1                                  \
);

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
