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

#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& sf,
    const UList<Type>& f,
    const direction d
);


template<class Type>
void T(Field<Type>& f1, const UList<Type>& f2);


template<class Type, int r>
void pow
(
    Field<typename powProduct<Type, r>::type>& f,
    const UList<Type>& vf
);


template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
    const UList<Type>& f, typename powProduct<Type, r>::type
  = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
    const tmp<Field<Type> >& tf, typename powProduct<Type, r>::type
  = pTraits<typename powProduct<Type, r>::type>::zero
);


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& f,
    const UList<Type>& vf
);

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const UList<Type>& f);

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const tmp<Field<Type> >& tf);


template<class Type>
void magSqr(Field<scalar>& sf, const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > magSqr(const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > magSqr(const tmp<Field<Type> >& tf);


template<class Type>
void mag(Field<scalar>& sf, const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > mag(const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > mag(const tmp<Field<Type> >& tf);


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& cf, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const tmp<Field<Type> >& tf);


template<class Type>
void cmptMag(Field<Type>& cf, const UList<Type>& f);

template<class Type>
tmp<Field<Type> > cmptMag(const UList<Type>& f);

template<class Type>
tmp<Field<Type> > cmptMag(const tmp<Field<Type> >& tf);


#define BINARY_FUNCTION(Func)                                                 \
                                                                              \
template<class Type>                                                          \
void Func(Field<Type>& f, const UList<Type>& f1, const UList<Type>& f2);      \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > Func(const UList<Type>& f1, const UList<Type>& f2);         \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > Func(const UList<Type>& f1, const tmp<Field<Type> >& tf2);  \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > Func(const tmp<Field<Type> >& tf1, const UList<Type>& f2);  \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > Func                                                        \
(                                                                             \
    const tmp<Field<Type> >& tf1,                                             \
    const tmp<Field<Type> >& tf2                                              \
);                                                                            \
                                                                              \
template<class Type>                                                          \
void Func(Field<Type>& f, const UList<Type>& f1, const Type& s);              \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > Func(const UList<Type>& f1, const Type& s);                 \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > Func(const tmp<Field<Type> >& tf1, const Type& s);

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(scale)

#undef BINARY_FUNCTION


#define TMP_UNARY_FUNCTION(ReturnType, Func)                                  \
                                                                              \
template<class Type>                                                          \
ReturnType Func(const tmp<Field<Type> >& tf1);

template<class Type>
Type max(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, sum)


template<class Type>
scalar sumProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
scalar sumSqr(const UList<Type>& f);

TMP_UNARY_FUNCTION(scalar, sumSqr)

template<class Type>
scalar sumMag(const UList<Type>& f);

TMP_UNARY_FUNCTION(scalar, sumMag)

template<class Type>
Type average(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, average)


#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                      \
                                                                              \
template<class Type>                                                          \
ReturnType gFunc(const UList<Type>& f);                                       \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION

template<class Type>
scalar gSumProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
Type gAverage(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(Op, OpFunc)                                            \
                                                                              \
template<class Type>                                                          \
void OpFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const UList<Type>& f1                                                     \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const UList<Type>& f1                                                     \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const tmp<Field<Type> >& tf1                                              \
);

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, Op, OpFunc)                          \
                                                                              \
template<class Type>                                                          \
void OpFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
);


#define BINARY_OPERATOR_FR(Type1, Type2, Op, OpFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_FT(Type1, Type2, Op, OpFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_RF(Type1, Type2, Op, OpFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
);

#define BINARY_OPERATOR_TF(Type1, Type2, Op, OpFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
);

#define BINARY_OPERATOR_RT(Type1, Type2, Op, OpFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_TR(Type1, Type2, Op, OpFunc)                          \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_RR(Type1, Type2, Op, OpFunc)                          \
    BINARY_OPERATOR_FF(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_FR(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_RF(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_RT(Type1, Type2, Op, OpFunc)

#define BINARY_OPERATOR_RN(Type1, Type2, Op, OpFunc)                          \
    BINARY_OPERATOR_FF(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_FT(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_RF(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_RT(Type1, Type2, Op, OpFunc)

#define BINARY_OPERATOR_NR(Type1, Type2, Op, OpFunc)                          \
    BINARY_OPERATOR_FF(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_FR(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_TF(Type1, Type2, Op, OpFunc)                              \
    BINARY_OPERATOR_TR(Type1, Type2, Op, OpFunc)

    //BINARY_OPERATOR_RR(Type, Type, +, add)
    //BINARY_OPERATOR_RR(Type, Type, -, subtract)
BINARY_OPERATOR_RN(Type, scalar, *, multiply)
BINARY_OPERATOR_NR(scalar, Type, *, multiply)
BINARY_OPERATOR_RN(Type, scalar, /, divide)

#undef BINARY_OPERATOR_RR
#undef BINARY_OPERATOR_RN
#undef BINARY_OPERATOR_NR
#undef BINARY_OPERATOR_FF
#undef BINARY_OPERATOR_FR
#undef BINARY_OPERATOR_TF
#undef BINARY_OPERATOR_TR
#undef BINARY_OPERATOR_FT
#undef BINARY_OPERATOR_RF
#undef BINARY_OPERATOR_RT


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(TYPE, Op, OpFunc)                             \
                                                                              \
template<class Type>                                                          \
void OpFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const TYPE& s,                                                            \
    const UList<Type>& f1                                                     \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const TYPE& s,                                                            \
    const UList<Type>& f1                                                     \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const TYPE& s,                                                            \
    const tmp<Field<Type> >& tf1                                              \
);


#define BINARY_TYPE_OPERATOR_FS(TYPE, Op, OpFunc)                             \
                                                                              \
template<class Type>                                                          \
void OpFunc                                                                   \
(                                                                             \
    Field<Type>& f,                                                           \
    const UList<Type>& f1,                                                    \
    const TYPE& s                                                             \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const UList<Type>& f1,                                                    \
    const TYPE& s                                                             \
);                                                                            \
                                                                              \
template<class Type>                                                          \
tmp<Field<Type> > operator Op                                                 \
(                                                                             \
    const tmp<Field<Type> >& tf1,                                             \
    const TYPE& s                                                             \
);


#define BINARY_TYPE_OPERATOR(TYPE, Op, OpFunc)                                \
    BINARY_TYPE_OPERATOR_SF(TYPE, Op, OpFunc)                                 \
    BINARY_TYPE_OPERATOR_FS(TYPE, Op, OpFunc)

    //BINARY_TYPE_OPERATOR(Type, +, add)
    //BINARY_TYPE_OPERATOR(Type, -, subtract)

BINARY_TYPE_OPERATOR(scalar, *, multiply)
BINARY_TYPE_OPERATOR_FS(scalar, /, divide)

#undef BINARY_TYPE_OPERATOR
#undef BINARY_TYPE_OPERATOR_SF
#undef BINARY_TYPE_OPERATOR_FS


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                 \
                                                                              \
template<class Type1, class Type2>                                            \
void OpFunc                                                                   \
(                                                                             \
    Field<typename product<Type1, Type2>::type>& f,                           \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
);                                                                            \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const UList<Type1>& f1, const UList<Type2>& f2);                  \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const UList<Type1>& f1, const tmp<Field<Type2> >& tf2);           \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<Field<Type1> >& tf1, const UList<Type2>& f2);           \
                                                                              \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<Field<Type1> >& tf1, const tmp<Field<Type2> >& tf2);    \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
void OpFunc                                                                   \
(                                                                             \
    Field<typename product<Type, Form>::type>& f,                             \
    const UList<Type>& f1,                                                    \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator Op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs);   \
                                                                              \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator Op(const tmp<Field<Type> >&tf1,const VectorSpace<Form,Cmpt,nCmpt>&vs);\
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
void OpFunc                                                                   \
(                                                                             \
    Field<typename product<Form, Type>::type>& f,                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const UList<Type>& f1                                                     \
);                                                                            \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1);   \
                                                                              \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&vs,const tmp<Field<Type> >&tf1);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
