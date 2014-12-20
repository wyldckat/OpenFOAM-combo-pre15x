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
    Generic PointPatchField type.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<template<class> class PatchField, class PointPatch, class Type>
inline void component
(
    PointPatchField<PatchField, PointPatch, typename Field<Type>::cmptType>& sf,
    const PointPatchField<PatchField, PointPatch, Type>& f,
    const direction d
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void T
(
    PointPatchField<PatchField, PointPatch, Type>& f1,
    const PointPatchField<PatchField, PointPatch, Type>& f2
){}


template<template<class> class PatchField, class PointPatch, class Type, int r>
inline void pow
(
    Field<typename powProduct<Type, r>::type>& f,
    const PointPatchField<PatchField, PointPatch, Type>& vf
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void sqr
(
    Field<typename outerProduct<Type, Type>::type>& f,
    const PointPatchField<PatchField, PointPatch, Type>& vf
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void magSqr
(
    PointPatchField<PatchField, PointPatch, scalar>& sf,
    const PointPatchField<PatchField, PointPatch, Type>& f
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void mag
(
    PointPatchField<PatchField, PointPatch, scalar>& sf,
    const PointPatchField<PatchField, PointPatch, Type>& f
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void cmptAv
(
    PointPatchField<PatchField, PointPatch, typename Field<Type>::cmptType>& cf,
    const PointPatchField<PatchField, PointPatch, Type>& f
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void cmptMag
(
    PointPatchField<PatchField, PointPatch, Type>& cf,
    const PointPatchField<PatchField, PointPatch, Type>& f
){}


#define BINARY_FUNCTION(func)                                                 \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void func                                                              \
(                                                                             \
    PointPatchField<PatchField, PointPatch, Type>& f,                         \
    const PointPatchField<PatchField, PointPatch, Type>& f1,                  \
    const PointPatchField<PatchField, PointPatch, Type>& f2                   \
){}                                                                           \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void func                                                              \
(                                                                             \
    PointPatchField<PatchField, PointPatch, Type>& f,                         \
    const PointPatchField<PatchField, PointPatch, Type>& f1,                  \
    const Type& s                                                             \
){}

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(scale)

#undef BINARY_FUNCTION


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

#define UNARY_OPERATOR(op, opFunc)                                            \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField<PatchField, PointPatch, Type>& f,                         \
    const PointPatchField<PatchField, PointPatch, Type>& f1                   \
){}

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField<PatchField, PointPatch, Type>& f,                         \
    const PointPatchField<PatchField, PointPatch, Type1>& f1,                 \
    const PointPatchField<PatchField, PointPatch, Type2>& f2                  \
){}

#define BINARY_OPERATOR_R(Type1, Type2, op, opFunc)                           \
    BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)

BINARY_OPERATOR_R(Type, Type, +, add)
BINARY_OPERATOR_R(Type, Type, -, subtract)
BINARY_OPERATOR_FF(scalar, Type, *, multiply)
BINARY_OPERATOR_FF(Type, scalar, /, divide)

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
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TYPE& s,                                                            \
    const PointPatchField<PatchField, PointPatch, Type>& f1                   \
){}


#define BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)                             \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField<PatchField, PointPatch, Type>& f,                         \
    const PointPatchField<PatchField, PointPatch, Type>& f1,                  \
    const TYPE& s                                                             \
){}


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
template                                                                      \
<                                                                             \
    template<class> class PatchField,                                         \
    class PointPatch,                                                         \
    class Type1,                                                              \
    class Type2                                                               \
>                                                                             \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField                                                           \
    <PatchField, PointPatch, typename product<Type1, Type2>::type>& f,        \
    const PointPatchField<PatchField, PointPatch, Type1>& f1,                 \
    const PointPatchField<PatchField, PointPatch, Type2>& f2                  \
){}                                                                           \
                                                                              \
template                                                                      \
<                                                                             \
    template<class> class PatchField,                                         \
    class PointPatch,                                                         \
    class Type,                                                               \
    class Form,                                                               \
    class Cmpt,                                                               \
    int nCmpt                                                                 \
>                                                                             \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField                                                           \
    <PatchField, PointPatch, typename product<Type, Form>::type>& f,          \
    const PointPatchField<PatchField, PointPatch, Type>& f1,                  \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
){}                                                                           \
                                                                              \
template                                                                      \
<                                                                             \
    template<class> class PatchField,                                         \
    class PointPatch,                                                         \
    class Form,                                                               \
    class Cmpt,                                                               \
    int nCmpt,                                                                \
    class Type                                                                \
>                                                                             \
inline void opFunc                                                            \
(                                                                             \
    PointPatchField                                                           \
    <PatchField, PointPatch, typename product<Form, Type>::type>& f,          \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const PointPatchField<PatchField, PointPatch, Type>& f1                   \
){}

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


template<template<class> class PatchField, class PointPatch>
void hdual
(
    PointPatchField<PatchField, PointPatch, vector>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void hdual
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, vector>&
){}

template<template<class> class PatchField, class PointPatch>
void diag
(
    PointPatchField<PatchField, PointPatch, vector>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void tr
(
    PointPatchField<PatchField, PointPatch, scalar>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void dev
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void dev2
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void det
(
    PointPatchField<PatchField, PointPatch, scalar>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void inv
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void hinv
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void symm
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void skew
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void eigenValues
(
    PointPatchField<PatchField, PointPatch, vector>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void eigenVectors
(
    PointPatchField<PatchField, PointPatch, tensor>&,
    const PointPatchField<PatchField, PointPatch, tensor>&
){}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
