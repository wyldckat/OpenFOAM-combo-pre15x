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
    Generic TetPointPatchField type.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<template<class> class PatchField, class PointPatch, class Type>
inline void component
(
    TetPointPatchField<PatchField, PointPatch, typename Field<Type>::cmptType>& sf,
    const TetPointPatchField<PatchField, PointPatch, Type>& f,
    const direction d
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void T
(
    TetPointPatchField<PatchField, PointPatch, Type>& f1,
    const TetPointPatchField<PatchField, PointPatch, Type>& f2
){}


template<template<class> class PatchField, class PointPatch, class Type, int r>
inline void pow
(
    Field<typename powProduct<Type, r>::type>& f,
    const TetPointPatchField<PatchField, PointPatch, Type>& vf
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void sqr
(
    Field<typename outerProduct<Type, Type>::type>& f,
    const TetPointPatchField<PatchField, PointPatch, Type>& vf
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void magSqr
(
    TetPointPatchField<PatchField, PointPatch, scalar>& sf,
    const TetPointPatchField<PatchField, PointPatch, Type>& f
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void mag
(
    TetPointPatchField<PatchField, PointPatch, scalar>& sf,
    const TetPointPatchField<PatchField, PointPatch, Type>& f
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void cmptAv
(
    TetPointPatchField<PatchField, PointPatch, typename Field<Type>::cmptType>& cf,
    const TetPointPatchField<PatchField, PointPatch, Type>& f
){}


template<template<class> class PatchField, class PointPatch, class Type>
inline void cmptMag
(
    TetPointPatchField<PatchField, PointPatch, Type>& cf,
    const TetPointPatchField<PatchField, PointPatch, Type>& f
){}


#define BINARY_FUNCTION(func)                                                 \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void func                                                              \
(                                                                             \
    TetPointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1,                  \
    const TetPointPatchField<PatchField, PointPatch, Type>& f2                   \
){}                                                                           \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void func                                                              \
(                                                                             \
    TetPointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1,                  \
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
    TetPointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1                   \
){}

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void opFunc                                                            \
(                                                                             \
    TetPointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TetPointPatchField<PatchField, PointPatch, Type1>& f1,                 \
    const TetPointPatchField<PatchField, PointPatch, Type2>& f2                  \
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
    TetPointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TYPE& s,                                                            \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1                   \
){}


#define BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)                             \
                                                                              \
template<template<class> class PatchField, class PointPatch, class Type>      \
inline void opFunc                                                            \
(                                                                             \
    TetPointPatchField<PatchField, PointPatch, Type>& f,                         \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1,                  \
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
    TetPointPatchField                                                           \
    <PatchField, PointPatch, typename product<Type1, Type2>::type>& f,        \
    const TetPointPatchField<PatchField, PointPatch, Type1>& f1,                 \
    const TetPointPatchField<PatchField, PointPatch, Type2>& f2                  \
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
    TetPointPatchField                                                           \
    <PatchField, PointPatch, typename product<Type, Form>::type>& f,          \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1,                  \
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
    TetPointPatchField                                                           \
    <PatchField, PointPatch, typename product<Form, Type>::type>& f,          \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const TetPointPatchField<PatchField, PointPatch, Type>& f1                   \
){}

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


template<template<class> class PatchField, class PointPatch>
void hdual
(
    TetPointPatchField<PatchField, PointPatch, vector>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void hdual
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, vector>&
){}

template<template<class> class PatchField, class PointPatch>
void diag
(
    TetPointPatchField<PatchField, PointPatch, vector>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void tr
(
    TetPointPatchField<PatchField, PointPatch, scalar>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void dev
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void dev2
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void det
(
    TetPointPatchField<PatchField, PointPatch, scalar>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void inv
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void hinv
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void symm
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void skew
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void eigenValues
(
    TetPointPatchField<PatchField, PointPatch, vector>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}

template<template<class> class PatchField, class PointPatch>
void eigenVectors
(
    TetPointPatchField<PatchField, PointPatch, tensor>&,
    const TetPointPatchField<PatchField, PointPatch, tensor>&
){}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
