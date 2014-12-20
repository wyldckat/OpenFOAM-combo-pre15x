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

\*---------------------------------------------------------------------------*/

#ifndef FieldFunctionM_C
#define FieldFunctionM_C

#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func)                               \
                                                                              \
void Func(Field<ReturnType>& f, const UList<Type1>& f1)                       \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F(ReturnType, f, =, ::Foam::Func, Type1, f1)           \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > Func(const UList<Type1>& f1)                          \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f1.size()));             \
    TFOR_ALL_F_OP_FUNC_F(ReturnType, tf(), =, ::Foam::Func, Type1, f1)        \
    return tf;                                                                \
}


#define UNARY_FUNCTION_N(ReturnType, Type1, Func)                             \
                                                                              \
UNARY_FUNCTION(ReturnType, Type1, Func)                                       \
                                                                              \
tmp<Field<ReturnType> > Func(const tmp<Field<Type1> >& tf1)                   \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(tf1().size()));          \
    TFOR_ALL_F_OP_FUNC_F(ReturnType, tf(), =, ::Foam::Func, Type1, tf1())     \
    tf1.clear();                                                              \
    return tf;                                                                \
}


#define UNARY_FUNCTION_R(ReturnType, Type1, Func)                             \
                                                                              \
UNARY_FUNCTION(ReturnType, Type1, Func)                                       \
                                                                              \
tmp<Field<ReturnType> > Func(const tmp<Field<Type1> >& tf1)                   \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    TFOR_ALL_F_OP_FUNC_F(ReturnType, tf(), =, ::Foam::Func, Type1, tf())      \
    return tf;                                                                \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                         \
                                                                              \
void OpFunc                                                                   \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1                                                    \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_OP_F(ReturnType, f, =, Op, Type1, f1)                       \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f1.size()));             \
    OpFunc(tf(), f1);                                                         \
    return tf;                                                                \
}                                                                             \


#define UNARY_OPERATOR_N(ReturnType, Type1, Op, OpFunc)                       \
                                                                              \
UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                                 \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(tf1().size()));          \
    OpFunc(tf(), tf1());                                                      \
    tf1.clear();                                                              \
    return tf;                                                                \
}


#define UNARY_OPERATOR_R(ReturnType, Type1, Op, OpFunc)                       \
                                                                              \
UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                                 \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    OpFunc(tf(), tf());                                                       \
    return tf;                                                                \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION_FF(ReturnType, Type1, Type2, Func)                    \
                                                                              \
void Func                                                                     \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F_F                                                    \
        (ReturnType, f, =, ::Foam::Func, Type, f1, Type, f2)                  \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f1.size()));             \
    Func(tf(), f1, f2);                                                       \
    return tf;                                                                \
}

#define BINARY_FUNCTION_FT(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf = Func(f1, tf2());                             \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_FUNCTION_FR(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf2.ptr());                                    \
    Func(tf(), f1, tf());                                                     \
    return tf;                                                                \
}



#define BINARY_FUNCTION_TF(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf = Func(tf1(), f2);                             \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_FUNCTION_TT(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf = Func(tf1(), tf2());                          \
    tf1.clear();                                                              \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_FUNCTION_TR(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf2.ptr());                                    \
    Func(tf(), tf1(), tf());                                                  \
    tf1.clear();                                                              \
    return tf;                                                                \
}


#define BINARY_FUNCTION_RF(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    Func(tf(), tf(), f2);                                                     \
    return tf;                                                                \
}

#define BINARY_FUNCTION_RT(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    Func(tf(), tf(), tf2());                                                  \
    tf2.clear();                                                              \
    return tf;                                                                \
}


#define BINARY_FUNCTION_NN(ReturnType, Type1, Type2, Func)                    \
    BINARY_FUNCTION_FF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_FT(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_TF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_TT(ReturnType, Type1, Type2, Func)

#define BINARY_FUNCTION_NR(ReturnType, Type1, Type2, Func)                    \
    BINARY_FUNCTION_FF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_FR(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_TF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_TR(ReturnType, Type1, Type2, Func)

#define BINARY_FUNCTION_RN(ReturnType, Type1, Type2, Func)                    \
    BINARY_FUNCTION_FF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_FT(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_RF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_RT(ReturnType, Type1, Type2, Func)

#define BINARY_FUNCTION_RR(ReturnType, Type1, Type2, Func)                    \
    BINARY_FUNCTION_FF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_FR(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_RF(ReturnType, Type1, Type2, Func)                        \
    BINARY_FUNCTION_RT(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)               \
                                                                              \
void Func                                                                     \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const Type1& s1,                                                          \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_S_F                                                    \
        (ReturnType, f, =, ::Foam::Func, Type1, s1, Type2, f2)                \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const Type1& s1,                                                          \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f2.size()));             \
    Func(tf(), s1, f2);                                                       \
    return tf;                                                                \
}


#define BINARY_TYPE_FUNCTION_ST(ReturnType, Type1, Type2, Func)               \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(tf2().size()));          \
    Func(tf(), s1, tf2());                                                    \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_TYPE_FUNCTION_SR(ReturnType, Type1, Type2, Func)               \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf2.ptr());                                    \
    Func(tf(), s1, tf());                                                     \
    return tf;                                                                \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)               \
                                                                              \
void Func                                                                     \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F_S                                                    \
        (ReturnType, f, =, ::Foam::Func, Type1, f1, Type2, s2)                \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f1.size()));             \
    Func(tf(), f1, s2);                                                       \
    return tf;                                                                \
}

#define BINARY_TYPE_FUNCTION_TS(ReturnType, Type1, Type2, Func)               \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(tf1().size()));          \
    Func(tf(), tf1(), s2);                                                    \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_TYPE_FUNCTION_RS(ReturnType, Type1, Type2, Func)               \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    Func(tf(), tf(), s2);                                                     \
    return tf;                                                                \
}


#define BINARY_TYPE_FUNCTION_NN(ReturnType, Type1, Type2, Func)               \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_ST(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_TS(ReturnType, Type1, Type2, Func)

#define BINARY_TYPE_FUNCTION_NR(ReturnType, Type1, Type2, Func)               \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SR(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_TS(ReturnType, Type1, Type2, Func)

#define BINARY_TYPE_FUNCTION_RN(ReturnType, Type1, Type2, Func)               \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_ST(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_RS(ReturnType, Type1, Type2, Func)

#define BINARY_TYPE_FUNCTION_RR(ReturnType, Type1, Type2, Func)               \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SR(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_RS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR_FF(ReturnType, Type1, Type2, Op, OpFunc)              \
                                                                              \
void OpFunc                                                                   \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_F_OP_F(ReturnType, f, =, Type1, f1, Op, Type2, f2)          \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f1.size()));             \
    OpFunc(tf(), f1, f2);                                                     \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FT(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf = f1 Op tf2();                                 \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_FR(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf2.ptr());                                    \
    OpFunc(tf(), f1, tf());                                                   \
    return tf;                                                                \
}



#define BINARY_OPERATOR_TF(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf = tf1() Op f2;                                 \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TT(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf = tf1() Op tf2();                              \
    tf1.clear();                                                              \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_OPERATOR_TR(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf2.ptr());                                    \
    OpFunc(tf(), tf1(), tf());                                                \
    tf1.clear();                                                              \
    return tf;                                                                \
}


#define BINARY_OPERATOR_RF(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    OpFunc(tf(), tf(), f2);                                                   \
    return tf;                                                                \
}

#define BINARY_OPERATOR_RT(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    OpFunc(tf(), tf(), tf2());                                                \
    tf2.clear();                                                              \
    return tf;                                                                \
}


#define BINARY_OPERATOR_NN(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_OPERATOR_FF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_FT(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_TF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_TT(ReturnType, Type1, Type2, Op, OpFunc)

#define BINARY_OPERATOR_NR(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_OPERATOR_FF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_FR(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_TF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_TR(ReturnType, Type1, Type2, Op, OpFunc)

#define BINARY_OPERATOR_RN(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_OPERATOR_FF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_FT(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_RF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_RT(ReturnType, Type1, Type2, Op, OpFunc)

#define BINARY_OPERATOR_RR(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_OPERATOR_FF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_FR(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_RF(ReturnType, Type1, Type2, Op, OpFunc)                  \
    BINARY_OPERATOR_RT(ReturnType, Type1, Type2, Op, OpFunc)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
void OpFunc                                                                   \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const Type1& s1,                                                          \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_S_OP_F(ReturnType, f, =, Type1, s1, Op, Type2, f2)          \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const Type1& s1,                                                          \
    const UList<Type2>& f2                                                    \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f2.size()));             \
    OpFunc(tf(), s1, f2);                                                     \
    return tf;                                                                \
}


#define BINARY_TYPE_OPERATOR_ST(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(tf2().size()));          \
    OpFunc(tf(), s1, tf2());                                                  \
    tf2.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_TYPE_OPERATOR_SR(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf2.ptr());                                    \
    OpFunc(tf(), s1, tf());                                                   \
    return tf;                                                                \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
void OpFunc                                                                   \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    TFOR_ALL_F_OP_F_OP_S(ReturnType, f, =, Type1, f1, Op, Type2, s2)          \
}                                                                             \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(f1.size()));             \
    OpFunc(tf(), f1, s2);                                                     \
    return tf;                                                                \
}

#define BINARY_TYPE_OPERATOR_TS(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(new Field<ReturnType>(tf1().size()));          \
    OpFunc(tf(), tf1(), s2);                                                  \
    tf1.clear();                                                              \
    return tf;                                                                \
}

#define BINARY_TYPE_OPERATOR_RS(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
)                                                                             \
{                                                                             \
    tmp<Field<ReturnType> > tf(tf1.ptr());                                    \
    OpFunc(tf(), tf(), s2);                                                   \
    return tf;                                                                \
}


#define BINARY_TYPE_OPERATOR_NN(ReturnType, Type1, Type2, Op, OpFunc)         \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_ST(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_TS(ReturnType, Type1, Type2, Op, OpFunc)

#define BINARY_TYPE_OPERATOR_NR(ReturnType, Type1, Type2, Op, OpFunc)         \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_TS(ReturnType, Type1, Type2, Op, OpFunc)

#define BINARY_TYPE_OPERATOR_RN(ReturnType, Type1, Type2, Op, OpFunc)         \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_ST(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_RS(ReturnType, Type1, Type2, Op, OpFunc)

#define BINARY_TYPE_OPERATOR_RR(ReturnType, Type1, Type2, Op, OpFunc)         \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_RS(ReturnType, Type1, Type2, Op, OpFunc)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
