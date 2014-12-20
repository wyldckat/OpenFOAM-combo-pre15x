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
    High performance macro functions for Field<Type> algebra.  These expand
    using either array element access (for vector machines) or pointer
    dereferencing for scalar machines as appropriate.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(ReturnType, Type1, Func)                               \
                                                                              \
void Func(Field<ReturnType>& f, const UList<Type1>& f1);                      \
tmp<Field<ReturnType> > Func(const UList<Type1>& f1);

#define UNARY_FUNCTION_N(ReturnType, Type1, Func)                             \
                                                                              \
UNARY_FUNCTION(ReturnType, Type1, Func)                                       \
tmp<Field<ReturnType> > Func(const tmp<Field<Type1> >& tf1);

#define UNARY_FUNCTION_R(ReturnType, Type1, Func)                             \
                                                                              \
UNARY_FUNCTION(ReturnType, Type1, Func)                                       \
tmp<Field<ReturnType> > Func(const tmp<Field<Type1> >& tf1);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                         \
                                                                              \
void OpFunc(Field<ReturnType>& f, const UList<Type1>& f1);                    \
tmp<Field<ReturnType> > operator Op(const UList<Type1>& f1);

#define UNARY_OPERATOR_N(ReturnType, Type1, Op, OpFunc)                       \
                                                                              \
UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                                 \
tmp<Field<ReturnType> > operator Op(const tmp<Field<Type1> >& tf1);

#define UNARY_OPERATOR_R(ReturnType, Type1, Op, OpFunc)                       \
                                                                              \
UNARY_OPERATOR(ReturnType, Type1, Op, OpFunc)                                 \
tmp<Field<ReturnType> > operator Op(const tmp<Field<Type1> >& tf1);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION_FF(ReturnType, Type1, Type2, Func)                    \
                                                                              \
void Func                                                                     \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
);                                                                            \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
);

#define BINARY_FUNCTION_FT(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_FUNCTION_FR(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
);


#define BINARY_FUNCTION_TF(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
);

#define BINARY_FUNCTION_TT(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_FUNCTION_TR(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);


#define BINARY_FUNCTION_RF(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
);

#define BINARY_FUNCTION_RT(ReturnType, Type1, Type2, Func)                    \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);


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
);                                                                            \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const Type1& s1,                                                          \
    const UList<Type2>& f2                                                    \
);

#define BINARY_TYPE_FUNCTION_ST(ReturnType, Type1, Type2, Func)               \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_TYPE_FUNCTION_SR(ReturnType, Type1, Type2, Func)               \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
);


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)               \
                                                                              \
void Func                                                                     \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
);                                                                            \
                                                                              \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
);

#define BINARY_TYPE_FUNCTION_TS(ReturnType, Type1, Type2, Func)               \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
);

#define BINARY_TYPE_FUNCTION_RS(ReturnType, Type1, Type2, Func)               \
tmp<Field<ReturnType> > Func                                                  \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
);


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
);                                                                            \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const UList<Type2>& f2                                                    \
);

#define BINARY_OPERATOR_FT(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_FR(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const tmp<Field<Type2> >& tf2                                             \
);


#define BINARY_OPERATOR_TF(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
);

#define BINARY_OPERATOR_TT(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_OPERATOR_TR(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);


#define BINARY_OPERATOR_RF(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const UList<Type2>& f2                                                    \
);

#define BINARY_OPERATOR_RT(ReturnType, Type1, Type2, Op, OpFunc)              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const tmp<Field<Type2> >& tf2                                             \
);


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
);                                                                            \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const Type1& s1,                                                          \
    const UList<Type2>& f2                                                    \
);

#define BINARY_TYPE_OPERATOR_ST(ReturnType, Type1, Type2, Op, OpFunc)         \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
);

#define BINARY_TYPE_OPERATOR_SR(ReturnType, Type1, Type2, Op, OpFunc)         \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const Type1& s1,                                                          \
    const tmp<Field<Type2> >& tf2                                             \
);


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)         \
                                                                              \
void OpFunc                                                                   \
(                                                                             \
    Field<ReturnType>& f,                                                     \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
);                                                                            \
                                                                              \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const UList<Type1>& f1,                                                   \
    const Type2& s2                                                           \
);

#define BINARY_TYPE_OPERATOR_TS(ReturnType, Type1, Type2, Op, OpFunc)         \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
);

#define BINARY_TYPE_OPERATOR_RS(ReturnType, Type1, Type2, Op, OpFunc)         \
tmp<Field<ReturnType> > operator Op                                           \
(                                                                             \
    const tmp<Field<Type1> >& tf1,                                            \
    const Type2& s2                                                           \
);


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


// ************************************************************************* //
