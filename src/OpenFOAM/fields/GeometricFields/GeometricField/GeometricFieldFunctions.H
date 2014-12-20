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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void component
(
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >& gcf,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const direction d
);

template<class Type, template<class> class PatchField, class GeoMesh>
void T
(
     GeometricField<Type, PatchField, GeoMesh>& gf,
     const GeometricField<Type, PatchField, GeoMesh>& gf1
);

template<class Type, template<class> class PatchField, class GeoMesh, int r>
void pow
(
    GeometricField<typename powProduct<Type, r>::type, PatchField, GeoMesh>& gf,
    const GeometricField<Type, PatchField, GeoMesh>& gf1
);

template<class Type, template<class> class PatchField, class GeoMesh, int r>
tmp<GeometricField<typename powProduct<Type, r>::type, PatchField, GeoMesh> >
pow
(
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    typename powProduct<Type, r>::type
);

template<class Type, template<class> class PatchField, class GeoMesh, int r>
tmp<GeometricField<typename powProduct<Type, r>::type, PatchField, GeoMesh> >
pow
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf,
    typename powProduct<Type, r>::type
);

template<class Type, template<class> class PatchField, class GeoMesh>
void sqr
(
    GeometricField
    <typename outerProduct<Type, Type>::type, PatchField, GeoMesh>& gf,
    const GeometricField<Type, PatchField, GeoMesh>& gf1
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        PatchField,
        GeoMesh
    >
>
sqr(const GeometricField<Type, PatchField, GeoMesh>& gf);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename outerProduct<Type, Type>::type,
        PatchField,
        GeoMesh
    >
>
sqr(const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf);

template<class Type, template<class> class PatchField, class GeoMesh>
void magSqr
(
    GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const GeometricField<Type, PatchField, GeoMesh>& gf
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > magSqr
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > magSqr
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
);

template<class Type, template<class> class PatchField, class GeoMesh>
void mag
(
    GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const GeometricField<Type, PatchField, GeoMesh>& gf
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > mag
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > mag
(
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf
);

template<class Type, template<class> class PatchField, class GeoMesh>
void cmptAv
(
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >& gcf,
    const GeometricField<Type, PatchField, GeoMesh>& gf
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
cmptAv(const GeometricField<Type, PatchField, GeoMesh>& gf);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
cmptAv(const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf);

#define BINARY_FUNCTION(func)                                                 \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void func                                                                     \
(                                                                             \
    GeometricField<Type, PatchField, GeoMesh>& gf,                            \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type, PatchField, GeoMesh>& gf2                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const GeometricField<Type, PatchField, GeoMesh>& gf2                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf2               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const GeometricField<Type, PatchField, GeoMesh>& gf2                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf2               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void func                                                                     \
(                                                                             \
    GeometricField<Type, PatchField, GeoMesh>& gf,                            \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type>& dt                                               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Type>& dt                                               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const Type& t                                                             \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const dimensioned<Type>& dt                                               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > func                          \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const Type& t                                                             \
);

BINARY_FUNCTION(max)
BINARY_FUNCTION(min)
BINARY_FUNCTION(scale)

#undef BINARY_FUNCTION


#define UNARY_FUNCTION_WITH_BOUNDARY(returnType, func, gFunc)                 \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
dimensioned<returnType> func                                                  \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf                       \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
dimensioned<returnType> func                                                  \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);

UNARY_FUNCTION_WITH_BOUNDARY(Type, max, gMax)
UNARY_FUNCTION_WITH_BOUNDARY(Type, min, gMin)

#undef UNARY_FUNCTION_WITH_BOUNDARY


#define UNARY_FUNCTION(returnType, func, gFunc)                               \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
dimensioned<returnType> func                                                  \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf                       \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
dimensioned<returnType> func                                                  \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);

UNARY_FUNCTION(Type, sum, gSum)
UNARY_FUNCTION(scalar, sumMag, gSumMag)
UNARY_FUNCTION(Type, average, gAverage)

#undef UNARY_FUNCTION


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

#define UNARY_OPERATOR(op, opFunc)                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void opFunc                                                                   \
(                                                                             \
    GeometricField<Type, PatchField, GeoMesh>& gf,                            \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);

UNARY_OPERATOR(-, negate)

#undef UNARY_OPERATOR


#define BINARY_OPERATOR_FF(Type1, Type2, op, opFunc)                          \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void opFunc                                                                   \
(                                                                             \
    GeometricField<Type, PatchField, GeoMesh>& gf,                            \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);

#define BINARY_OPERATOR_FTR(Type1, Type2, op, opFunc)                         \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2              \
);

#define BINARY_OPERATOR_FT(Type1, Type2, op, opFunc)                          \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2              \
);

#define BINARY_OPERATOR_TRF(Type1, Type2, op, opFunc)                         \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);

#define BINARY_OPERATOR_TF(Type1, Type2, op, opFunc)                          \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);

#define BINARY_OPERATOR_TRT(Type1, Type2, op, opFunc)                         \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2              \
);

#define BINARY_OPERATOR_TTR(Type1, Type2, op, opFunc)                         \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2              \
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
template<class Type, template<class> class PatchField, class GeoMesh>         \
void opFunc                                                                   \
(                                                                             \
    GeometricField<Type, PatchField, GeoMesh>& gf,                            \
    const dimensioned<TYPE>& dt,                                              \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const dimensioned<TYPE>& dt,                                              \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const TYPE& t,                                                            \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const dimensioned<TYPE>& dt,                                              \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const TYPE& t,                                                            \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);


#define BINARY_TYPE_OPERATOR_FT(TYPE, op, opFunc)                             \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
void opFunc                                                                   \
(                                                                             \
    GeometricField<Type, PatchField, GeoMesh>& gf,                            \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<TYPE>& dt                                               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<TYPE>& dt                                               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const TYPE& t                                                             \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const dimensioned<TYPE>& dt                                               \
);                                                                            \
                                                                              \
template<class Type, template<class> class PatchField, class GeoMesh>         \
tmp<GeometricField<Type, PatchField, GeoMesh> > operator op                   \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const TYPE& t                                                             \
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
template                                                                      \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>   \
void opFunc                                                                   \
(                                                                             \
    GeometricField                                                            \
    <typename product<Type1, Type2>::type, PatchField, GeoMesh>& gf,          \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);                                                                            \
                                                                              \
template                                                                      \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>   \
tmp                                                                           \
<                                                                             \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh> \
>                                                                             \
operator op                                                                   \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);                                                                            \
                                                                              \
template                                                                      \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>   \
tmp                                                                           \
<                                                                             \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh> \
>                                                                             \
operator op                                                                   \
(                                                                             \
    const GeometricField<Type1, PatchField, GeoMesh>& gf1,                    \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2              \
);                                                                            \
                                                                              \
template                                                                      \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>   \
tmp                                                                           \
<                                                                             \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh> \
>                                                                             \
operator op                                                                   \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,             \
    const GeometricField<Type2, PatchField, GeoMesh>& gf2                     \
);                                                                            \
                                                                              \
template                                                                      \
<class Type1, class Type2, template<class> class PatchField, class GeoMesh>   \
tmp                                                                           \
<                                                                             \
    GeometricField<typename product<Type1, Type2>::type, PatchField, GeoMesh> \
>                                                                             \
operator op                                                                   \
(                                                                             \
    const tmp<GeometricField<Type1, PatchField, GeoMesh> >& tgf1,             \
    const tmp<GeometricField<Type2, PatchField, GeoMesh> >& tgf2              \
);                                                                            \
                                                                              \
template                                                                      \
<class Form, class Type, template<class> class PatchField, class GeoMesh>     \
void opFunc                                                                   \
(                                                                             \
    GeometricField                                                            \
    <typename product<Type, Form>::type, PatchField, GeoMesh>& gf,            \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Form>& dvs                                              \
);                                                                            \
                                                                              \
template                                                                      \
<class Form, class Type, template<class> class PatchField, class GeoMesh>     \
tmp<GeometricField<typename product<Type, Form>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const dimensioned<Form>& dvs                                              \
);                                                                            \
                                                                              \
template                                                                      \
<                                                                             \
    class Form,                                                               \
    class Cmpt,                                                               \
    int nCmpt,                                                                \
    class Type, template<class> class PatchField,                             \
    class GeoMesh                                                             \
>                                                                             \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1,                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template                                                                      \
<class Form, class Type, template<class> class PatchField, class GeoMesh>     \
tmp<GeometricField<typename product<Type, Form>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const dimensioned<Form>& dvs                                              \
);                                                                            \
                                                                              \
template                                                                      \
<                                                                             \
    class Form,                                                               \
    class Cmpt,                                                               \
    int nCmpt,                                                                \
    class Type, template<class> class PatchField,                             \
    class GeoMesh                                                             \
>                                                                             \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1,              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
                                                                              \
template                                                                      \
<class Form, class Type, template<class> class PatchField, class GeoMesh>     \
void opFunc                                                                   \
(                                                                             \
    GeometricField                                                            \
    <typename product<Form, Type>::type, PatchField, GeoMesh>& gf,            \
    const dimensioned<Form>& dvs,                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template                                                                      \
<class Form, class Type, template<class> class PatchField, class GeoMesh>     \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const dimensioned<Form>& dvs,                                             \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template                                                                      \
<                                                                             \
    class Form,                                                               \
    class Cmpt,                                                               \
    int nCmpt,                                                                \
    class Type, template<class> class PatchField,                             \
    class GeoMesh                                                             \
>                                                                             \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const GeometricField<Type, PatchField, GeoMesh>& gf1                      \
);                                                                            \
                                                                              \
template                                                                      \
<class Form, class Type, template<class> class PatchField, class GeoMesh>     \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const dimensioned<Form>& dvs,                                             \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);                                                                            \
                                                                              \
template                                                                      \
<                                                                             \
    class Form,                                                               \
    class Cmpt,                                                               \
    int nCmpt,                                                                \
    class Type, template<class> class PatchField,                             \
    class GeoMesh                                                             \
>                                                                             \
tmp<GeometricField<typename product<Form, Type>::type, PatchField, GeoMesh> > \
operator op                                                                   \
(                                                                             \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& tgf1               \
);

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
