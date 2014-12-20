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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    scalarFieldField

Description
    Specialisation of FieldField<T> for scalar.

SourceFiles
    scalarFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef scalarFieldField_H
#define scalarFieldField_H

#include "FieldField.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void stabilise
(
    FieldField<Field, scalar>& Res,
    const FieldField<Field, scalar>& sf,
    const scalar s
);

tmp<FieldField<Field, scalar> > stabilise
(
    const FieldField<Field, scalar>&,
    const scalar s
);

tmp<FieldField<Field, scalar> > stabilise
(
    const tmp<FieldField<Field, scalar> >&,
    const scalar s
);


template<template<class> class Field>
void divide
(
    FieldField<Field, scalar>&,
    const scalar,
    const FieldField<Field, scalar>&
);


template<template<class> class Field>
tmp<FieldField<Field, scalar> > operator/
(
    const scalar,
    const FieldField<Field, scalar>&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > operator/
(
    const scalar,
    const tmp<FieldField<Field, scalar> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Field>
void pow
(
    FieldField<Field, scalar>&,
    const FieldField<Field, scalar>&,
    const FieldField<Field, scalar>&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const FieldField<Field, scalar>&,
    const FieldField<Field, scalar>&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const tmp<FieldField<Field, scalar> >&,
    const FieldField<Field, scalar>&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const FieldField<Field, scalar>&,
    const tmp<FieldField<Field, scalar> >&
);

tmp<FieldField<Field, scalar> > pow
(
    const tmp<FieldField<Field, scalar> >&,
    const tmp<FieldField<Field, scalar> >&
);


template<template<class> class Field>
void pow
(
    FieldField<Field, scalar>&,
    const FieldField<Field, scalar>&,
    const scalar&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const FieldField<Field, scalar>&,
    const scalar&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const tmp<FieldField<Field, scalar> >&,
    const scalar&
);


template<template<class> class Field>
void pow
(
    FieldField<Field, scalar>&,
    const scalar&,
    const FieldField<Field, scalar>&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const scalar&,
    const FieldField<Field, scalar>&
);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > pow
(
    const scalar&,
    const tmp<FieldField<Field, scalar> >&
);


#define transFunc(func)                                                        \
template<template<class> class Field>                                          \
void func(FieldField<Field, scalar>& Res, const FieldField<Field, scalar>& sf);\
                                                                               \
template<template<class> class Field>                                          \
tmp<FieldField<Field, scalar> > func(const FieldField<Field, scalar>&);        \
                                                                               \
template<template<class> class Field>                                          \
tmp<FieldField<Field, scalar> > func(const tmp<FieldField<Field, scalar> >&);

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


#define transFunc(func)                                                        \
template<template<class> class Field>                                          \
void func                                                                      \
(                                                                              \
    FieldField<Field, scalar>& Res,                                            \
    const int n,                                                               \
    const FieldField<Field, scalar>& sf);                                      \
                                                                               \
template<template<class> class Field>                                          \
tmp<FieldField<Field, scalar> > func                                           \
(                                                                              \
    const int n,                                                               \
    const FieldField<Field, scalar>&                                           \
);                                                                             \
                                                                               \
template<template<class> class Field>                                          \
tmp<FieldField<Field, scalar> > func                                           \
(                                                                              \
    const int n,                                                               \
    const tmp<FieldField<Field, scalar> >&                                     \
);

transFunc(jn)
transFunc(yn)

#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "scalarFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
