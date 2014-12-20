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

Class
    GeometricScalarField

Description
    Scalar specific part of the implementation of GeometricField.

SourceFiles
    GeometricScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricScalarField_H
#define GeometricScalarField_H

#include "GeometricField.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void stabilise
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const dimensioned<scalar>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > stabilise
(
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const dimensioned<scalar>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > stabilise
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&,
    const dimensioned<scalar>&
);


template<template<class> class PatchField, class GeoMesh>
void divide
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const dimensioned<scalar>&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const dimensioned<scalar>&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const dimensioned<scalar>&,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const scalar&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const scalar&,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const dimensioned<scalar>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const dimensioned<scalar>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&,
    const dimensioned<scalar>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>&,
    const scalar&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&,
    const scalar&
);


template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const dimensioned<scalar>&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const dimensioned<scalar>&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const dimensioned<scalar>&,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const scalar&,
    const GeometricField<scalar, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const scalar&,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&
);


#define transFunc(func)                                                        \
template<template<class> class PatchField, class GeoMesh>                      \
void func                                                                      \
(                                                                              \
    GeometricField<scalar, PatchField, GeoMesh>& Res,                          \
    const GeometricField<scalar, PatchField, GeoMesh>& sf                      \
);                                                                             \
                                                                               \
template<template<class> class PatchField, class GeoMesh>                      \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                         \
(                                                                              \
    const GeometricField<scalar, PatchField, GeoMesh>&                         \
);                                                                             \
                                                                               \
template<template<class> class PatchField, class GeoMesh>                      \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                         \
(                                                                              \
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&                   \
);

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
template<template<class> class PatchField, class GeoMesh>                      \
void func                                                                      \
(                                                                              \
    GeometricField<scalar, PatchField, GeoMesh>& Res,                          \
    const int n,                                                               \
    const GeometricField<scalar, PatchField, GeoMesh>& sf                      \
);                                                                             \
                                                                               \
template<template<class> class PatchField, class GeoMesh>                      \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                         \
(                                                                              \
    const int n,                                                               \
    const GeometricField<scalar, PatchField, GeoMesh>&                         \
);                                                                             \
                                                                               \
template<template<class> class PatchField, class GeoMesh>                      \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                         \
(                                                                              \
    const int n,                                                               \
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >&                   \
);

transFunc(jn)
transFunc(yn)

#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GeometricScalarField.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
