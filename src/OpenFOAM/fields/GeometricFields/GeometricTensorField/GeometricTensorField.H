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
    TensorField

Description
    Tensor specific part of the implementation of GeometricField.

SourceFiles
    GeometricTensorField.C


\*---------------------------------------------------------------------------*/

#ifndef GeometricTensorField_H
#define GeometricTensorField_H

#include "GeometricField.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void hdual
(
    GeometricField<vector, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > operator*
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > operator*
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void hdual
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<vector, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > operator*
(
    const GeometricField<vector, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > operator*
(
    const tmp<GeometricField<vector, PatchField, GeoMesh> >&
);


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void diag
(
    GeometricField<vector, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > diag
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > diag
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void tr
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > tr
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > tr
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void dev
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void dev2
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev2
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev2
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void det
(
    GeometricField<scalar, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > det
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > det
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void inv
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > inv
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > inv
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void hinv
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > hinv
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > hinv
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void symm
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > symm
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > symm
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void skew
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > skew
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > skew
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void eigenValues
(
    GeometricField<vector, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > eigenValues
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > eigenValues
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


template<template<class> class PatchField, class GeoMesh>
void eigenVectors
(
    GeometricField<tensor, PatchField, GeoMesh>&,
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > eigenVectors
(
    const GeometricField<tensor, PatchField, GeoMesh>&
);

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > eigenVectors
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GeometricTensorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
