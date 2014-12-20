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
    tensorFieldField

Description
    Specialisation of FieldField<T> for tensor.

SourceFiles
    tensorFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef tensorFieldField_H
#define tensorFieldField_H

#include "FieldField.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Field>
void hdual
(
    FieldField<Field, vector>&,
    const FieldField<Field, tensor>&
);

template<template<class> class Field>
tmp<FieldField<Field, vector> > operator*
(
    const tmp<FieldField<Field, tensor> >&
);

template<template<class> class Field>
void hdual
(
    FieldField<Field, tensor>&,
    const FieldField<Field, vector>&
);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > operator*
(
    const tmp<FieldField<Field, vector> >&
);

template<template<class> class Field>
void diag(FieldField<Field, vector>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, vector> > diag(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void tr(FieldField<Field, scalar>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > tr(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void dev(FieldField<Field, tensor>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > dev(const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > dev(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void dev2(FieldField<Field, tensor>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > dev2(const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > dev2(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void det(FieldField<Field, scalar>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, scalar> > det(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void inv(FieldField<Field, tensor>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > inv(const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > inv(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void hinv(FieldField<Field, tensor>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > hinv(const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > hinv(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void symm(FieldField<Field, tensor>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > symm(const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > symm(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void skew(FieldField<Field, tensor>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > skew(const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > skew(const tmp<FieldField<Field, tensor> >&);

template<template<class> class Field>
void eigenValues(FieldField<Field, vector>&, const FieldField<Field, tensor>&);

template<template<class> class Field>
tmp<FieldField<Field, vector> > eigenValues
(
    const tmp<FieldField<Field, tensor> >&
);

template<template<class> class Field>
void eigenVectors
(
    FieldField<Field, tensor>&,
    const FieldField<Field, tensor>&
);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > eigenVectors
(
    const FieldField<Field, tensor>&
);

template<template<class> class Field>
tmp<FieldField<Field, tensor> > eigenVectors
(
    const tmp<FieldField<Field, tensor> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tensorFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
