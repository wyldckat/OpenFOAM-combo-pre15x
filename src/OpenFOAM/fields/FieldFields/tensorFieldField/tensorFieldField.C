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

Description
    Specialisation of FieldField<T> for tensor.

\*---------------------------------------------------------------------------*/

#include "tensorFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

template<template<class> class Field>
void hdual(FieldField<Field, vector>& vf, const FieldField<Field, tensor>& tf)
{
    forAll(vf, i)
    {
        hdual(vf[i], tf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, vector> > operator*
(
    const tmp<FieldField<Field, tensor> >& tf
)
{
    tmp<FieldField<Field, vector> > hDual
    (
        FieldField<Field, vector>::NewCalculatedType(tf())
    );
    hdual(hDual(), tf);
    tf.clear();
    return hDual;
}


template<template<class> class Field>
void hdual(FieldField<Field, tensor>& vf, const FieldField<Field, vector>& tf)
{
    forAll(vf, i)
    {
        hdual(vf[i], tf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, tensor> > operator*
(
    const tmp<FieldField<Field, vector> >& tf
)
{
    tmp<FieldField<Field, tensor> > hDual
    (
        FieldField<Field, tensor>::NewCalculatedType(tf())
    );
    hdual(hDual(), tf);
    tf.clear();
    return hDual;
}


// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

#define tensorFunc(func)                                                      \
template<template<class> class Field>                                         \
void func(FieldField<Field, tensor>& tf, const FieldField<Field, tensor>& tf1)\
{                                                                             \
    forAll(tf, i)                                                             \
    {                                                                         \
        func(tf[i], tf1[i]);                                                  \
    }                                                                         \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, tensor> > func(const FieldField<Field, tensor>& tf)     \
{                                                                             \
    tmp<FieldField<Field, tensor> > result                                    \
    (                                                                         \
        FieldField<Field, tensor>::NewCalculatedType(tf)                      \
    );                                                                        \
    func(result(), tf);                                                       \
    return result;                                                            \
}                                                                             \
                                                                              \
template<template<class> class Field>                                         \
tmp<FieldField<Field, tensor> > func                                          \
(                                                                             \
    const tmp<FieldField<Field, tensor> >& tf                                 \
)                                                                             \
{                                                                             \
    tmp<FieldField<Field, tensor> > result(tf.ptr());                         \
    func(result(), result());                                                 \
    return result;                                                            \
}


template<template<class> class Field>
void diag(FieldField<Field, vector>& vf, const FieldField<Field, tensor>& tf)
{
    forAll(vf, i)
    {
        diag(vf[i], tf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, vector> > diag(const tmp<FieldField<Field, tensor> >& tf)
{
    tmp<FieldField<Field, vector> > result
    (
        FieldField<Field, vector>::NewCalculatedType(tf())
    );
    diag(result(), tf);
    tf.clear();
    return result;
}

template<template<class> class Field>
void tr(FieldField<Field, scalar>& sf, const FieldField<Field, tensor>& tf)
{
    forAll(sf, i)
    {
        tr(sf[i], tf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > tr(const tmp<FieldField<Field, tensor> >& tf)
{
    tmp<FieldField<Field, scalar> > result
    (
        FieldField<Field, scalar>::NewCalculatedType(tf())
    );
    tr(result(), tf);
    tf.clear();
    return result;
}


tensorFunc(dev)
tensorFunc(dev2)

template<template<class> class Field>
void det(FieldField<Field, scalar>& sf, const FieldField<Field, tensor>& tf)
{
    forAll(sf, i)
    {
        det(sf[i], tf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, scalar> > det(const tmp<FieldField<Field, tensor> >& tf)
{
    tmp<FieldField<Field, scalar> > result
    (
        FieldField<Field, scalar>::NewCalculatedType(tf())
    );
    det(result(), tf);
    tf.clear();
    return result;
}


template<template<class> class Field>
void inv(FieldField<Field, tensor>& tf, const FieldField<Field, tensor>& tf1)
{
    forAll(tf, i)
    {
        inv(tf[i], tf1[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, tensor> > inv(const FieldField<Field, tensor>& tf)
{
    tmp<FieldField<Field, tensor> > result
    (
        FieldField<Field, tensor>::NewCalculatedType(tf())
    );
    inv(result(), tf);
    return result;
}

template<template<class> class Field>
tmp<FieldField<Field, tensor> > inv(const tmp<FieldField<Field, tensor> >& tf)
{
    tmp<FieldField<Field, tensor> > result(tf.ptr());
    inv(result(), result());
    return result;
}


tensorFunc(hinv)
tensorFunc(symm)
tensorFunc(skew)
template<template<class> class Field>
void eigenValues
(
    FieldField<Field, vector>& vf,
    const FieldField<Field, tensor>& tf
)
{
    forAll(vf, i)
    {
        eigenValues(vf[i], tf[i]);
    }
}

template<template<class> class Field>
tmp<FieldField<Field, vector> > eigenValues
(
    const tmp<FieldField<Field, tensor> >& tf
)
{
    tmp<FieldField<Field, vector> > result
    (
        FieldField<Field, vector>::NewCalculatedType(tf())
    );
    eigenValues(result(), tf);
    tf.clear();
    return result;
}


tensorFunc(eigenVectors)

#undef tensorFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
