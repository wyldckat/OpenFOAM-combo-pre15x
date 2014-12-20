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
    Scalar specific part of the implementation of GeometricField.

\*---------------------------------------------------------------------------*/

#include "GeometricScalarField.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void stabilise
(
    GeometricField<scalar, PatchField, GeoMesh>& result,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    stabilise(result.internalField(), gsf.internalField(), ds.value());
    stabilise(result.boundaryField(), gsf.boundaryField(), ds.value());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > stabilise
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tRes
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "stabilise(" + gsf.name() + ',' + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            ds.dimensions() + gsf.dimensions()
        )
    );

    stabilise(tRes(), gsf, ds);

    return tRes;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > stabilise
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const dimensioned<scalar>& ds
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tRes
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "stabilise(" + gsf.name() + ',' + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf,
            ds.dimensions() + gsf.dimensions()
        )
    );

    stabilise(tRes(), tRes(), ds);

    return tRes;
}


template<template<class> class PatchField, class GeoMesh>
void divide
(
    GeometricField<scalar, PatchField, GeoMesh>& result,
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    divide(result.internalField(), ds.value(), gsf.internalField());
    divide(result.boundaryField(), ds.value(), gsf.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tRes
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                ds.name() + "/" + gsf.name(),
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            ds.dimensions()/gsf.dimensions()
        )
    );

    divide(tRes(), ds, gsf);

    return tRes;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const dimensioned<scalar>& ds,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tRes
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                ds.name() + "/" + gsf.name(),
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf,
            ds.dimensions()/gsf.dimensions()
        )
    );

    divide(tRes(), ds, tRes());

    return tRes;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const scalar& s,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    return dimensionedScalar(s)/gsf;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > operator/
(
    const scalar& s,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    return dimensionedScalar(s)/tgsf;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>& Pow,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    pow(Pow.internalField(), gsf1.internalField(), gsf2.internalField());
    pow(Pow.boundaryField(), gsf1.boundaryField(), gsf2.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "pow(" + gsf1.name() + ", " + gsf2.name() + ')',
                gsf1.instance(),
                gsf1.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf1.mesh(),
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", 1.0, gsf2.dimensions())
            )
        )
    );

    pow(tPow(), gsf1, gsf2);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf1,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1 = tgsf1();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "pow(" + gsf1.name() + ", " + gsf2.name() + ')',
                gsf1.instance(),
                gsf1.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf1,
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", 1.0, gsf2.dimensions())
            )
        )
    );

    pow(tPow(), tPow(), gsf2);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "pow(" + gsf1.name() + ", " + gsf2.name() + ')',
                gsf1.instance(),
                gsf1.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf2,
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", 1.0, gsf2.dimensions())
            )
        )
    );

    pow(tPow(), gsf1, tPow());

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf1,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf2
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1 = tgsf1();
    const GeometricField<scalar, PatchField, GeoMesh>& gsf2 = tgsf2();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "pow(" + gsf1.name() + ", " + gsf2.name() + ')',
                gsf1.instance(),
                gsf1.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf1,
            pow
            (
                gsf1.dimensions(),
                dimensionedScalar("1", 1.0, gsf2.dimensions())
            )
        )
    );

    pow(tPow(), tPow(), gsf2);

    tgsf2.clear();

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>& tPow,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensioned<scalar>& ds
)
{
    pow(tPow.internalField(), gsf.internalField(), ds.value());
    pow(tPow.boundaryField(), gsf.boundaryField(), ds.value());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const dimensionedScalar& ds
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "pow(" + gsf.name() + ", " + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            pow(gsf.dimensions(), ds)
        )
    );

    pow(tPow(), gsf, ds);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const dimensionedScalar& ds
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "pow(" + gsf.name() + ", " + ds.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf,
            pow(gsf.dimensions(), ds)
        )
    );

    pow(tPow(), tPow(), ds);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const GeometricField<scalar, PatchField, GeoMesh>& gsf,
    const scalar& s
)
{
    return pow(gsf, dimensionedScalar(s));
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf,
    const scalar& s
)
{
    return pow(tgsf, dimensionedScalar(s));
}


template<template<class> class PatchField, class GeoMesh>
void pow
(
    GeometricField<scalar, PatchField, GeoMesh>& tPow,
    const dimensioned<scalar>& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    pow(tPow.internalField(), ds.value(), gsf.internalField());
    pow(tPow.boundaryField(), ds.value(), gsf.boundaryField());
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const dimensionedScalar& ds,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "pow(" + ds.name() + ", " + gsf.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gsf.mesh(),
            pow(ds, gsf.dimensions())
        )
    );

    pow(tPow(), ds, gsf);

    return tPow;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const dimensionedScalar& ds,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();

    tmp<GeometricField<scalar, PatchField, GeoMesh> > tPow
    (
        GeometricField<scalar, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "pow(" + ds.name() + ", " + gsf.name() + ')',
                gsf.instance(),
                gsf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgsf,
            pow(ds, gsf.dimensions())
        )
    );

    pow(tPow(), ds, gsf);

    return tPow;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const scalar& s,
    const GeometricField<scalar, PatchField, GeoMesh>& gsf
)
{
    return pow(dimensionedScalar(s), gsf);
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > pow
(
    const scalar& s,
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf
)
{
    return pow(dimensionedScalar(s), tgsf);
}


#define transFunc(func)                                                     \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
void func                                                                   \
(                                                                           \
    GeometricField<scalar, PatchField, GeoMesh>& gsf,                       \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1                 \
)                                                                           \
{                                                                           \
    func(gsf.internalField(), gsf1.internalField());                        \
    func(gsf.boundaryField(), gsf1.boundaryField());                        \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf                  \
)                                                                           \
{                                                                           \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        new GeometricField<scalar, PatchField, GeoMesh>                     \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            gsf.mesh(),                                                     \
            func(gsf.dimensions())                                          \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), gsf);                                                     \
                                                                            \
    return tFunc;                                                           \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf           \
)                                                                           \
{                                                                           \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();        \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        GeometricField<scalar, PatchField, GeoMesh>::New                    \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            tgsf,                                                           \
            func(gsf.dimensions())                                          \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), tFunc());                                                 \
                                                                            \
    return tFunc;                                                           \
}

transFunc(pow3)
transFunc(pow4)
transFunc(sqrt)
transFunc(sign)
transFunc(pos)
transFunc(neg)

#undef transFunc


#define transFunc(func)                                                     \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
void func                                                                   \
(                                                                           \
    GeometricField<scalar, PatchField, GeoMesh>& gsf,                       \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1                 \
)                                                                           \
{                                                                           \
    func(gsf.internalField(), gsf1.internalField());                        \
    func(gsf.boundaryField(), gsf1.boundaryField());                        \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf                  \
)                                                                           \
{                                                                           \
    if (!gsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (#func"(const GeometricField<scalar, PatchField, GeoMesh>& gsf)")   \
            << "gsf not dimensionless"                                      \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        new GeometricField<scalar, PatchField, GeoMesh>                     \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            gsf.mesh(),                                                     \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), gsf);                                                     \
                                                                            \
    return tFunc;                                                           \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf           \
)                                                                           \
{                                                                           \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();        \
                                                                            \
    if (!gsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (#func"(const tmp<GeometricField<scalar, PatchField, GeoMesh> >& gsf)")\
            << " : gsf not dimensionless"                                   \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        GeometricField<scalar, PatchField, GeoMesh>::New                    \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            tgsf,                                                           \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), tFunc());                                                 \
                                                                            \
    return tFunc;                                                           \
}

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


#define transFunc(func)                                                     \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
void func                                                                   \
(                                                                           \
    GeometricField<scalar, PatchField, GeoMesh>& gsf,                       \
    const int n,                                                            \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf1                 \
)                                                                           \
{                                                                           \
    func(gsf.internalField(), n, gsf1.internalField());                     \
    func(gsf.boundaryField(), n, gsf1.boundaryField());                     \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const int n,                                                            \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf                  \
)                                                                           \
{                                                                           \
    if (!gsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (                                                                   \
            #func"(const int n, "                                           \
            "const GeometricField<scalar, PatchField, GeoMesh>& gsf)"       \
        )   << "gsf not dimensionless"                                      \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        new GeometricField<scalar, PatchField, GeoMesh>                     \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            gsf.mesh(),                                                     \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), n, gsf);                                                  \
                                                                            \
    return tFunc;                                                           \
}                                                                           \
                                                                            \
template<template<class> class PatchField, class GeoMesh>                   \
tmp<GeometricField<scalar, PatchField, GeoMesh> > func                      \
(                                                                           \
    const int n,                                                            \
    const tmp<GeometricField<scalar, PatchField, GeoMesh> >& tgsf           \
)                                                                           \
{                                                                           \
    const GeometricField<scalar, PatchField, GeoMesh>& gsf = tgsf();        \
                                                                            \
    if (!gsf.dimensions().dimensionless())                                  \
    {                                                                       \
        FatalErrorIn                                                        \
        (                                                                   \
            #func"(const int n, "                                           \
            "const tmp<GeometricField<scalar, PatchField, GeoMesh> >& gsf)" \
        )   << " : gsf not dimensionless"                                   \
            << abort(FatalError);                                           \
    }                                                                       \
                                                                            \
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tFunc                 \
    (                                                                       \
        GeometricField<scalar, PatchField, GeoMesh>::New                    \
        (                                                                   \
            IOobject                                                        \
            (                                                               \
                #func "(" + gsf.name() + ')',                               \
                gsf.instance(),                                             \
                gsf.db(),                                                   \
                IOobject::NO_READ,                                          \
                IOobject::NO_WRITE                                          \
            ),                                                              \
            tgsf,                                                           \
            dimless                                                         \
        )                                                                   \
    );                                                                      \
                                                                            \
    func(tFunc(), n, tFunc());                                              \
                                                                            \
    return tFunc;                                                           \
}

transFunc(jn)
transFunc(yn)

#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
