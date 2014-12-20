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
    Tensor specific part of the implementation of GeometricField.

\*---------------------------------------------------------------------------*/

#include "GeometricTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void hdual
(
    GeometricField<vector, PatchField, GeoMesh>& Hdual,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    hdual(Hdual.internalField(), gtf.internalField());
    hdual(Hdual.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > operator*
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<vector, PatchField, GeoMesh> > tHdual
    (
        new GeometricField<vector, PatchField, GeoMesh>
        (
            IOobject
            (
                "*" + gtf.name(),
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    hdual(tHdual(), gtf);

    return tHdual;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > operator*
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    tmp<GeometricField<vector, PatchField, GeoMesh> > tHdual = *tgtf();
    tgtf.clear();
    return tHdual;
}


template<template<class> class PatchField, class GeoMesh>
void hdual
(
    GeometricField<tensor, PatchField, GeoMesh>& Hdual,
    const GeometricField<vector, PatchField, GeoMesh>& gtf
)
{
    hdual(Hdual.internalField(), gtf.internalField());
    hdual(Hdual.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > operator*
(
    const GeometricField<vector, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tHdual
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "*" + gtf.name(),
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    hdual(tHdual(), gtf);

    return tHdual;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > operator*
(
    const tmp<GeometricField<vector, PatchField, GeoMesh> >& tgtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tHdual = *tgtf();
    tgtf.clear();
    return tHdual;
}


// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<template<class> class PatchField, class GeoMesh>
void tr
(
    GeometricField<scalar, PatchField, GeoMesh>& Tr,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tr(Tr.internalField(), gtf.internalField());
    tr(Tr.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > tr
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tTr
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "tr(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    tr(tTr(), gtf);

    return tTr;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > tr
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tTr = tr(tgtf());
    tgtf.clear();
    return tTr;
}


template<template<class> class PatchField, class GeoMesh>
void dev
(
    GeometricField<tensor, PatchField, GeoMesh>& Dev,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    dev(Dev.internalField(), gtf.internalField());
    dev(Dev.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tDev
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "dev(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    dev(tDev(), gtf);

    return tDev;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tDev
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "dev(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            gtf.dimensions()
        )
    );

    dev(tDev(), tDev());

    return tDev;
}


template<template<class> class PatchField, class GeoMesh>
void dev2
(
    GeometricField<tensor, PatchField, GeoMesh>& Dev2,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    dev2(Dev2.internalField(), gtf.internalField());
    dev2(Dev2.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev2
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tDev2
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "dev2(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    dev2(tDev2(), gtf);

    return tDev2;
}


template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > dev2
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tDev2
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "dev2(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            gtf.dimensions()
        )
    );

    dev2(tDev2(), tDev2());

    return tDev2;
}


template<template<class> class PatchField, class GeoMesh>
void det
(
    GeometricField<scalar, PatchField, GeoMesh>& Det,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    det(Det.internalField(), gtf.internalField());
    det(Det.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > det
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tDet
    (
        new GeometricField<scalar, PatchField, GeoMesh>
        (
            IOobject
            (
                "det(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            pow(gtf.dimensions(), tensor::dim)
        )
    );

    det(tDet(), gtf);

    return tDet;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<scalar, PatchField, GeoMesh> > det
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    tmp<GeometricField<scalar, PatchField, GeoMesh> > tDet = det(tgtf());
    tgtf.clear();
    return tDet;
}


template<template<class> class PatchField, class GeoMesh>
void inv
(
    GeometricField<tensor, PatchField, GeoMesh>& Inv,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    inv(Inv.internalField(), gtf.internalField());
    inv(Inv.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > inv
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tInv
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "inv(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            dimless/gtf.dimensions()
        )
    );

    inv(tInv(), gtf);

    return tInv;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > inv
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tInv
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "inv(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            dimless/gtf.dimensions()
        )
    );

    inv(tInv(), tInv());

    return tInv;
}



template<template<class> class PatchField, class GeoMesh>
void hinv
(
    GeometricField<tensor, PatchField, GeoMesh>& Hinv,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    hinv(Hinv.internalField(), gtf.internalField());
    hinv(Hinv.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > hinv
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tHinv
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "hinv(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            dimless/gtf.dimensions()
        )
    );

    hinv(tHinv(), gtf);

    return tHinv;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > hinv
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tHinv
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "hinv(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            dimless/gtf.dimensions()
        )
    );

    hinv(tHinv(), tHinv());

    return tHinv;
}


template<template<class> class PatchField, class GeoMesh>
void symm
(
    GeometricField<tensor, PatchField, GeoMesh>& Symm,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    symm(Symm.internalField(), gtf.internalField());
    symm(Symm.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > symm
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tSymm
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "symm(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    symm(tSymm(), gtf);

    return tSymm;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > symm
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tSymm
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "symm(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            gtf.dimensions()
        )
    );

    symm(tSymm(), tSymm());

    return tSymm;
}


template<template<class> class PatchField, class GeoMesh>
void skew
(
    GeometricField<tensor, PatchField, GeoMesh>& Skew,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    skew(Skew.internalField(), gtf.internalField());
    skew(Skew.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > skew
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tSkew
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "skew(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            gtf.dimensions()
        )
    );

    skew(tSkew(), gtf);

    return tSkew;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > skew
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tSkew
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "skew(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            gtf.dimensions()
        )
    );

    skew(tSkew(), tSkew());

    return tSkew;
}


template<template<class> class PatchField, class GeoMesh>
void eigenValues
(
    GeometricField<vector, PatchField, GeoMesh>& EigenValues,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    eigenValues(EigenValues.internalField(), gtf.internalField());
    eigenValues(EigenValues.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > eigenValues
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<vector, PatchField, GeoMesh> > tEigenValues
    (
        new GeometricField<vector, PatchField, GeoMesh>
        (
            IOobject
            (
                "eigenValues(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            pow(gtf.dimensions(), tensor::dim)
        )
    );

    eigenValues(tEigenValues(), gtf);

    return tEigenValues;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<vector, PatchField, GeoMesh> > eigenValues
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    tmp<GeometricField<vector, PatchField, GeoMesh> > tEigenValues =
        eigenValues(tgtf());
    tgtf.clear();
    return tEigenValues;
}


template<template<class> class PatchField, class GeoMesh>
void eigenVectors
(
    GeometricField<tensor, PatchField, GeoMesh>& EigenVectors,
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    eigenVectors(EigenVectors.internalField(), gtf.internalField());
    eigenVectors(EigenVectors.boundaryField(), gtf.boundaryField());
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > eigenVectors
(
    const GeometricField<tensor, PatchField, GeoMesh>& gtf
)
{
    tmp<GeometricField<tensor, PatchField, GeoMesh> > tEigenVectors
    (
        new GeometricField<tensor, PatchField, GeoMesh>
        (
            IOobject
            (
                "eigenVectors(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gtf.mesh(),
            dimless
        )
    );

    eigenVectors(tEigenVectors(), gtf);

    return tEigenVectors;
}

template<template<class> class PatchField, class GeoMesh>
tmp<GeometricField<tensor, PatchField, GeoMesh> > eigenVectors
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& tgtf
)
{
    const GeometricField<tensor, PatchField, GeoMesh>& gtf = tgtf();

    tmp<GeometricField<tensor, PatchField, GeoMesh> > tEigenVectors
    (
        GeometricField<tensor, PatchField, GeoMesh>::New
        (
            IOobject
            (
                "eigenVectors(" + gtf.name() + ')',
                gtf.instance(),
                gtf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtf,
            dimless
        )
    );

    eigenVectors(tEigenVectors(), tEigenVectors());

    return tEigenVectors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
