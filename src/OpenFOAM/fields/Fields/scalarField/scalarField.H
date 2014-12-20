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
    scalarField

Description
    Specialisation of Field<T> for scalar.

SourceFiles
    scalarField.C

\*---------------------------------------------------------------------------*/

#ifndef scalarField_H
#define scalarField_H

#include "Field.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<scalar> scalarField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<scalarField> scalarField::component(const direction) const;

template<>
void component
(
    scalarField& sf,
    const UList<scalar>& f,
    const direction
);

template<>
void scalarField::replace(const direction, const UList<scalar>& sf);

void stabilise(scalarField& Res, const UList<scalar>& sf, const scalar s);
tmp<scalarField> stabilise(const UList<scalar>&, const scalar s);
tmp<scalarField> stabilise(const tmp<scalarField>&, const scalar s);

void divide(scalarField&, const scalar, const UList<scalar>&);
tmp<scalarField> operator/(const scalar, const UList<scalar>&);
tmp<scalarField> operator/(const scalar, const tmp<scalarField>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void pow(scalarField&, const UList<scalar>&, const UList<scalar>&);
tmp<scalarField> pow(const UList<scalar>&, const UList<scalar>&);
tmp<scalarField> pow(const tmp<scalarField>&, const UList<scalar>&);
tmp<scalarField> pow(const UList<scalar>&, const tmp<scalarField>&);
tmp<scalarField> pow(const tmp<scalarField>&, const tmp<scalarField>&);

void pow(scalarField&, const UList<scalar>&, const scalar&);
tmp<scalarField> pow(const UList<scalar>&, const scalar&);
tmp<scalarField> pow(const tmp<scalarField>&, const scalar&);

void pow(scalarField&, const scalar&, const UList<scalar>&);
tmp<scalarField> pow(const scalar&, const UList<scalar>&);
tmp<scalarField> pow(const scalar&, const tmp<scalarField>&);

void atan2(scalarField&, const UList<scalar>&, const UList<scalar>&);
tmp<scalarField> atan2(const UList<scalar>&, const UList<scalar>&);
tmp<scalarField> atan2(const tmp<scalarField>&, const UList<scalar>&);
tmp<scalarField> atan2(const UList<scalar>&, const tmp<scalarField>&);
tmp<scalarField> atan2(const tmp<scalarField>&, const tmp<scalarField>&);


#define transFunc(func)                                \
void func(scalarField& Res, const UList<scalar>& sf);  \
tmp<scalarField> func(const UList<scalar>&);           \
tmp<scalarField> func(const tmp<scalarField>&);

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


#define transFunc(func)                                             \
void func(scalarField& Res, const int n, const UList<scalar>& sf);  \
tmp<scalarField> func(const int n, const UList<scalar>&);           \
tmp<scalarField> func(const int n, const tmp<scalarField>&);

transFunc(jn)
transFunc(yn)

#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
