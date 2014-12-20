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
    tensorField

Description
    Specialisation of Field<T> for tensor.

SourceFiles
    tensorField.C

\*---------------------------------------------------------------------------*/

#ifndef tensorField_H
#define tensorField_H

#include "scalarField.H"
#include "vectorField.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<tensor> tensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void hdual(vectorField&, const UList<tensor>&);
tmp<vectorField> operator*(const tmp<tensorField>&);

void hdual(tensorField&, const UList<vector>&);
tmp<tensorField> operator*(const tmp<vectorField>&);

void diag(vectorField&, const UList<tensor>&);
tmp<vectorField> diag(const tmp<tensorField>&);

void tr(scalarField&, const UList<tensor>&);
tmp<scalarField> tr(const tmp<tensorField>&);

void dev(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> dev(const UList<tensor>&);
tmp<tensorField> dev(const tmp<tensorField>&);

void dev2(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> dev2(const UList<tensor>&);
tmp<tensorField> dev2(const tmp<tensorField>&);

void det(scalarField&, const UList<tensor>&);
tmp<scalarField> det(const tmp<tensorField>&);

void inv(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> inv(const UList<tensor>&);
tmp<tensorField> inv(const tmp<tensorField>&);

void hinv(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> hinv(const UList<tensor>&);
tmp<tensorField> hinv(const tmp<tensorField>&);

void symm(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> symm(const UList<tensor>&);
tmp<tensorField> symm(const tmp<tensorField>&);

void skew(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> skew(const UList<tensor>&);
tmp<tensorField> skew(const tmp<tensorField>&);

void eigenValues(vectorField&, const UList<tensor>&);
tmp<vectorField> eigenValues(const tmp<tensorField>&);

void eigenVectors(Field<tensor>&, const UList<tensor>&);
tmp<tensorField> eigenVectors(const UList<tensor>&);
tmp<tensorField> eigenVectors(const tmp<tensorField>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
