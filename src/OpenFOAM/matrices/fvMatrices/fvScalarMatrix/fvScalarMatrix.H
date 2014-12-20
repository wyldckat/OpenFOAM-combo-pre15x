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
    fvScalarMatrix

Description

SourceFiles
    fvScalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvScalarMatrix_H
#define fvScalarMatrix_H

#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<>
void fvMatrix<scalar>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction,
    const scalar value
);

template<>
lduMatrix::solverPerformance fvMatrix<scalar>::solve(Istream&);

// Return the matrix residual
template<>
tmp<scalarField> fvMatrix<scalar>::residual() const;

// H operator
template<>
tmp<volScalarField> fvMatrix<scalar>::H() const;

// operator&
template<>
tmp<scalarField> fvMatrix<scalar>::operator&(const scalarField&) const;

template<>
tmp<scalarField> fvMatrix<scalar>::operator&(const tmp<scalarField>&) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
