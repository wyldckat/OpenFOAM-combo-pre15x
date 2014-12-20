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
     Finite-Volume scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "fvScalarMatrix.H"
#include "zeroGradientFvPatchFields.H"
#include "lduCoupledInterfacePtrsList.H"

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
)
{
    internalCoeffs_[patchi][facei] +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

    boundaryCoeffs_[patchi][facei] = value;
}


template<>
lduMatrix::solverPerformance fvMatrix<scalar>::solve(Istream& solverControls)
{
    if (debug)
    {
        Info<< "fvMatrix<scalar>::solve(Istream& solverControls) : "
               "solving fvMatrix<scalar>"
            << endl;
    }

    scalarField saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    scalarField totalSource = source_;
    addBoundarySource(totalSource, 0);

    lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

    forAll (interfaces, patchI)
    {
        interfaces[patchI] = &psi_.boundaryField()[patchI];
    }

    // Solver call
    lduMatrix::solverPerformance solverPerf = lduMatrix::solver::New
    (
        psi_.name(),
        psi_.internalField(),
        *this,
        totalSource,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces,
        0,                         // dummy solving component
        solverControls
    )->solve();

    diag() = saveDiag;

    psi_.correctBoundaryConditions();

    return solverPerf;
}


// Return the matrix residual
template<>
tmp<scalarField> fvMatrix<scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);

    lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

    forAll (interfaces, patchI)
    {
        interfaces[patchI] = &psi_.boundaryField()[patchI];
    }

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.internalField(),
            source_ - boundaryDiag*psi_.internalField(),
            boundaryCoeffs_,
            interfaces,
            0
        )
    );

    addBoundarySource(tres());

    return tres;
}


// H operator
template<>
tmp<volScalarField> fvMatrix<scalar>::H() const
{
    tmp<volScalarField> tHphi
    (
        new volScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& Hphi = tHphi();

    Hphi.internalField() = (lduMatrix::H(psi_.internalField()) + source_);
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
