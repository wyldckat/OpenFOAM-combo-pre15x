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
     Finite-Area scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "faScalarMatrix.H"
#include "zeroGradientFaPatchFields.H"
#include "lduCoupledInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<>
void faMatrix<scalar>::setComponentReference
(
    const label patchI,
    const label edgeI,
    const direction,
    const scalar value
)
{
    const labelList::subList faceLabels =
        psi_.mesh().boundary()[patchI].edgeFaces();

    internalCoeffs_[patchI][edgeI] +=
        diag()[faceLabels[edgeI]];

    boundaryCoeffs_[patchI][edgeI] = value;
}


template<>
lduMatrix::solverPerformance faMatrix<scalar>::solve(Istream& solverControls)
{
#   ifdef DEBUGfaMatrix
    if (debug)
    {
        Info<< "faMatrix<scalar>::solve(Istream& solverControls) : "
               "solving faMatrix<scalar>"
            << endl;
    }
#   endif

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

    solverPerf.print();

    diag() = saveDiag;

    psi_.correctBoundaryConditions();

    return solverPerf;
}


// Return the matrix residual
template<>
tmp<scalarField> faMatrix<scalar>::residual() const
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
tmp<areaScalarField> faMatrix<scalar>::H() const
{
    tmp<areaScalarField> tHphi
    (
        new areaScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField Hphi = tHphi();

    Hphi.internalField() = (lduMatrix::H(psi_.internalField()) + source_);
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


// operator&
template<>
tmp<scalarField> faMatrix<scalar>::operator&(const scalarField& psi) const
{
    tmp<scalarField > tApsi(new scalarField(psi.size()));
    scalarField& Apsi = tApsi();

    // There is no boundaries on psi, so no interface coupling
    // 
    lduMatrix::Amul
    (
        Apsi,
        psi,
        FieldField<Field, scalar>(0),
        lduCoupledInterfacePtrsList(0),
        0
    );

    Apsi -= source_;

    scalarField boundaryDiag(psi.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);
    Apsi += boundaryDiag*psi;

    scalarField boundarySource(psi.size(), 0.0);
    addBoundarySource(boundarySource);

    Apsi -= boundarySource;

    return tApsi;
}

template<>
tmp<scalarField> faMatrix<scalar>::operator&(const tmp<scalarField>& tpsi) const
{
    tmp<scalarField > tHpsi = operator&(tpsi());
    tpsi.clear();
    return tHpsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
