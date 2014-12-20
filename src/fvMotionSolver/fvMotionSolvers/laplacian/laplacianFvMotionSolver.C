/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "laplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        laplacianFvMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laplacianFvMotionSolver::laplacianFvMotionSolver
(
    const polyMesh& mesh,
    Istream& msData
)
:
    fvMotionSolver(mesh),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fvMesh_
    ),
    pointMotionU_
    (
        IOobject
        (
            "pointMotionU",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vpi_.interpolate(cellMotionU_, cellMotionU_.boundaryField().types())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(*this, lookup("diffusivity"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laplacianFvMotionSolver::~laplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::laplacianFvMotionSolver::curPoints() const
{
    // The points have moved so before interpolation update 
    // volPointInterpolation accordingly
    vpi_.movePoints();

    vpi_.interpolate(cellMotionU_, pointMotionU_);

    tmp<pointField> tcurPoints
    (
        fvMesh_.points()
      + fvMesh_.time().deltaT().value()*pointMotionU_.internalField()
    );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::laplacianFvMotionSolver::solve()
{
    diffusivityPtr_->correct();
    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellMotionU_,
            "laplacian(diffusivity,cellMotionU)"
        )
    );
}


// ************************************************************************* //
