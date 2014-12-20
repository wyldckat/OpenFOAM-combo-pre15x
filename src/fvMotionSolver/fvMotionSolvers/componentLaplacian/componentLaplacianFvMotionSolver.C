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

#include "componentLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(componentLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMotionSolver,
        componentLaplacianFvMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::componentLaplacianFvMotionSolver::componentLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    Istream& msData
)
:
    fvMotionSolver(mesh),
    cmptName_(msData),
    cmpt_(0),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU" + cmptName_,
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
            "pointMotionU" + cmptName_,
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
{
    if (cmptName_ == "x")
    {
        cmpt_ = vector::X;
    }
    else if (cmptName_ == "y")
    {
        cmpt_ = vector::Y;
    }
    else if (cmptName_ == "z")
    {
        cmpt_ = vector::Z;
    }
    else
    {
        FatalErrorIn
        (
            "componentLaplacianFvMotionSolver::componentLaplacianFvMotionSolver"
            "(const polyMesh& mesh, Istream& msData)"
        )   << "Given component name " << cmptName_ << " should be x, y or z"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::componentLaplacianFvMotionSolver::~componentLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::componentLaplacianFvMotionSolver::curPoints() const
{
    // The points have moved so before interpolation update
    // volPointInterpolation accordingly
    vpi_.movePoints();

    vpi_.interpolate(cellMotionU_, pointMotionU_);

    tmp<pointField> tcurPoints(new pointField(fvMesh_.points()));

    tcurPoints().replace
    (
        cmpt_,
        tcurPoints().component(cmpt_)
      + fvMesh_.time().deltaT().value()*pointMotionU_.internalField()
    );

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::componentLaplacianFvMotionSolver::solve()
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
