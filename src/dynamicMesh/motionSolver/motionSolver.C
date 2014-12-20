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
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

\*---------------------------------------------------------------------------*/

#include "motionSolver.H"
#include "motionDiff.H"
#include "twoDPointCorrector.H"
#include "tetFem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::motionSolver::motionSolver(const polyMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "motionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    tetMesh_(mesh),
    motionU_
    (
        IOobject
        (
            "motionU",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        tetMesh_
    ),
    gammaPtr_(NULL),
    correct2DPtr_(NULL),
    firstMotion_(true)
{
    // Set 2-D motion correction

    Switch twoDMotion(lookup("twoDMotion"));

    if (twoDMotion)
    {
        correct2DPtr_ = new twoDPointCorrector(mesh);
    }

    Switch frozen(lookup("frozenDiffusion"));

    if (frozen)
    {
        gammaPtr_ =
            new elementScalarField(motionDiff::New(tetMesh_, *this)->gamma());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolver::~motionSolver()
{
    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(correct2DPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::elementScalarField& Foam::motionSolver::gamma()
{
    if (!gammaPtr_)
    {
        freezeGamma();
    }

    return *gammaPtr_;
}


bool Foam::motionSolver::freezeGamma()
{
    if (!gammaPtr_)
    {
        gammaPtr_ =
            new elementScalarField
            (
                motionDiff::New(tetMesh_, *this)->gamma()
            );

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::motionSolver::releaseGamma()
{
    if (gammaPtr_)
    {
        deleteDemandDrivenData(gammaPtr_);
        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::pointField> Foam::motionSolver::newPoints()
{
    solve();

    return curPoints();
}


Foam::tmp<Foam::pointField> Foam::motionSolver::curPoints() const
{
    // Process current point positions
    tmp<pointField> tcurPoints
    (
        new pointField
        (
            tetMesh_().points()
          + vectorField
            (
                vectorField::subField
                (
                    motionU_.internalField(),
                    tetMesh_().nPoints()
                )
            )*tetMesh_.time().deltaT().value()
        )
    );

    // Correct for 2-D motion
    if (correct2DPtr_)
    {
        Info << "Correcting 2-D mesh motion ";
        correct2DPtr_->correctPoints(tcurPoints());
        Info << " ...done" << endl;
    }

    return tcurPoints;
}


void Foam::motionSolver::solve()
{
    // Solve for mesh motion
    bool frozen = true;

    // Check if the diffusion field is frozen; if not, calculate it
    if (!gammaPtr_)
    {
        frozen = false;

        gammaPtr_ =
            new elementScalarField(motionDiff::New(tetMesh_, *this)->gamma());
    }

    tetFemVectorMatrix motionEqn
    (
        tetFem::laplacian
        (
            *gammaPtr_,
            motionU_
        )
    );

    if (firstMotion_)
    {
        firstMotion_ = false;

        // In the first solution, solve the motion twice to avoid relative
        // tolerance problem
        for (label i = 0; i < 2; i++)
        {
            motionEqn.solve();
        }
    }
    else
    {
        motionEqn.solve();
    }

    // If diffusion is not frozen, delete it
    if (!frozen)
    {
        deleteDemandDrivenData(gammaPtr_);
    }
}

bool Foam::motionSolver::read()
{
    if (regIOobject::read())
    {
        Switch twoDMotion(lookup("twoDMotion"));

        if (twoDMotion)
        {
            deleteDemandDrivenData(correct2DPtr_);
            correct2DPtr_ =
                new twoDPointCorrector(tetMesh_());
        }

        Switch frozen(lookup("frozenDiffusion"));

        if (frozen)
        {
            if (!gammaPtr_)
            {
                // Freeze diffusivity
                freezeGamma();
            }
        }
        else
        {
            // Release diffisivity (if it has been frozen)
            releaseGamma();
        }

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::motionSolver::updateTopology()
{
    tetMesh_.updateTopology();
    if (correct2DPtr_) correct2DPtr_->updateTopology();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
