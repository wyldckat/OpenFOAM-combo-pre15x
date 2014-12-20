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
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based pseudoSolid motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

\*---------------------------------------------------------------------------*/

#include "pseudoSolidMotionSolver.H"
#include "motionDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "twoDPointCorrector.H"
#include "tetFem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudoSolidMotionSolver, 0);
    addToRunTimeSelectionTable(motionSolver, pseudoSolidMotionSolver, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pseudoSolidMotionSolver::pseudoSolidMotionSolver(const polyMesh& mesh)
:
    laplaceMotionSolver(mesh)
{
    const dictionary& pseudoSolidDic = subDict("pseudoSolid");

    nu_ = readScalar(pseudoSolidDic.lookup("poissonsRatio"));

    nCorrectors_ =  readInt(pseudoSolidDic.lookup("nCorrectors"));

    convergenceTolerance_ = 
	readScalar(pseudoSolidDic.lookup("convergenceTolerance"));
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudoSolidMotionSolver::~pseudoSolidMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pseudoSolidMotionSolver::solve()
{
    // Solve for mesh motion

    if (!frozen_ && !firstMotion_)
    {
        Info << "Correct mesh motion diffusion field." << endl;

        diffusion()->correct();
    }

    int iCorr = 0;
    scalar initialResidual = 0;

    do
    {
        Info << "Correction: " << ++iCorr << endl;

        tetFemVectorMatrix motionEqn
        (
            tetFem::laplacian(diffusion()->gamma(), motionU())
          + tetFem::laplacianTranspose(diffusion()->gamma(), motionU())
          + tetFem::laplacianTrace
            (
                (2*nu_/(1+2*nu_))*diffusion()->gamma(), 
                motionU()
            )
        );

        initialResidual = motionEqn.solve().initialResidual();

        Info << "Initial residual: " << initialResidual << endl;
    }
    while (initialResidual > convergenceTolerance_ && iCorr < nCorrectors_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
