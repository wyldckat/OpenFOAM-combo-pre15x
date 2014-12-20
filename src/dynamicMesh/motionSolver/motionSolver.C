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
    Virtual base class for mesh motion solver.

\*---------------------------------------------------------------------------*/

#include "motionSolver.H"
#include "Time.H"
#include "polyMesh.H"
#include "twoDPointCorrector.H"

#include <dlfcn.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSolver::motionSolver(const polyMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
	    )
    ),
    mesh_(mesh),
    correct2DPtr_(NULL)
{
    // Set 2-D motion correction
    Switch twoDMotion(lookup("twoDMotion"));

    if (twoDMotion)
    {
        correct2DPtr_ = new twoDPointCorrector(mesh);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::motionSolver> Foam::motionSolver::New
(
    const polyMesh& mesh
)
{
	IOdictionary solverDict
	(
	    IOobject
	    (
            "dynamicMeshDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
	    )
	);

    word solverTypeName(solverDict.lookup("solver"));

    Info << "Selecting motion solver: " << solverTypeName << endl;

    if (solverDict.found("solverLib"))
    {
        string solverLibName = solverDict.lookup("solverLib");

        if (solverLibName.size())
        {
            label nSolvers = 0;

            if (dictionaryConstructorTablePtr_)
            {
                nSolvers = dictionaryConstructorTablePtr_->size();
            }

            void* solverLibPtr =
                dlopen(solverLibName.c_str(), RTLD_LAZY|RTLD_GLOBAL);

            if (!solverLibPtr)
            {
                FatalErrorIn
                (
                    "motionSolver::New(const polyMesh& mesh)"
                )   << "could not load " << dlerror()
                    << exit(FatalError);
            }
            else if 
            (
               !dictionaryConstructorTablePtr_
             || dictionaryConstructorTablePtr_->size() <= nSolvers
            )
            {
                WarningIn
                (
                    "motionSolver::New(const polyMesh& mesh)"
                )   << "library " << solverLibName
                    << " did not introduce any new solvers"
                    << endl << endl;
            }
        }
    }

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "motionSolver::New(const polyMesh& mesh)"
        )   << "solver table is empty"
            << exit(FatalError);
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "motionSolver::New(const polyMesh& mesh)"
        )   << "Unknown solver type " << solverTypeName
            << endl << endl
            << "Valid solver types are: " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<motionSolver>(cstrIter()(mesh));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSolver::~motionSolver()
{
    deleteDemandDrivenData(correct2DPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::motionSolver::newPoints()
{
    solve();
    return curPoints();
}


void Foam::motionSolver::twoDCorrectPoints(pointField& p) const
{
    // Correct for 2-D motion
    if (correct2DPtr_)
    {
        correct2DPtr_->correctPoints(p);
    }
}


void Foam::motionSolver::updateMesh()
{
    if (correct2DPtr_) correct2DPtr_->updateMesh();
}


// ************************************************************************* //
