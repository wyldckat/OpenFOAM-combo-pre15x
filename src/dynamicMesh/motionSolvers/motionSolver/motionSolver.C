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
#include "twoDPointCorrector.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSolver, 0);
    defineRunTimeSelectionTable(motionSolver, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    totDisplacementPtr_(NULL),
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
    word solverTypeName = 
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
	).lookup("solver");

    Info << "Selection motion solver: " << solverTypeName << endl;

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
    deleteDemandDrivenData(totDisplacementPtr_);
    deleteDemandDrivenData(correct2DPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


void Foam::motionSolver::updateTopology()
{
    tetMesh_.updateTopology();
    if (correct2DPtr_) correct2DPtr_->updateTopology();
}


Foam::tmp<Foam::elementScalarField>
Foam::motionSolver::distortionEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "distortionEnergy",
                tetMesh().time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    elementTensorField gradU = 
        tetFec::elementGrad(motionU_)*tetMesh().time().deltaT();
 
    Ud = (1.0/2.0)*((gradU && gradU) + (gradU && gradU.T()))
       - (1.0/3.0)*tr(gradU)*tr(gradU);
    
    return tUd;
}


Foam::tmp<Foam::elementScalarField>
Foam::motionSolver::deformationEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "deforamationEnergy",
                tetMesh().time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    elementTensorField gradU = 
	tetFec::elementGrad(motionU_)*tetMesh().time().deltaT();

    scalar nu = 0.25;

    Ud = (1.0/2.0)*((gradU&&gradU) + (gradU&&gradU.T()));
       + (nu/(1.0+2.0*nu))*tr(gradU)*tr(gradU);
    
    return tUd;
}



Foam::tmp<Foam::elementScalarField>
Foam::motionSolver::totDistortionEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "totDeformationEnergy",
                tetMesh().time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    if(!needTotDisplacement())
    {
        FatalErrorIn
        (
            "motionSolver::totDeformationEnergy()"
	)   << "Total displacement field is not stored "
	    << "in motionSolver object." << endl
	    << exit(FatalError);
    }

    elementTensorField gradU = tetFec::elementGrad(totDisplacement());

    Ud = (1.0/2.0)*((gradU&&gradU) + (gradU&&gradU.T()))
       - (1.0/3.0)*tr(gradU)*tr(gradU);

    return tUd;
}



Foam::tmp<Foam::elementScalarField>
Foam::motionSolver::totDeformationEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "totDistortionEnergy",
                tetMesh().time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tetMesh(),
            dimensionedScalar("0.0", dimless, 0.0)
        )
    );

    elementScalarField& Ud = tUd();

    if(!needTotDisplacement())
    {
        FatalErrorIn
        (
            "motionSolver::totDistortionEnergy()"
	)   << "Total displacement field is not stored." << endl
	    << exit(FatalError);
    }

    elementTensorField gradU = tetFec::elementGrad(totDisplacement());

    scalar nu = 0.25;

    Ud = (1.0/2.0)*((gradU && gradU) + (gradU && gradU.T()));
       + (nu/(1.0+2.0*nu))*tr(gradU)*tr(gradU);
    
    return tUd;
}


// ************************************************************************* //
