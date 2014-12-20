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

#include "tetDecompositionMotionSolver.H"
#include "tetFec.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetDecompositionMotionSolver, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetDecompositionMotionSolver::tetDecompositionMotionSolver
(
    const polyMesh& mesh
)
:
    motionSolver(mesh),
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
    totDisplacementPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetDecompositionMotionSolver::~tetDecompositionMotionSolver()
{
    deleteDemandDrivenData(totDisplacementPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::tetDecompositionMotionSolver::curPoints() const
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

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::tetDecompositionMotionSolver::updateMesh()
{
    notImplemented("tetDecompositionMotionSolver::updateMesh()");
}


void Foam::tetDecompositionMotionSolver::updateTetTopology
(
    const tetPolyMeshMapperFaceDecomp& mapper
)
{
    tetMesh_.updateMesh(mapper);
    motionSolver::updateMesh();
}


Foam::tmp<Foam::elementScalarField>
Foam::tetDecompositionMotionSolver::distortionEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "distortionEnergy",
                tetMesh().time().timeName(),
                mesh(),
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
Foam::tetDecompositionMotionSolver::deformationEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "deforamationEnergy",
                tetMesh().time().timeName(),
                mesh(),
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
Foam::tetDecompositionMotionSolver::totDistortionEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "totDeformationEnergy",
                tetMesh().time().timeName(),
                mesh(),
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
            "tetDecompositionMotionSolver::totDeformationEnergy()"
	)   << "Total displacement field is not stored "
	    << "in tetDecompositionMotionSolver object." << endl
	    << exit(FatalError);
    }

    elementTensorField gradU = tetFec::elementGrad(totDisplacement());

    Ud = (1.0/2.0)*((gradU&&gradU) + (gradU&&gradU.T()))
       - (1.0/3.0)*tr(gradU)*tr(gradU);

    return tUd;
}



Foam::tmp<Foam::elementScalarField>
Foam::tetDecompositionMotionSolver::totDeformationEnergy() const
{
    tmp<elementScalarField> tUd
    (
        new elementScalarField
        (
            IOobject
            (
                "totDistortionEnergy",
                tetMesh().time().timeName(),
                mesh(),
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
            "tetDecompositionMotionSolver::totDistortionEnergy()"
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
