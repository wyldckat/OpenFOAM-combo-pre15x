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

\*---------------------------------------------------------------------------*/

#include "leastSquaresVectors.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::leastSquaresVectors(const fvMesh& mesh)
:
    MeshObject<leastSquaresVectors>(mesh),
    pVectors_(NULL),
    nVectors_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{
    deleteDemandDrivenData(pVectors_);
    deleteDemandDrivenData(nVectors_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// LeastSquaresVectors using inverse-distance-squared weighting which gives
// the same gradient as Gauss linear for orthogonal grids

void Foam::leastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    pVectors_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectors_;

    nVectors_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectors_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Build the d-vectors
    surfaceVectorField d = mesh_.Sf()/(mesh_.magSf()*mesh_.deltaCoeffs());

    if (!mesh_.orthogonal())
    {
        d -= mesh_.correctionVectors()/mesh_.deltaCoeffs();
    }

    // Set up temporary storage for the dd tensor (before inversion)
    tensorField dd(mesh_.nCells(), tensor::zero);

    forAll(owner, faceI)
    {
        tensor wdd = 1.0/magSqr(d[faceI])*d[faceI]*d[faceI];

        dd[owner[faceI]] += wdd;

        // Yes, it is += because both vectors have the "wrong" sign
        dd[neighbour[faceI]] += wdd;
    }

    // Visit the boundaries. Coupled boundaries are taken into account
    // in the construction of d vectors.  
    forAll(d.boundaryField(), patchI)
    {
        const fvPatchVectorField& pd = d.boundaryField()[patchI];

        const labelList::subList FaceCells = pd.patch().faceCells();

        forAll(pd, patchFaceI)
        {
            dd[FaceCells[patchFaceI]] +=
                (1.0/magSqr(pd[patchFaceI]))*pd[patchFaceI]*pd[patchFaceI];
        }
    }

    // Invert the dd tensor
    tensorField invDd = inv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, faceI)
    {
        lsP[faceI] =
            (1.0/magSqr(d[faceI]))*(invDd[owner[faceI]] & d[faceI]);

        lsN[faceI] =
            ((-1.0)/magSqr(d[faceI]))*(invDd[neighbour[faceI]] & d[faceI]);
    }

    forAll(lsP.boundaryField(), patchI)
    {
        const fvPatchVectorField& pd = d.boundaryField()[patchI];

        fvPatchVectorField& patchLsP = lsP.boundaryField()[patchI];

        const fvPatch& p = patchLsP.patch();
        const labelList::subList FaceCells = p.faceCells();

        forAll(p, patchFaceI)
        {
            patchLsP[patchFaceI] =
                (1.0/magSqr(pd[patchFaceI]))
               *(invDd[FaceCells[patchFaceI]] & pd[patchFaceI]);
        }
    }

    if (debug)
    {
        Info<< "leastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::pVectors() const
{
    if (!pVectors_)
    {
        makeLeastSquaresVectors();
    }

    return (*pVectors_);
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::nVectors() const
{
    if (!nVectors_)
    {
        makeLeastSquaresVectors();
    }

    return (*nVectors_);
}


bool Foam::leastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(pVectors_);
    deleteDemandDrivenData(nVectors_);

    return true;
}


// ************************************************************************* //
