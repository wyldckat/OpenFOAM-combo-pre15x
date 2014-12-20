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
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresVectors::~leastSquaresVectors()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
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

    const fvMesh& mesh = mesh_;

    pVectorsPtr_ = new surfaceVectorField
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
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
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
    surfaceVectorField& lsN = *nVectorsPtr_;

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


    forAll(owner, facei)
    {
        tensor wdd = 1.0/magSqr(d[facei])*d[facei]*d[facei];

        dd[owner[facei]] += wdd;

        // Yes, it is += because both vectors have the "wrong" sign
        dd[neighbour[facei]] += wdd;
    }

    // Visit the boundaries. Coupled boundaries are taken into account
    // in the construction of d vectors.  
    forAll(d.boundaryField(), patchi)
    {
        const fvPatchVectorField& pd = d.boundaryField()[patchi];

        const labelList::subList faceCells = pd.patch().faceCells();

        forAll(pd, patchFacei)
        {
            label faceCelli = faceCells[patchFacei];

            dd[faceCelli] +=
                (1.0/magSqr(pd[patchFacei]))*pd[patchFacei]*pd[patchFacei];
        }
    }


    // Invert the dd tensor
    tensorField invDd = inv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, facei)
    {
        lsP[facei] =
            (1.0/magSqr(d[facei]))*(invDd[owner[facei]] & d[facei]);

        lsN[facei] =
            ((-1.0)/magSqr(d[facei]))*(invDd[neighbour[facei]] & d[facei]);
    }

    forAll(lsP.boundaryField(), patchi)
    {
        const fvPatchVectorField& pd = d.boundaryField()[patchi];

        fvPatchVectorField& patchLsP = lsP.boundaryField()[patchi];

        const fvPatch& p = patchLsP.patch();
        const labelList::subList FaceCells = p.faceCells();

        forAll(p, patchFacei)
        {
            patchLsP[patchFacei] =
                (1.0/magSqr(pd[patchFacei]))
               *(invDd[FaceCells[patchFacei]] & pd[patchFacei]);
        }
    }


    // For 3D meshes check the determinant of the dd tensor and switch to
    // Gauss if it is less than 2
    if (mesh.nD() == 3)
    {
        label nBadCells = 0;

        const cellList& cells = mesh.cells();
        const scalarField& V = mesh.V();
        const surfaceVectorField& Sf = mesh.Sf();
        const surfaceScalarField& w = mesh.weights();

        forAll (dd, celli)
        {
            if (det(dd[celli]) < 3)
            {
                nBadCells++;

                const cell& c = cells[celli];

                forAll(c, cellFacei)
                {
                    label facei = c[cellFacei];

                    if (mesh.isInternalFace(facei))
                    {
                        if (celli == owner[facei])
                        {
                            lsP[facei] = (1 - w[facei])*Sf[facei]/V[celli];
                        }
                        else
                        {
                            lsN[facei] = -w[facei]*Sf[facei]/V[celli];
                        }
                    }
                    else
                    {
                        label patchi = mesh.boundaryMesh().whichPatch(facei);

                        if (mesh.boundary()[patchi].size())
                        {
                            label patchFacei = 
                                facei - mesh.boundaryMesh()[patchi].start();
                        
                            lsP.boundaryField()[patchi][patchFacei] = 
                                Sf.boundaryField()[patchi][patchFacei]/V[celli];
                        }
                    }
                }
            }
        }

        if (debug)
        {
            InfoIn("leastSquaresVectors::makeLeastSquaresVectors()")
                << "number of bad cells switched to Gauss = " << nBadCells
                << endl;
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
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::leastSquaresVectors::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


bool Foam::leastSquaresVectors::movePoints()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);

    return true;
}


// ************************************************************************* //
