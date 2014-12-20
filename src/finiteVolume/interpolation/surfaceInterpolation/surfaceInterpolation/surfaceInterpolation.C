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
    Cell to face interpolation scheme. Included in fvMesh.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "demandDrivenData.H"
#include "coupledFvPatch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceInterpolation, 0);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void surfaceInterpolation::clearOut()
{
    deleteDemandDrivenData(weightingFactors_);
    deleteDemandDrivenData(differenceFactors_);
    deleteDemandDrivenData(correctionVectors_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

surfaceInterpolation::surfaceInterpolation(const fvMesh& fvm)
:
    fvSchemes(static_cast<const objectRegistry&>(fvm)),
    fvSolution(static_cast<const objectRegistry&>(fvm)),
    mesh_(fvm),
    weightingFactors_(NULL),
    differenceFactors_(NULL),
    orthogonal_(false),
    correctionVectors_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

surfaceInterpolation::~surfaceInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceScalarField& surfaceInterpolation::weights() const
{
    if (!weightingFactors_)
    {
        makeWeights();
    }

    return (*weightingFactors_);
}


const surfaceScalarField& surfaceInterpolation::deltaCoeffs() const
{
    if (!differenceFactors_)
    {
        makeDeltaCoeffs();
    }

    return (*differenceFactors_);
}


bool surfaceInterpolation::orthogonal() const
{
    if (orthogonal_ == false && !correctionVectors_)
    {
        makeCorrectionVectors();
    }

    return orthogonal_;
}


const surfaceVectorField& surfaceInterpolation::correctionVectors() const
{
    if (orthogonal())
    {
        FatalErrorIn("surfaceInterpolation::correctionVectors()")
            << "cannot return correctionVectors; mesh is orthogonal"
            << abort(FatalError);
    }

    return (*correctionVectors_);
}


// Do what is neccessary if the mesh has moved
bool surfaceInterpolation::movePoints()
{
    deleteDemandDrivenData(weightingFactors_);
    deleteDemandDrivenData(differenceFactors_);

    orthogonal_ = false;
    deleteDemandDrivenData(correctionVectors_);

    return true;
}


void surfaceInterpolation::makeWeights() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeWeights() : "
            << "Constructing weighting factors for face interpolation"
            << endl;
    }


    weightingFactors_ = new surfaceScalarField
    (
        IOobject
        (
            "weightingFactors",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless
    );
    surfaceScalarField& weightingFactors = *weightingFactors_;


    // Set local references to mesh data
    const surfaceVectorField& faceCentres = mesh_.Cf();
    const volVectorField& cellCentres = mesh_.C();
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();
    const surfaceVectorField& areas = mesh_.Sf();
    const surfaceScalarField& magAreas = mesh_.magSf();


    // Setup temporary storage for cell center to face center distances
    scalarField ownerCellFaceDistances(areas.size());
    scalarField neighbourCellFaceDistances(areas.size());

    forAll(owner, faceI)
    {
        // Note: added mag in the dot-product in the numerator.  For
        // all valid meshes, the non-orthogonality will be less that
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  Meshes with this problem will
        // be caught in mesh checking.  
        ownerCellFaceDistances[faceI] =
        mag
        (
            areas[faceI] & (faceCentres[faceI] - cellCentres[owner[faceI]])
        )/magAreas[faceI];

        neighbourCellFaceDistances[faceI] =
        mag
        (
            areas[faceI] &
            (
                cellCentres[neighbour[faceI]] - faceCentres[faceI]
            )
        )/magAreas[faceI];
    }


    weightingFactors.internalField() =
        neighbourCellFaceDistances
       /(
            ownerCellFaceDistances + neighbourCellFaceDistances
        );

    forAll(mesh_.boundary(), patchI)
    {
        mesh_.boundary()[patchI].makeWeights
        (
            weightingFactors.boundaryField()[patchI]
        );
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeWeights() : "
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
}


void surfaceInterpolation::makeDeltaCoeffs() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    differenceFactors_ = new surfaceScalarField
    (
        IOobject
        (
            "differenceFactors_",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless/dimLength
    );
    surfaceScalarField& DeltaCoeffs = *differenceFactors_;


    // Set local references to mesh data
    const volVectorField& cellCentres = mesh_.C();
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();
    const surfaceVectorField& areas = mesh_.Sf();
    const surfaceScalarField& magAreas = mesh_.magSf();

    forAll(owner, faceI)
    {
        vector delta =
            cellCentres[neighbour[faceI]] - cellCentres[owner[faceI]];
        vector unitArea = areas[faceI]/magAreas[faceI];

        // Standard cell-centre distance form
        //DeltaCoeffs[faceI] = (unitArea & delta)/magSqr(delta);

        // Slightly under-relaxed form
        //DeltaCoeffs[faceI] = 1.0/mag(delta);

        // More under-relaxed form
        DeltaCoeffs[faceI] = 1.0/(mag(unitArea & delta) + VSMALL);
    }

    forAll(DeltaCoeffs.boundaryField(), patchI)
    {
        mesh_.boundary()[patchI].makeDeltaCoeffs
        (
            DeltaCoeffs.boundaryField()[patchI]
        );
    }
}


void surfaceInterpolation::makeCorrectionVectors() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeCorrectionVectors() : "
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    correctionVectors_ = new surfaceVectorField
    (
        IOobject
        (
            "correctionVectors",
            mesh_.pointsInstance(),
            mesh_
        ),
        mesh_,
        dimless
    );
    surfaceVectorField& CorrVecs = *correctionVectors_;

    // Set local references to mesh data
    const volVectorField& cellCentres = mesh_.C();
    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();
    const surfaceVectorField& areas = mesh_.Sf();
    const surfaceScalarField& magAreas = mesh_.magSf();
    const surfaceScalarField& DeltaCoeffs = deltaCoeffs();

    forAll(owner, faceI)
    {
        vector unitArea = areas[faceI]/magAreas[faceI];
        vector delta =
            cellCentres[neighbour[faceI]] - cellCentres[owner[faceI]];

        CorrVecs[faceI] = unitArea - delta*DeltaCoeffs[faceI];
    }

    // Boundary correction vectors set to zero for boundary patches
    // and calculated consistently with internal corrections for
    // coupled patches

    forAll(CorrVecs.boundaryField(), patchI)
    {
        fvPatchVectorField& patchCorrVecs = CorrVecs.boundaryField()[patchI];

        if (!patchCorrVecs.coupled())
        {
            patchCorrVecs = vector::zero;
        }
        else
        {
            const fvPatchScalarField& patchDeltaCoeffs
                = DeltaCoeffs.boundaryField()[patchI];

            const fvPatch& p = patchCorrVecs.patch();

            vectorField patchDeltas = mesh_.boundary()[patchI].delta();

            forAll(p, patchFaceI)
            {
                vector unitArea =
                    areas.boundaryField()[patchI][patchFaceI]
                    /magAreas.boundaryField()[patchI][patchFaceI];

                const vector& delta = patchDeltas[patchFaceI];

                patchCorrVecs[patchFaceI] =
                    unitArea - delta*patchDeltaCoeffs[patchFaceI];
            }
        }
    }

    scalar NonOrthogCoeff = 0.0;

    if (magAreas.size() > 0)
    {
        NonOrthogCoeff =
            asin
            (
                min
                (
                    (sum(magAreas*mag(CorrVecs))/sum(magAreas)).value(),
                    1.0
                )
            )*180.0/mathematicalConstant::pi;
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeCorrectionVectors() : "
            << "non-orthogonality coefficient = " << NonOrthogCoeff << " deg."
            << endl;
    }

    //NonOrthogCoeff = 0.0;

    if (NonOrthogCoeff < 0.1)
    {
        orthogonal_ = true;
        deleteDemandDrivenData(correctionVectors_);
    }
    else
    {
        orthogonal_ = false;
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeCorrectionVectors() : "
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }

    //mesh_.checkFaceDotProduct(true, NULL);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
