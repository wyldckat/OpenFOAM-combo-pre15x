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
    Cell to face interpolation scheme. Included in fvMesh.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "demandDrivenData.H"
#include "coupledFvPatch.H"
#include "physicalConstants.H"

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
    deleteDemandDrivenData(skewCorrectionVectors_);
    deleteDemandDrivenData(leastSquarePvectors_);
    deleteDemandDrivenData(leastSquareNvectors_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

surfaceInterpolation::surfaceInterpolation(const fvMesh& fvm)
:
    fvSchemes((const objectRegistry&)fvm),
    fvSolution((const objectRegistry&)fvm),
    mesh_(fvm),
    weightingFactors_(NULL),
    differenceFactors_(NULL),
    orthogonal_(false),
    correctionVectors_(NULL),
    skew_(true),
    skewCorrectionVectors_(NULL),
    leastSquarePvectors_(NULL),
    leastSquareNvectors_(NULL)
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


bool surfaceInterpolation::skew() const
{
    if (skew_ == true && !skewCorrectionVectors_)
    {
        makeSkewCorrectionVectors();
    }

    return skew_;
}


const surfaceVectorField& surfaceInterpolation::skewCorrectionVectors() const
{
    if (!skew())
    {
        FatalErrorIn("surfaceInterpolation::skewCorrectionVectors()")
            << "cannot return correctionVectors; mesh is now skewed"
            << abort(FatalError);
    }

    return (*skewCorrectionVectors_);
}


const surfaceVectorField& surfaceInterpolation::leastSquarePvectors() const
{
    if (!leastSquarePvectors_)
    {
        makeLeastSquareVectors();
    }

    return (*leastSquarePvectors_);
}


const surfaceVectorField& surfaceInterpolation::leastSquareNvectors() const
{
    if (!leastSquareNvectors_)
    {
        makeLeastSquareVectors();
    }

    return (*leastSquareNvectors_);
}


// Do what is neccessary if the mesh has moved
bool surfaceInterpolation::movePoints()
{
    deleteDemandDrivenData(weightingFactors_);
    deleteDemandDrivenData(differenceFactors_);

    orthogonal_ = false;
    deleteDemandDrivenData(correctionVectors_);

    skew_ = true;
    deleteDemandDrivenData(skewCorrectionVectors_);

    deleteDemandDrivenData(leastSquarePvectors_);
    deleteDemandDrivenData(leastSquareNvectors_);

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
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimless
    );
    surfaceScalarField& weightingFactors = *weightingFactors_;


    // Set local references to mesh data
    const surfaceVectorField& faceCentres = mesh().Cf();
    const volVectorField& cellCentres = mesh().C();
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
    const surfaceVectorField& areas = mesh().Sf();
    const surfaceScalarField& magAreas = mesh().magSf();


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
#          ifdef BAD_MESH_STABILISATION
         + VSMALL
#          endif
        );

    forAll(mesh().boundary(), patchI)
    {
        mesh().boundary()[patchI].makeWeights
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
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimless/dimLength
    );
    surfaceScalarField& DeltaCoeffs = *differenceFactors_;


    // Set local references to mesh data
    const volVectorField& cellCentres = mesh().C();
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
    const surfaceVectorField& areas = mesh().Sf();
    const surfaceScalarField& magAreas = mesh().magSf();

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
        DeltaCoeffs[faceI] = 
            1.0
           /(
               mag(unitArea & delta)
#              ifdef BAD_MESH_STABILISATION
             + VSMALL
#              endif
           );
    }

    forAll(DeltaCoeffs.boundaryField(), patchI)
    {
        mesh().boundary()[patchI].makeDeltaCoeffs
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
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimless
    );
    surfaceVectorField& CorrVecs = *correctionVectors_;

    // Set local references to mesh data
    const volVectorField& cellCentres = mesh().C();
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();
    const surfaceVectorField& areas = mesh().Sf();
    const surfaceScalarField& magAreas = mesh().magSf();

#   ifndef NEW_NON_ORTH
    const surfaceScalarField& DeltaCoeffs = deltaCoeffs();
#   endif

    forAll(owner, faceI)
    {
        vector unitArea = areas[faceI]/magAreas[faceI];
        vector delta =
            cellCentres[neighbour[faceI]] - cellCentres[owner[faceI]];

#       ifndef NEW_NON_ORTH
        CorrVecs[faceI] = unitArea - delta*DeltaCoeffs[faceI];
#       else
        CorrVecs[faceI] =
            unitArea
          - delta*
            (
                (unitArea & delta)
               /(
                   magSqr(delta)
#                  ifdef BAD_MESH_STABILISATION
                 + VSMALL
#                  endif
               )
            );
#       endif
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
#           ifndef NEW_NON_ORTH
            const fvPatchScalarField& patchDeltaCoeffs
                = DeltaCoeffs.boundaryField()[patchI];
#           endif

            const fvPatch& p = patchCorrVecs.patchMesh();

            vectorField patchDeltas = mesh().boundary()[patchI].delta();

            forAll(p, patchFaceI)
            {
                vector unitArea =
                    areas.boundaryField()[patchI][patchFaceI]
                    /magAreas.boundaryField()[patchI][patchFaceI];

                const vector& delta = patchDeltas[patchFaceI];

#               ifndef NEW_NON_ORTH
                patchCorrVecs[patchFaceI] =
                    unitArea - delta*patchDeltaCoeffs[patchFaceI];
#               else
                patchCorrVecs[patchFaceI] =
                    unitArea
                  - delta*
                    (
                        (unitArea & delta)
                       /(
                            magSqr(delta)
#                           ifdef BAD_MESH_STABILISATION
                          + VSMALL
#                           endif
                       )
                    );
#               endif
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
            )*180.0/physicalConstant::pi;
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
}


void surfaceInterpolation::makeSkewCorrectionVectors() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeSkewCorrectionVectors() : "
            << "Constructing skew correction vectors"
            << endl;
    }

    skewCorrectionVectors_ = new surfaceVectorField
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimless
    );
    surfaceVectorField& SkewCorrVecs = *skewCorrectionVectors_;

    // Set local references to mesh data
    const volVectorField& C = mesh().C();
    const surfaceVectorField& Cf = mesh().Cf();
    const surfaceVectorField& Sf = mesh().Sf();

    const unallocLabelList& owner = mesh().owner();

    // Build the d-vectors
    surfaceVectorField d = Sf/(mesh().magSf()*deltaCoeffs());

    if (!orthogonal())
    {
        d -= correctionVectors()/deltaCoeffs();
    }

    forAll(owner, faceI)
    {
        vector Cpf = Cf[faceI] - C[owner[faceI]];

        SkewCorrVecs[faceI] =
            Cpf - ((Sf[faceI] & Cpf)/(Sf[faceI] & d[faceI]))*d[faceI];
    }


    forAll(SkewCorrVecs.boundaryField(), patchI)
    {
        fvPatchVectorField& patchSkewCorrVecs =
            SkewCorrVecs.boundaryField()[patchI];

        if (!patchSkewCorrVecs.coupled())
        {
            patchSkewCorrVecs = vector::zero;
        }
        else
        {
            const fvPatch& p = patchSkewCorrVecs.patchMesh();
            const labelList::subList FaceCells = p.faceCells();
            const vectorField& patchFaceCentres = Cf.boundaryField()[patchI];
            const vectorField& patchSf = Sf.boundaryField()[patchI];
            const vectorField& patchD = d.boundaryField()[patchI];

            forAll(p, patchFaceI)
            {
                vector Cpf =
                    patchFaceCentres[patchFaceI] - C[FaceCells[patchFaceI]];

                patchSkewCorrVecs[patchFaceI] =
                    Cpf
                  - (
                        (patchSf[patchFaceI] & Cpf)/
                        (patchSf[patchFaceI] & patchD[patchFaceI])
                    )*patchD[patchFaceI];
            }
        }
    }

    scalar skewCoeff = 0.0;

    if (Sf.internalField().size() > 0)
    {
        skewCoeff = max(mag(SkewCorrVecs)/mag(d)).value();
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeSkewCorrectionVectors() : "
            << "skew coefficient = " << skewCoeff << endl;
    }

    //skewCoeff = 0.0;

    if (skewCoeff < 1e-5)
    {
        skew_ = false;
        deleteDemandDrivenData(skewCorrectionVectors_);
    }
    else
    {
        skew_ = true;
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeSkewCorrectionVectors() : "
            << "Finished constructing skew correction vectors"
            << endl;
    }
}

/* This version using inverse-distance-squared weighting doesn't work
   well for high aspect ratio cells and has been replaced by the unweighted
   version below.

void surfaceInterpolation::makeLeastSquareVectors() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeLeastSquareVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    leastSquarePvectors_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquareP",
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *leastSquarePvectors_;

    leastSquareNvectors_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquareN",
            mesh().PointsInstance(),
            mesh()
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *leastSquareNvectors_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // Build the d-vectors
    surfaceVectorField d = mesh().Sf()/(mesh().magSf()*deltaCoeffs());

    if (!orthogonal())
    {
        d -= correctionVectors()/deltaCoeffs();
    }

    // Set up temporary storage for the dd tensor (before inversion)
    tensorField dd(mesh().nCells(), tensor::zero);

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

        const labelList::subList FaceCells = pd.patchMesh().faceCells();

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

        const fvPatch& p = patchLsP.patchMesh();
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
        Info<< "surfaceInterpolation::makeLeastSquareVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}
*/

void surfaceInterpolation::makeLeastSquareVectors() const
{
    if (debug)
    {
        Info<< "surfaceInterpolation::makeLeastSquareVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    leastSquarePvectors_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquareP",
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *leastSquarePvectors_;

    leastSquareNvectors_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquareN",
            mesh().pointsInstance(),
            mesh()
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *leastSquareNvectors_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // Build the d-vectors
    surfaceVectorField d = mesh().Sf()/(mesh().magSf()*deltaCoeffs());

    if (!orthogonal())
    {
        d -= correctionVectors()/deltaCoeffs();
    }

    // Set up temporary storage for the dd tensor (before inversion)
    tensorField dd(mesh().nCells(), tensor::zero);

    forAll(owner, faceI)
    {
        tensor wdd = d[faceI]*d[faceI];

        dd[owner[faceI]] += wdd;

        // Yes, it is += because both vectors have the "wrong" sign
        dd[neighbour[faceI]] += wdd;
    }

    // Visit the boundaries. Coupled boundaries are taken into account
    // in the construction of d vectors.  
    forAll(d.boundaryField(), patchI)
    {
        const fvPatchVectorField& pd = d.boundaryField()[patchI];

        const labelList::subList FaceCells = pd.patchMesh().faceCells();

        forAll(pd, patchFaceI)
        {
            dd[FaceCells[patchFaceI]] += pd[patchFaceI]*pd[patchFaceI];
        }
    }

    // Invert the dd tensor
    tensorField invDd = inv(dd);


    // Revisit all faces and calculate the lsP and lsN vectors
    forAll(owner, faceI)
    {
        lsP[faceI] = invDd[owner[faceI]] & d[faceI];
        lsN[faceI] = -(invDd[neighbour[faceI]] & d[faceI]);
    }

    forAll(lsP.boundaryField(), patchI)
    {
        const fvPatchVectorField& pd = d.boundaryField()[patchI];

        fvPatchVectorField& patchLsP = lsP.boundaryField()[patchI];

        const fvPatch& p = patchLsP.patchMesh();
        const labelList::subList FaceCells = p.faceCells();

        forAll(p, patchFaceI)
        {
            patchLsP[patchFaceI] =
                invDd[FaceCells[patchFaceI]] & pd[patchFaceI];
        }
    }

    if (debug)
    {
        Info<< "surfaceInterpolation::makeLeastSquareVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
