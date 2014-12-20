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

Class
    volPointInterpolation

Description

\*---------------------------------------------------------------------------*/

#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "volFields.H"
#include "emptyFvPatch.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "parallelInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPointInterpolation, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void volPointInterpolation::makeWeights() const
{
    if (pointWeightsPtr_)
    {
        FatalErrorIn("volPointInterpolation::makeWeights() const")
            << "weighting factors already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    const pointField& points = fvMesh_.points();
    const labelListList& pointCells = fvMesh_.pointCells();
    const vectorField& cellCentres = fvMesh_.cellCentres();

    // Allocate storage for weighting factors
    pointWeightsPtr_ = new FieldField<Field, scalar>(points.size());
    FieldField<Field, scalar>& weightingFactors = *pointWeightsPtr_;

    forAll(weightingFactors, pointi)
    {
        weightingFactors.hook(new scalarField(pointCells[pointi].size()));
    }


    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        forAll(pointCells[pointi], pointCelli)
        {
            label celli = pointCells[pointi][pointCelli];

            weightingFactors[pointi][pointCelli] = 1.0/
                mag
                (
                    points[pointi] - cellCentres[celli]
                );
        }
    }

    pointScalarField volPointSumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            fvMesh_.polyMesh::instance(),
            fvMesh_
        ),
        pointMesh_,
        dimensionedScalar("zero", dimless, 0)
    );

    // Calculate weighting factors using inverse distance weighting
    forAll(points, pointi)
    {
        forAll(pointCells[pointi], pointCelli)
        {
            volPointSumWeights[pointi] += weightingFactors[pointi][pointCelli];
        }
    }

    // Update coupled boundaries
    // Work-around for cyclic parallels.  
    if (Pstream::parRun() && !vMesh().parallelData().cyclicParallel())
    {
        forAll (volPointSumWeights.boundaryField(), patchI)
        {
            if (volPointSumWeights.boundaryField()[patchI].coupled())
            {
                volPointSumWeights.boundaryField()[patchI].initAddField();
            }
        }

        forAll (volPointSumWeights.boundaryField(), patchI)
        {
            if (volPointSumWeights.boundaryField()[patchI].coupled())
            {
                volPointSumWeights.boundaryField()[patchI].addField
                (
                    volPointSumWeights.internalField()
                );
            }
        }
    }

    forAll(points, pointi)
    {
        forAll(pointCells[pointi], pointCelli)
        {
            weightingFactors[pointi][pointCelli] /= volPointSumWeights[pointi];
        }
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


void volPointInterpolation::makeBoundaryAddressing() const
{
    if (boundaryPointsPtr_)
    {
        FatalErrorIn
        (
            "void volPointInterpolation::"
            "makeBoundaryAddressing() const"
        )   << "boundary addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeBoundaryAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    // Go through all patches and mark up the external edge points
    labelHashSet pointsCorrectionMap(max(vMesh().nPoints()/10, 100));

    const fvBoundaryMesh& bm = vMesh().boundary();

    forAll (bm, patchI)
    {
        // If the patch is empty, skip it
        // If the patch is coupled, and there are no cyclic parallels, skip it
        if
        (
            bm[patchI].type() != emptyFvPatch::typeName
//         && !bm[patchI].coupled()
        && !(
                bm[patchI].coupled()
             && Pstream::parRun()
             && !vMesh().parallelData().cyclicParallel()
            )
        )
        {
            const labelList& bp = bm[patchI].patch().boundaryPoints();

            const labelList& meshPoints = bm[patchI].patch().meshPoints();

            forAll (bp, pointI)
            {
                pointsCorrectionMap.insert(meshPoints[bp[pointI]]);
            }
        }
    }

    // Grab the points to correct
    boundaryPointsPtr_ = new labelList(pointsCorrectionMap.toc());
    labelList& ptc = *boundaryPointsPtr_;

    sort(ptc);

    if (debug)
    {
        Info<< "volPointInterpolation::makeBoundaryAddressing() : "
            << "finished constructing boundary addressing"
            << endl;
    }
}


void volPointInterpolation::makeBoundaryWeights() const
{
    if (pointBoundaryWeightsPtr_)
    {
        FatalErrorIn
        (
            "void volPointInterpolation::"
            "makeBoundaryWeights() const"
        )   << "boundary weighting factors already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeBoundaryWeights() : "
            << "constructing boundary weighting factors"
            << endl;
    }

    const labelList& ptc = boundaryPoints();

    // Calculate the correction vectors and weighting factors

    pointBoundaryWeightsPtr_ = new FieldField<Field, scalar>(ptc.size());
    FieldField<Field, scalar>& w = *pointBoundaryWeightsPtr_;

    const labelListList& pf = vMesh().pointFaces();

    const volVectorField& centres = vMesh().C();

    const fvBoundaryMesh& bm = vMesh().boundary();

    pointScalarField volPointSumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            vMesh().polyMesh::instance(),
            vMesh()
        ),
        pMesh(),
        dimensionedScalar("zero", dimless, 0)
    );

    forAll (ptc, pointI)
    {
        const label curPoint = ptc[pointI];

        const labelList& curFaces = pf[curPoint];

        w.hook(new scalarField(curFaces.size()));
        scalarField& curWeights = w[pointI];

        label nFacesAroundPoint = 0;

        const vector& pointLoc = vMesh().points()[curPoint];

        // Go through all the faces
        forAll (curFaces, faceI)
        {
            if (!vMesh().isInternalFace(curFaces[faceI]))
            {
                // This is a boundary face.  If not in the empty patch
                // or coupled calculate the extrapolation vector
                label patchID =
                    vMesh().boundaryMesh().whichPatch(curFaces[faceI]);

                if
                (
                    bm[patchID].type() != emptyFvPatch::typeName
                && !(
                        bm[patchID].coupled()
                    && Pstream::parRun()
                    && !vMesh().parallelData().cyclicParallel()
                    )
//                 && !bm[patchID].coupled()
                )
                {
                    curWeights[nFacesAroundPoint] =
                        1.0/mag
                        (
                            pointLoc
                          - centres.boundaryField()[patchID]
                            [
                                bm[patchID].patch().whichFace(curFaces[faceI])
                            ]
                        );

                    nFacesAroundPoint++;
                }
            }
        }

        // Reset the sizes of the local weights
        curWeights.setSize(nFacesAroundPoint);

        // Collect the sum of weights for parallel correction
        volPointSumWeights[curPoint] += sum(curWeights);
    }

    // Do parallel correction of weights

    // Update coupled boundaries
    // Work-around for cyclic parallels.  
    if (Pstream::parRun() && !vMesh().parallelData().cyclicParallel())
    {
        forAll (volPointSumWeights.boundaryField(), patchI)
        {
            if (volPointSumWeights.boundaryField()[patchI].coupled())
            {
                volPointSumWeights.boundaryField()[patchI].initAddField();
            }
        }

        forAll (volPointSumWeights.boundaryField(), patchI)
        {
            if (volPointSumWeights.boundaryField()[patchI].coupled())
            {
                volPointSumWeights.boundaryField()[patchI].addField
                (
                    volPointSumWeights.internalField()
                );
            }
        }
    }

    // Re-scale the weights for the current point
    forAll (ptc, pointI)
    {
        w[pointI] /= volPointSumWeights[ptc[pointI]];
    }

    if (debug)
    {
        Info<< "volPointInterpolation::makeBoundaryWeights() : "
            << "finished constructing boundary weighting factors"
            << endl;
    }
}


// Do what is neccessary if the mesh has changed topology
void volPointInterpolation::clearAddressing() const
{
    deleteDemandDrivenData(patchInterpolatorsPtr_);
    deleteDemandDrivenData(boundaryPointsPtr_);
}


// Do what is neccessary if the mesh has moved
void volPointInterpolation::clearGeom() const
{
    deleteDemandDrivenData(pointWeightsPtr_);
    deleteDemandDrivenData(pointBoundaryWeightsPtr_);
}


const PtrList<primitivePatchInterpolation>&
volPointInterpolation::patchInterpolators() const
{
    if (!patchInterpolatorsPtr_)
    {
        const fvBoundaryMesh& bdry = fvMesh_.boundary();

        patchInterpolatorsPtr_ =
            new PtrList<primitivePatchInterpolation>(bdry.size());

        forAll (bdry, patchI)
        {
            patchInterpolatorsPtr_->hook
            (
                new primitivePatchInterpolation(bdry[patchI].patch())
            );
        }
    }

    return *patchInterpolatorsPtr_;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

volPointInterpolation::volPointInterpolation
(
    const fvMesh& vm,
    const pointMesh& pm
)
:
    fvMesh_(vm),
    pointMesh_(pm),
    pointWeightsPtr_(NULL),
    patchInterpolatorsPtr_(NULL),
    boundaryPointsPtr_(NULL),
    pointBoundaryWeightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

volPointInterpolation::~volPointInterpolation()
{
    clearAddressing();
    clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return point weights
const FieldField<Field, scalar>& volPointInterpolation::pointWeights() const
{
    // If weighting factor array not present then construct
    if (!pointWeightsPtr_)
    {
        makeWeights();
    }

    return *pointWeightsPtr_;
}


// Return boundary points
const labelList& volPointInterpolation::boundaryPoints() const
{
    // If weighting factor array not present then construct
    if (!boundaryPointsPtr_)
    {
        makeBoundaryAddressing();
    }

    return *boundaryPointsPtr_;
}


// Return boundary point weights
const FieldField<Field, scalar>&
volPointInterpolation::pointBoundaryWeights() const
{
    // If weighting factor array not present then construct
    if (!pointBoundaryWeightsPtr_)
    {
        makeBoundaryWeights();
    }

    return *pointBoundaryWeightsPtr_;
}


// Do what is neccessary if the mesh has moved
void volPointInterpolation::updateTopology()
{
    clearAddressing();
    clearGeom();
}


// Do what is neccessary if the mesh has moved
bool volPointInterpolation::movePoints()
{
    clearGeom();

    if (patchInterpolatorsPtr_)
    {
        PtrList<primitivePatchInterpolation>& pi = *patchInterpolatorsPtr_;

        forAll (pi, patchI)
        {
            pi[patchI].movePoints();
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
