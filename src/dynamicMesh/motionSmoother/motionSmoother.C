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

\*---------------------------------------------------------------------------*/

#include "motionSmoother.H"
#include "twoDPointCorrector.H"
#include "faceSet.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(motionSmoother, 0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::motionSmoother::checkFld(const pointScalarField& fld)
{
    forAll(fld, pointI)
    {
        const scalar val = fld[pointI];

        if ((val > -GREAT) && (val < GREAT))
        {}
        else
        {
            FatalErrorIn("motionSmoother::checkFld")
                << "Problem : point:" << pointI << " value:" << val
                << abort(FatalError);
        }
    }
}


Foam::labelHashSet Foam::motionSmoother::getPoints
(
    const labelHashSet& faceLabels
) const
{
    labelHashSet usedPoints(mesh_.nPoints()/100);

    for
    (
        labelHashSet::const_iterator iter = faceLabels.begin();
        iter != faceLabels.end();
        ++iter
    )
    {
        const face& f = mesh_.faces()[iter.key()];

        forAll(f, fp)
        {
            usedPoints.insert(f[fp]);
        }
    }

    return usedPoints;
}

// Smooth on selected points (usually patch points)
void Foam::motionSmoother::minSmooth
(
    const labelList& meshPoints,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    const labelListList& pointEdges = mesh_.pointEdges();
    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();

    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        // get min of current value and avg of neighbour values
        const labelList& pEdges = pointEdges[pointI];

        // unweighted average
        scalar nbrAvg = avg(fld, edges, points, pEdges, pointI);

        newFld[pointI] = min(fld[pointI], 0.5*fld[pointI] + 0.5*nbrAvg);
    }
}


// Smooth on all internal points
void Foam::motionSmoother::minSmooth
(
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    const labelListList& pointEdges = mesh_.pointEdges();
    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();

    forAll(fld, pointI)
    {
        if (isInternalPoint(pointI))
        {
            // get min of current value and avg of neighbour values
            const labelList& pEdges = pointEdges[pointI];

            // unweighted average
            scalar nbrAvg = avg(fld, edges, points, pEdges, pointI);

            newFld[pointI] = min(fld[pointI], 0.5*fld[pointI] + 0.5*nbrAvg);
        }
    }

    newFld.correctBoundaryConditions();
}


// Scale on selected points
void Foam::motionSmoother::scaleField
(
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    for
    (
        labelHashSet::const_iterator iter = pointLabels.begin();
        iter != pointLabels.end();
        ++iter
    )
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] *= scale;
        }
    }
    fld.correctBoundaryConditions();
}


// Scale on selected points (usually patch points)
void Foam::motionSmoother::scaleField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        if (pointLabels.found(pointI))
        {
            fld[pointI] *= scale;
        }
    }
}


bool Foam::motionSmoother::isInternalPoint(const label pointI) const
{
    return isInternalPoint_.get(pointI) == 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::motionSmoother::motionSmoother
(
    polyMesh& mesh,
    pointMesh& pMesh,
    indirectPrimitivePatch& pp,
    const labelList& adaptPatchIDs,
    const scalar reduction,
    const label nSmoothScale,
    const scalar maxNonOrtho,
    const scalar minVol,
    const scalar maxConcave,
    const scalar minArea
)
:
    mesh_(mesh),
    pMesh_(pMesh),
    pp_(pp),
    adaptPatchIDs_(adaptPatchIDs),
    reduction_(reduction),
    nSmoothScale_(nSmoothScale),
    maxNonOrtho_(maxNonOrtho),
    minVol_(minVol),
    maxConcave_(maxConcave),
    minArea_(minArea),
    displacement_
    (
        IOobject
        (
            "displacement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
    ),
    scale_
    (
        IOobject
        (
            "scale",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(mesh_.points()),
    isInternalPoint_(mesh_.nPoints(), 1),
    twoDMotion_(false),
    correct2DPtr_(NULL)
{
    const pointBoundaryMesh& patches = pMesh_.boundary();

    // Get twoDcorrector stuff
    IOdictionary motionProps
    (
        IOobject
        (
            "motionProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    twoDMotion_ = motionProps.lookup("twoDMotion");

    if (twoDMotion_)
    {
        correct2DPtr_ = new twoDPointCorrector(mesh_);
    }

    // Determine internal points. Note that for twoD there are no internal
    // points so we use the points of adaptPatchIDs instead

    if (twoDMotion_)
    {
        const labelList& meshPoints = pp.meshPoints();

        forAll(meshPoints, i)
        {
            isInternalPoint_.set(meshPoints[i], 0);
        }
    }
    else
    {
        forAll(patches, patchI)
        {
            const pointPatch& pp = patches[patchI];

            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, i)
            {
                isInternalPoint_.set(meshPoints[i], 0);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSmoother::~motionSmoother()
{
    deleteDemandDrivenData(correct2DPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::motionSmoother::mesh() const
{
    return mesh_;
}


const Foam::indirectPrimitivePatch& Foam::motionSmoother::patch() const
{
    return pp_;
}


Foam::scalar Foam::motionSmoother::maxNonOrtho() const
{
    return maxNonOrtho_;
}


Foam::scalar Foam::motionSmoother::minVol() const
{
    return minVol_;
}


Foam::scalar Foam::motionSmoother::maxConcave() const
{
    return maxConcave_;
}


Foam::scalar Foam::motionSmoother::minArea() const
{
    return minArea_;
}


Foam::pointVectorField& Foam::motionSmoother::displacement()
{
    return displacement_;
}


const Foam::pointScalarField& Foam::motionSmoother::scale() const
{
    return scale_;
}


const Foam::pointField& Foam::motionSmoother::oldPoints() const
{
    return oldPoints_;
}


void Foam::motionSmoother::correct()
{
    oldPoints_ = mesh_.points();

    scale_ = 1.0;

    // No need to update twoDmotion corrector since only holds edge labels
    // which will remain the same as before. So unless the mesh was distorted
    // severely outside of motionSmoother there will be no need.
}


Foam::tmp<Foam::scalarField> Foam::motionSmoother::movePoints
(
    pointField& newPoints
)
{
    // Correct for 2-D motion
    if (correct2DPtr_)
    {
        Info << "Correct-ing 2-D mesh motion";

        // We do not want to move 3D planes so project all points onto those
        const pointField& oldPoints = mesh_.points();
        const edgeList& edges = mesh_.edges();

        const labelList& neIndices = twoDCorrector().normalEdgeIndices();
        const vector& pn = twoDCorrector().planeNormal();

        forAll(neIndices, i)
        {
            const edge& e = edges[neIndices[i]];

            point& pStart = newPoints[e.start()];

            pStart += pn*(pn & (oldPoints[e.start()] - pStart));

            point& pEnd = newPoints[e.end()];

            pEnd += pn*(pn & (oldPoints[e.end()] - pEnd));
        }

        // Correct tangentially
        correct2DPtr_->correctPoints(newPoints);
        Info << " ...done" << endl;
    }    

    tmp<scalarField> tsweptVol = mesh_.movePoints(newPoints);

    pp_.movePoints(newPoints);

    pMesh_.movePoints();

    return tsweptVol;
}


bool Foam::motionSmoother::scaleMesh
(
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    if (!smoothMesh && adaptPatchIDs_.size() == 0)
    {
        FatalErrorIn("motionSmoother::scaleMesh(const bool")
            << "You specified both no movement on the internal mesh points"
            << " (smoothMesh = false)" << nl
            << "and no movement on the patch (adaptPatchIDs is empty)" << nl
            << "Hence nothing to adapt."
            << exit(FatalError);
    }

    cpuTime timer;

    // Set newPoints as old + scale*displacement
    pointField newPoints
    (
        oldPoints_
      + scale_.internalField()*displacement_.internalField()
    );
    Info<< "Calculated new mesh position in = "
        << timer.cpuTimeIncrement() << " s\n" << nl << endl;

    {
        Info<< "scale:" << "  min:" << min(scale_.internalField())
            << "  max:" << max(scale_.internalField()) << endl;

        pointScalarField magDisp(mag(displacement_));

        Info<< "disp:" << "  min:" << min(magDisp.internalField())
            << "  max:" << max(magDisp.internalField()) << endl;
    }



    // Move
    movePoints(newPoints);
    Info<< "Moved mesh in = "
        << timer.cpuTimeIncrement() << " s\n" << nl << endl;

    // Check
    faceSet wrongFaces(mesh_, "wrongFaces", mesh_.nFaces()/100+100);
    checkMesh(wrongFaces);

    Info<< "Checked mesh in = "
        << timer.cpuTimeIncrement() << " s\n" << nl << endl;

    if (wrongFaces.size() <= nAllowableErrors)
    {
        return true;
    }
    else
    {
        Info<< "Writing faceSet with incorrect faces to " << wrongFaces.name()
            << endl;
        wrongFaces.write();

        // Find out points used by wrong faces and scale displacement.
        labelHashSet usedPoints(getPoints(wrongFaces));

        if (adaptPatchIDs_.size() != 0)
        {
            // Scale conflicting patch points
            scaleField(pp_.meshPoints(), usedPoints, reduction_, scale_);
        }
        if (smoothMesh)
        {
            // Scale conflicting internal points
            scaleField(usedPoints, reduction_, scale_);
        }

        for (label i = 0; i < nSmoothScale_; i++)
        {
            if (adaptPatchIDs_.size() != 0)
            {
                // Smooth patch values
                pointScalarField oldScale(scale_);
                minSmooth
                (
                    pp_.meshPoints(),
                    oldScale,
                    scale_
                );
                checkFld(scale_);
            }
            if (smoothMesh)
            {
                // Smooth internal values
                pointScalarField oldScale(scale_);
                minSmooth(oldScale, scale_);
                checkFld(scale_);
            }
        }
        Info<< "After smoothing :"
            << " min:" << Foam::gMin(scale_)
            << " max:" << Foam::gMax(scale_)
            << endl;

        Info<< nl << endl;

        return false;
    }
}


void Foam::motionSmoother::updateTopology()
{
    if (correct2DPtr_)
    {
        correct2DPtr_->updateTopology();
    }    
}


bool Foam::motionSmoother::checkMesh(labelHashSet& wrongFaces) const
{
    if (maxNonOrtho_ < 180.0-SMALL)
    {
        Info<< "Checking non orthogonality" << endl;

        label nOldSize = wrongFaces.size();
        mesh_.setOrthWarn(maxNonOrtho_);
        mesh_.checkFaceDotProduct(false, &wrongFaces);

        Info<< "Detected " << wrongFaces.size() - nOldSize
            << " faces with non-orthogonality > " << maxNonOrtho_ << " degrees"
            << endl;
    }

    if (minVol_ > -GREAT)
    {
        Info<< "Checking face pyramids" << endl;

        label nOldSize = wrongFaces.size();
        mesh_.checkFacePyramids(false, minVol_, &wrongFaces);
        Info<< "Detected additional " << wrongFaces.size() - nOldSize
            << " faces with illegal face pyramids" << endl;
    }

    if (maxConcave_ < 180.0-SMALL)
    {
        Info<< "Checking face angles" << endl;

        label nOldSize = wrongFaces.size();
        mesh_.checkFaceAngles(false, maxConcave_, &wrongFaces);
        Info<< "Detected additional " << wrongFaces.size() - nOldSize
            << " faces with concavity > " << maxConcave_ << " degrees"
            << endl;
    }

    if (minArea_ > -SMALL)
    {
        Info<< "Checking face areas" << endl;

        label nOldSize = wrongFaces.size();

        const scalarField magFaceAreas = mag(mesh_.faceAreas());

        forAll(magFaceAreas, faceI)
        {
            if (magFaceAreas[faceI] < minArea_)
            {
                wrongFaces.insert(faceI);
            }
        }
        Info<< "Detected additional " << wrongFaces.size() - nOldSize
            << " faces with area < " << minArea_ << " m^2" << endl;
    }

    return wrongFaces.size() > 0;
}


// ************************************************************************* //
