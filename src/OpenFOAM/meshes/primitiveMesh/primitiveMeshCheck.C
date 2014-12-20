/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "primitiveMesh.H"
#include "pyramidPointFaceRef.H"
#include "cell.H"
#include "mathematicalConstants.H"
#include "boolList.H"
#include "labelHashSet.H"
#include "ListOps.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar primitiveMesh::orthWarn_ = 70; // deg
scalar primitiveMesh::skewWarn_ = 2;
scalar primitiveMesh::aspectWarn_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool primitiveMesh::checkClosedBoundary(const bool report) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkClosedBoundary("
            << "const bool) const: "
            << "checking whether the boundary is closed" << endl;
    }

    // Loop through all boundary faces and sum up the face area vectors.
    // For a closed boundary, this should be zero in all vector components

    vector sumClosed(vector::zero);
    scalar sumMagClosedBoundary = 0;

    const vectorField& areas = faceAreas();

    for (label faceI = nInternalFaces(); faceI < areas.size(); faceI++)
    {
        sumClosed += areas[faceI];
        sumMagClosedBoundary += mag(areas[faceI]);
    }

    // Check the openness in the maximum direction (this is APPROXIMATE!)
    scalar maxOpen = max
    (
        sumClosed.component(vector::X),
        max
        (
            sumClosed.component(vector::Y),
            sumClosed.component(vector::Z)
        )
    );

    reduce(sumClosed, sumOp<vector>());
    reduce(maxOpen, maxOp<scalar>());

    if (maxOpen > closedTolerance_*max(1.0, sumMagClosedBoundary))
    {
        SeriousErrorIn
        (
            "primitiveMesh::checkClosedBoundary(const bool report) const"
        )   << "Possible hole in boundary description" << endl;

        Info<< "Boundary openness in x-direction = "
            << sumClosed.component(vector::X) << endl;

        Info<< "Boundary openness in y-direction = "
            << sumClosed.component(vector::Y) << endl;

        Info<< "Boundary openness in z-direction = "
            << sumClosed.component(vector::Z) << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Boundary openness in x-direction = "
                << sumClosed.component(vector::X) << endl;

            Info<< "Boundary openness in y-direction = "
                << sumClosed.component(vector::Y) << endl;

            Info<< "Boundary openness in z-direction = "
                << sumClosed.component(vector::Z) << endl;

            Info<< "Boundary closed (OK)." << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkClosedCells
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkClosedCells("
            << "const bool, labelHashSet*) const: "
            << "checking whether cells are closed" << endl;
    }

    // Check that all cells labels are valid
    const cellList& c = cells();

    label nErrorClosed = 0;

    forAll (c, cI)
    {
        const cell& curCell = c[cI];

        if (min(curCell) < 0 || max(curCell) > nFaces())
        {
            WarningIn
            (
                "bool primitiveMesh::checkClosedCells("
                "const bool, labelHashSet*) const"
            )   << "Cell " << cI << " contains face labels out of range: "
                << curCell << " Max face index = " << nFaces() << endl;

            if (setPtr)
            {
                setPtr->insert(cI);
            }

            nErrorClosed++;
        }
    }

    if (nErrorClosed > 0)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkClosedCells("
            "const bool report, labelHashSet*) const"
        )  << nErrorClosed << " cells with invalid face labels found"
            << endl;

        return true;
    }

    // Loop through cell faces and sum up the face area vectors for each cell.
    // This should be zero in all vector components

    vectorField sumClosed(nCells(), vector::zero);

    scalarField sumMagClosed(nCells(), 0);

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    const vectorField& areas = faceAreas();

    forAll (own, faceI)
    {
        // Add to owner
        sumClosed[own[faceI]] += areas[faceI];
        sumMagClosed[own[faceI]] += mag(areas[faceI]);
    }

    forAll (nei, faceI)
    {
        // Subtract from neighbour
        sumClosed[nei[faceI]] -= areas[faceI];
        sumMagClosed[nei[faceI]] += mag(areas[faceI]);
    }

    label nOpen = 0;
    scalar maxOpenCell = 0;

    label nAspect = 0;
    scalar maxAspectRatio = 0;

    const scalarField& vols = cellVolumes();

    // Check the sums
    forAll (sumClosed, cellI)
    {
        scalar maxOpen = max
        (
            sumClosed[cellI].component(vector::X),
            max
            (
                sumClosed[cellI].component(vector::Y),
                sumClosed[cellI].component(vector::Z)
            )
        );

        maxOpenCell = max(maxOpenCell, maxOpen);

        if (maxOpen > closedTolerance_*max(1.0, sumMagClosed[cellI]))
        {
            if (debug || report)
            {
                Pout<< "Cell " << cellI << " is not closed. "
                    << "Face area vectors sum up to " << mag(sumClosed[cellI])
                    << " directionwise " << sumClosed[cellI] << " or "
                    << mag(sumClosed[cellI])
                       /(mag(sumMagClosed[cellI]) + VSMALL)
                    << " of the area of all faces of the cell. " << endl
                    << "    Area magnitudes sum up to "
                    << sumMagClosed[cellI] << endl;
            }

            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nOpen++;
        }

        scalar aspectRatio =
            1.0/6.0*sumMagClosed[cellI]/pow(vols[cellI], 2.0/3.0);

        maxAspectRatio = max(maxAspectRatio, aspectRatio);

        if (aspectRatio > aspectWarn_)
        {
            if (debug || report)
            {
                Pout<< "High aspect ratio for cell " << cellI
                    << ": " << aspectRatio << endl;
            }

            nAspect++;
        }
    }

    reduce(nOpen, sumOp<label>());
    reduce(maxOpenCell, maxOp<scalar>());

    reduce(nAspect, sumOp<label>());
    reduce(maxAspectRatio, maxOp<scalar>());


    if (nOpen > 0)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkClosedCells("
            "const bool report, labelHashSet*) const"
        )   << nOpen << " open cells found. Max cell openness: "
            << maxOpenCell << endl;

        return true;
    }

    if (nAspect > 0)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkClosedCells("
            "const bool report, labelHashSet*) const"
        )   << nAspect << " high aspect ratio cells found.  "
            << "Max aspect ratio: " << maxAspectRatio
            << endl;

        return true;
    }

    if (debug || report)
    {
        Info<< "Max cell openness = " << maxOpenCell
            << "  Max aspect ratio = " << maxAspectRatio
            << ".  All cells OK.\n" << endl;
    }

    return false;
}


bool primitiveMesh::checkFaceAreas
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceAreas("
            << "const bool, labelHashSet*) const: "
            << "checking face area magnitudes" << endl;
    }

    const scalarField magFaceAreas = mag(faceAreas());

    const labelList& own = faceOwner();

    const labelList& nei = faceNeighbour();

    scalar minArea = GREAT;
    scalar maxArea = -GREAT;

    forAll (magFaceAreas, faceI)
    {
        if (magFaceAreas[faceI] < VSMALL)
        {
            if (debug || report)
            {
                if (isInternalFace(faceI))
                {
                    Pout<< "Zero or negative face area detected for "
                        << "internal face " << faceI << " between cells "
                        << own[faceI] << " and " << nei[faceI]
                        << ".  Face area magnitude = "
                        << magFaceAreas[faceI] << endl;
                }
                else
                {
                    Pout<< "Zero or negative face area detected for "
                        << "boundary face " << faceI << " next to cell "
                        << own[faceI] << ".  Face area magnitude = "
                        << magFaceAreas[faceI] << endl;
                }
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
            }
        }

        minArea = min(minArea, magFaceAreas[faceI]);
        maxArea = max(maxArea, magFaceAreas[faceI]);
    }

    reduce(minArea, minOp<scalar>());
    reduce(maxArea, maxOp<scalar>());

    if (minArea < VSMALL)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkFaceAreas("
            "const bool report, labelHashSet*) const"
        )   << "Zero or negative face area detected.  Minimum negative area: " 
            << minArea << ". This mesh is invalid"
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Minumum face area = " << minArea
                << ". Maximum face area = " << maxArea
                << ".  Face area magnitudes OK.\n" << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkCellVolumes
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellVolumes("
            << "const bool, labelHashSet*) const: "
            << "checking cell volumes" << endl;
    }

    const scalarField& vols = cellVolumes();

    scalar minVolume = GREAT;
    scalar maxVolume = -GREAT;

    label nNegVolCells = 0;

    forAll (vols, cellI)
    {
        if (vols[cellI] < VSMALL)
        {
            if (debug || report)
            {
                Pout<< "Zero or negative cell volume detected for cell "
                    << cellI << ".  Volume = " << vols[cellI] << endl;
            }

            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nNegVolCells++;
        }

        minVolume = min(minVolume, vols[cellI]);
        maxVolume = max(maxVolume, vols[cellI]);
    }

    reduce(minVolume, minOp<scalar>());
    reduce(maxVolume, maxOp<scalar>());
    reduce(nNegVolCells, sumOp<label>());

    if (minVolume < VSMALL)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkCellVolumes("
            "const bool report, labelHashSet*) const"
        )   << "Zero or negative cell volume detected.  "
            << "Minimum negative volume: " 
            << minVolume << ".\nNumber of negative volume cells: "
            << nNegVolCells << ".  This mesh is invalid"
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Min volume = " << minVolume
                << ". Max volume = " << maxVolume
                << ".  Total volume = " << sum(vols)
                << ".  Cell volumes OK.\n" << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkFaceDotProduct
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceDotProduct("
            << "const bool, labelHashSet*) const: "
            << "checking mesh non-orthogonality" << endl;
    }

    // for all internal faces check theat the d dot S product is positive
    const vectorField& centres = cellCentres();
    const vectorField& areas = faceAreas();

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(orthWarn_/180.0*mathematicalConstant::pi);

    scalar minDDotS = GREAT;

    scalar sumDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    forAll (nei, faceI)
    {
        vector d = centres[nei[faceI]] - centres[own[faceI]];
        const vector& s = areas[faceI];

        scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

        if (dDotS < severeNonorthogonalityThreshold)
        {
            if (dDotS > SMALL)
            {
                if (debug || report)
                {
                    // Severe non-orthogonality but mesh still OK
                    Pout<< "Severe non-orthogonality for face " << faceI
                        << " between cells " << own[faceI]
                        << " and " << nei[faceI]
                        << ": Angle = "
                        << ::acos(dDotS)/mathematicalConstant::pi*180.0
                        << " deg." << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                severeNonOrth++;
            }
            else
            {
                // Non-orthogonality greater than 90 deg
                WarningIn
                (
                    "primitiveMesh::checkFaceDotProduct"
                    "(const bool report, labelHashSet* setPtr) const"
                )   << "Severe non-orthogonality detected for face " << faceI
                    << " between cells " << own[faceI] << " and " << nei[faceI]
                    << ": Angle = " << ::acos(dDotS)/mathematicalConstant::pi*180.0
                    << " deg." << endl;

                errorNonOrth++;


                if (setPtr)
                {
                    setPtr->insert(faceI);
                }
            }
        }

        if (dDotS < minDDotS)
        {
            minDDotS = dDotS;
        }

        sumDDotS += dDotS;
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    label neiSize = nei.size();
    reduce(neiSize, sumOp<label>());

    // Only report if there are some internal faces
    if (neiSize > 0)
    {
        if (minDDotS < severeNonorthogonalityThreshold)
        {
            Info<< "Number of non-orthogonality errors: " << errorNonOrth
                << ". Number of severely non-orthogonal faces: "
                << severeNonOrth  << "." << endl;
        }
    }

    if (debug || report)
    {
        if (neiSize > 0)
        {
            Info<< "Mesh non-orthogonality Max: "
                << ::acos(minDDotS)/mathematicalConstant::pi*180.0
                << " average: " <<
                   ::acos(sumDDotS/neiSize)/mathematicalConstant::pi*180.0
                << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        SeriousErrorIn
        (
            "primitiveMesh::checkFaceDotProduct"
            "(const bool report, labelHashSet* setPtr) const"
        )   << "Error in non-orthogonality detected" << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Non-orthogonality check OK.\n" << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFacePyramids("
            << "const bool, const scalar, labelHashSet*) const: "
            << "checking face orientation" << endl;
    }

    // check whether face area vector points to the cell with higher label
    const vectorField& ctrs = cellCentres();

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    const faceList& f = faces();

    const pointField& p = points();

    label nErrorPyrs = 0;

    forAll (f, faceI)
    {
        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = pyramidPointFaceRef(f[faceI], ctrs[own[faceI]]).mag(p);

        if (pyrVol > -minPyrVol)
        {
            if (debug || report)
            {
                Pout<< "bool primitiveMesh::checkFacePyramids("
                    << "const bool, const scalar, labelHashSet*) const: "
                    << "face " << faceI << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVol
                    << " Face " << f[faceI] << " area: " << f[faceI].mag(p)
                    << " Owner cell: " << own[faceI] << endl
                    << "Owner cell vertex labels: "
                    << cells()[own[faceI]].labels(f)
                    << endl;
            }


            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nErrorPyrs++;
        }

        if (isInternalFace(faceI))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef(f[faceI], ctrs[nei[faceI]]).mag(p);

            if (pyrVol < minPyrVol)
            {
                if (debug || report)
                {
                    Pout<< "bool primitiveMesh::checkFacePyramids("
                        << "const bool, const scalar, labelHashSet*) const: "
                        << "face " << faceI << " points the wrong way. " << endl
                        << "Pyramid volume: " << -pyrVol
                        << " Face " << f[faceI] << " area: " << f[faceI].mag(p)
                        << " Neighbour cell: " << nei[faceI] << endl
                        << "Neighbour cell vertex labels: "
                        << cells()[nei[faceI]].labels(f)
                        << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                nErrorPyrs++;
            }
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkFacePyramids("
            "const bool, const scalar, labelHashSet*) const"
        )   << "Error in face pyramids: faces pointing the wrong way!"
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Face pyramids OK.\n" << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkFaceSkewness
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceSkewnesss("
            << "const bool, labelHashSet*) const: "
            << "checking face skewness" << endl;
    }

    // Warn if the skew correction vector is more than skewWarning times
    // larger than the face area vector

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    const vectorField& cellCtrs = cellCentres();
    const vectorField& faceCtrs = faceCentres();

    scalar maxSkew = 0;

    label nWarnSkew = 0;

    forAll(nei, faceI)
    {
        scalar dOwn = mag(faceCtrs[faceI] - cellCtrs[own[faceI]]);
        scalar dNei = mag(faceCtrs[faceI] - cellCtrs[nei[faceI]]);

        point faceIntersection =
            cellCtrs[own[faceI]]*dNei/(dOwn+dNei)
          + cellCtrs[nei[faceI]]*dOwn/(dOwn+dNei);

        scalar skewness = 
            mag(faceCtrs[faceI] - faceIntersection)
            /(mag(cellCtrs[nei[faceI]] - cellCtrs[own[faceI]]) + VSMALL);

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skewness > skewWarn_)
        {
            if (debug || report)
            {
                Pout<< " Severe skewness for face " << faceI
                    << " skewness = " << skewness << endl;
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nWarnSkew++;
        }

        if(skewness > maxSkew)
        {
            maxSkew = skewness;
        }
    }


    // Boundary faces: consider them to have only skewness error.
    // (i.e. treat as if mirror cell on other side)

    const vectorField& fAreas = faceAreas();

    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        vector faceNormal = fAreas[faceI];
        faceNormal /= mag(faceNormal) + VSMALL;

        vector dOwn = faceCtrs[faceI] - cellCtrs[own[faceI]];

        vector dWall = faceNormal*(faceNormal & dOwn);

        point faceIntersection = cellCtrs[own[faceI]] + dWall;

        scalar skewness =
            mag(faceCtrs[faceI] - faceIntersection)
            /(2*mag(dWall) + VSMALL);

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skewness > skewWarn_)
        {
            if (debug || report)
            {
                Pout<< " Severe skewness for boundary face " << faceI
                    << " skewness = " << skewness << endl;
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nWarnSkew++;
        }

        if(skewness > maxSkew)
        {
            maxSkew = skewness;
        }
    }


    reduce(maxSkew, maxOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if (nWarnSkew > 0)
    {
        WarningIn
        (
            "primitiveMesh::checkFaceSkewness"
            "(const bool report, labelHashSet* setPtr) const"
        )   << "Large face skewness detected.  Max skewness = " << 100*maxSkew
            << " percent.\nThis may impair the quality of the result." << nl
            << nWarnSkew << " highly skew faces detected."
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Max skewness = " << 100*maxSkew
                << " percent.  Face skewness OK.\n" << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkPoints
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkPoints"
            << "(const bool, labelHashSet*) const: "
            << "checking points" << endl;
    }

    label nFaceErrors = 0;
    label nCellErrors = 0;

    const labelListList& pf = pointFaces();

    forAll (pf, pointI)
    {
        if (pf[pointI].size() == 0)
        {
            WarningIn
            (
                "bool primitiveMesh::checkPoints"
                "(const bool, labelHashSet*) const"
            )   << "Point " << pointI << " not used by any faces." << endl;

            if (setPtr)
            {
                setPtr->insert(pointI);
            }

            nFaceErrors++;
        }
    }

    const labelListList& pc = pointCells();

    forAll (pc, pointI)
    {
        if (pc[pointI].size() == 0)
        {
            WarningIn
            (
                "bool primitiveMesh::checkPoints"
                "(const bool, labelHashSet*) const"
            )   << "Point " << pointI << " not used by any cells." << endl;

            if (setPtr)
            {
                setPtr->insert(pointI);
            }
            nCellErrors++;
        }
    }

    reduce(nFaceErrors, sumOp<label>());
    reduce(nCellErrors, sumOp<label>());

    if (nFaceErrors > 0 || nCellErrors > 0)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkPoints"
            "(const bool, labelHashSet*) const"
        )   << "Error in point usage detected: " << nFaceErrors
            << " unused points found in the mesh.  This mesh is invalid."
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info << "Point usage check OK.\n" << endl;
        }

        return false;
    }
}


// Check convexity of angles in a face. Allow a slight non-convexity.
// E.g. maxDeg = 10 allows for angles < 190 (or 10 degrees concavity)
// (if truly concave and points not visible from face centre the face-pyramid
//  check in checkMesh will fail)
bool primitiveMesh::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceAngles"
            << "(const bool, const scalar, labelHashSet*) const: "
            << "checking face angles" << endl;
    }

    if (maxDeg < -SMALL || maxDeg > 180+SMALL)
    {
        FatalErrorIn
        (
            "primitiveMesh::checkFaceAngles"
            "(const bool, const scalar, labelHashSet*)"
        )   << "maxDeg should be [0..180] but is now " << maxDeg
            << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(maxDeg/180.0*mathematicalConstant::pi);

    const pointField& p = points();
    const faceList& fcs = faces();
    vectorField faceNormals(faceAreas());
    faceNormals /= mag(faceNormals) + VSMALL;

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    label errorFaceI = -1;

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f[0]] - p[f[f.size()-1]]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + VSMALL;

        forAll(f, fp0)
        {
            // Get vertex after fp
            label fp1 = (fp0 + 1) % f.size();

            // Normalized vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + VSMALL;

            if (magEPrev > SMALL && magE10 > SMALL)
            {
                vector edgeNormal = ePrev ^ e10;
                scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormals[faceI]) < SMALL)
                    {
                        if (faceI != errorFaceI)
                        {
                            // Count only one error per face.
                            errorFaceI = faceI;
                            nConcave++;
                        }

                        if (setPtr)
                        {
                            setPtr->insert(faceI);
                        }

                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if (report)
    {
        if (maxEdgeSin > SMALL)
        {
            scalar maxConcaveDegr =
                Foam::asin(Foam::min(1.0, maxEdgeSin))
             * 180.0/mathematicalConstant::pi;

            Info<< "There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = " << maxConcaveDegr
                << " degrees.\n" << endl;
        }
        else
        {
            Info<< "All angles in faces are convex or less than "  << maxDeg
                << " degrees concave.\n" << endl;
        }
    }

    if (nConcave > 0)
    {
        WarningIn
        (
            "primitiveMesh::checkFaceAngles"
            "(const bool, const scalar, labelHashSet*)"
        )   << nConcave  << " face points with severe concave angle (> "
            << maxDeg << " deg) found.\n"
            << endl;

        return true;
    }
    else
    {
        return false;
    }
}


// Check warpage of faces. Is calculated as the difference between areas of
// individual triangles and the overall area of the face (which ifself is
// is the average of the areas of the individual triangles).
bool primitiveMesh::checkFaceFlatness
(
    const bool report,
    const scalar warnFlatness,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceFlatness"
            << "(const bool, const scalar, labelHashSet*) const: "
            << "checking face flatness" << endl;
    }

    if (warnFlatness < 0 || warnFlatness > 1)
    {
        FatalErrorIn
        (
            "primitiveMesh::checkFaceFlatness"
            "(const bool, const scalar, labelHashSet*)"
        )   << "warnFlatness should be [0..1] but is now " << warnFlatness
            << abort(FatalError);
    }


    const pointField& p = points();
    const faceList& fcs = faces();
    const pointField& fctrs = faceCentres();

    // Areas are calculated as the sum of areas. (see
    // primitiveMeshFaceCentresAndAreas.C)
    scalarField magAreas(mag(faceAreas()));

    label nWarped = 0;

    scalar minFlatness = GREAT;
    scalar sumFlatness = 0;
    label nSummed = 0;

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        if (f.size() > 3 && magAreas[faceI] > VSMALL)
        {
            const point& fc = fctrs[faceI];

            // Calculate the sum of magnitude of areas and compare to magnitude
            // of sum of areas.

            scalar sumA = 0.0;

            forAll(f, fp)
            {
                const point& thisPoint = p[f[fp]];
                const point& nextPoint = p[f.nextLabel(fp)];

                // Triangle around fc.
                vector n = 0.5*((nextPoint - thisPoint)^(fc - thisPoint));
                sumA += mag(n);
            }

            scalar flatness = magAreas[faceI] / (sumA+VSMALL);

            sumFlatness += flatness;
            nSummed++;

            minFlatness = min(minFlatness, flatness);

            if (flatness < warnFlatness)
            {
                nWarped++;

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }
            }
        }
    }


    reduce(nWarped, sumOp<label>());
    reduce(minFlatness, minOp<scalar>());

    reduce(nSummed, sumOp<label>());
    reduce(sumFlatness, sumOp<scalar>());

    if (report)
    {
        if (nSummed > 0)
        {
            Info<< "Face flatness (1 = flat, 0 = butterfly) : average = "
                << sumFlatness / nSummed << "  min = " << minFlatness << endl;
        }

        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with ratio between projected and actual area < "
                << warnFlatness
                << ".\nMinimum ratio (minimum flatness, maximum warpage) = "
                << minFlatness << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the ratio between projected"
                << " and actual area is > " << warnFlatness << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        WarningIn
        (
            "primitiveMesh::checkFaceFlatness"
            "(const bool, const scalar, labelHashSet*)"
        )   << nWarped  << " faces with severe warpage (flatness < "
            << warnFlatness << ") found.\n"
            << endl;

        return true;
    }
    else
    {
        return false;
    }
}


bool primitiveMesh::checkUpperTriangular
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkUpperTriangular("
            << "const bool, labelHashSet*) const: "
            << "checking face ordering" << endl;
    }

    // Check whether internal faces are ordered in the upper triangular order
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    const cellList& c = cells();
    const labelListList& cc = cellCells();

    label internal = nInternalFaces();

    labelList checkInternalFaces(internal, -1);

    label nChecks = 0;

    bool error = false;

    // Loop through faceCells once more and make sure that for internal cell
    // the first label is smaller
    for (label faceI = 0; faceI < internal; faceI++)
    {
        if (own[faceI] >= nei[faceI])
        {
            if (debug || report)
            {
                Pout<< "bool primitiveMesh::checkUpperTriangular("
                    << "const bool, labelHashSet*) const : " << endl
                    << "face " << faceI
                    << " has the owner label greater than neighbour:" << endl
                    << own[faceI] << tab << nei[faceI] << endl;
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            error  = true;
        }
    }

    // Loop through all cells. For each cell, find the face that is internal and
    // add it to the check list (upper triangular order).
    // Once the list is completed, check it against the faceCell list

    forAll (c, cellI)
    {
        const labelList& curFaces = c[cellI];

        // Using the fact that cell neighbour always appear
        // in the increasing order
        const labelList& curNbrs = cc[cellI];

        boolList usedNbr(curNbrs.size(), false);

        for (label nSweeps = 0; nSweeps < curNbrs.size(); nSweeps++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = c.size();

            forAll (curNbrs, nbrI)
            {
                if
                (
                    curNbrs[nbrI] > cellI
                 && !usedNbr[nbrI]
                 && curNbrs[nbrI] < minNei
                )
                {
                    nextNei = nbrI;
                    minNei = curNbrs[nbrI];
                }
            }

            if (nextNei > -1)
            {
                // Mark this neighbour as used
                usedNbr[nextNei] = true;

                forAll (curFaces, faceI)
                {
                    if (curFaces[faceI] < internal)
                    {
                        if (nei[curFaces[faceI]] == curNbrs[nextNei])
                        {
                            checkInternalFaces[nChecks] = curFaces[faceI];
                            nChecks++;

                            break;
                        }
                    }
                }
            }
        }
    }

    // Check list created. If everything is OK, the face label is equal to index
    forAll (checkInternalFaces, faceI)
    {
        if (checkInternalFaces[faceI] != faceI)
        {
            error = true;

            Pout<< "bool primitiveMesh::checkUpperTriangular(const bool"
                << ", labelHashSet*) const : " << endl
                << "face " << faceI << " out of position. Markup label: "
                << checkInternalFaces[faceI] << ". All subsequent faces will "
                << "also be out of position. Please check the mesh manually."
                << endl;

            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            break;
        }
    }

    reduce(error, orOp<bool>());

    if (error)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkUpperTriangular(const bool"
            ", labelHashSet*) const"
        )   << "Error in face ordering: faces not in upper triangular order!"
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Upper triangular ordering OK.\n" << endl;
        }

        return false;
    }
}


bool primitiveMesh::checkCellsZipUp
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellsZipUp("
            << "const bool, labelHashSet*) const: "
            << "checking topological cell openness" << endl;
    }

    label nOpenCells = 0;

    const faceList& f = faces();

    const cellList& c = cells();

    forAll (c, cellI)
    {
        const labelList& curFaces = c[cellI];

        const edgeList cellEdges = c[cellI].edges(f);

        labelList edgeUsage(cellEdges.size(), 0);

        forAll (curFaces, faceI)
        {
            edgeList curFaceEdges = f[curFaces[faceI]].edges();

            forAll (curFaceEdges, faceEdgeI)
            {
                const edge& curEdge = curFaceEdges[faceEdgeI];

                forAll (cellEdges, cellEdgeI)
                {
                    if (cellEdges[cellEdgeI] == curEdge)
                    {
                        edgeUsage[cellEdgeI]++;
                        break;
                    }
                }
            }
        }

        edgeList singleEdges(cellEdges.size());
        label nSingleEdges = 0;

        forAll (edgeUsage, edgeI)
        {
            if (edgeUsage[edgeI] == 1)
            {
                singleEdges[nSingleEdges] = cellEdges[edgeI];
                nSingleEdges++;
            }
            else if (edgeUsage[edgeI] != 2)
            {
                WarningIn
                (
                    "bool primitiveMesh::checkCellsZipUp("
                    "const bool, labelHashSet*) const"
                )   << "edge " << cellEdges[edgeI] << " in cell " << cellI
                    << " used " << edgeUsage[edgeI] << " times. " << endl
                    << "Should be 1 or 2 - serious error in mesh structure"
                    << endl;

                if (setPtr)
                {
                    setPtr->insert(cellI);
                }
            }
        }

        if (nSingleEdges > 0)
        {
            if (debug || report)
            {
                singleEdges.setSize(nSingleEdges);

                Pout<< "bool primitiveMesh::checkCellsZipUp(const bool"
                    << ", labelHashSet*) const : " << endl
                    << "Cell " << cellI << " has got " << nSingleEdges
                    << " unmatched edges: " << singleEdges << endl;
            }

            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nOpenCells++;
        }
    }

    reduce(nOpenCells, sumOp<label>());

    if (nOpenCells > 0)
    {
        WarningIn
        (
            "bool primitiveMesh::checkCellsZipUp("
            "const bool, labelHashSet*) const"
        )   << nOpenCells
            << " open cells found.  Please use the mesh zip-up tool. "
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Topological cell zip-up check OK.\n" << endl;
        }

        return false;
    }
}


// Vertices of face within point range and unique.
bool primitiveMesh::checkFaceVertices
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceVertices("
            << "const bool, labelHashSet*) const: "
            << "checking face vertices" << endl;
    }

    // Check that all vertex labels are valid
    const faceList& f = faces();

    label nErrorFaces = 0;

    forAll (f, fI)
    {
        const face& curFace = f[fI];

        if (min(curFace) < 0 || max(curFace) > nPoints())
        {
            WarningIn
            (
                "bool primitiveMesh::checkFaceVertices("
                "const bool, labelHashSet*) const"
            )   << "Face " << fI << " contains vertex labels out of range: "
                << curFace << " Max point index = " << nPoints()-1 << endl;

            if (setPtr)
            {
                setPtr->insert(fI);
            }

            nErrorFaces++;
        }

        // Uniqueness of vertices
        labelHashSet facePoints(2*curFace.size());

        forAll(curFace, fp)
        {
            bool inserted = facePoints.insert(curFace[fp]);

            if (!inserted)
            {
                WarningIn
                (
                    "bool primitiveMesh::checkFaceVertices("
                    "const bool, labelHashSet*) const"
                )   << "Face " << fI << " contains duplicate vertex labels: "
                    << curFace << endl;

                if (setPtr)
                {
                    setPtr->insert(fI);
                }

                nErrorFaces++;
            }
        }
    }

    reduce(nErrorFaces, sumOp<label>());

    if (nErrorFaces > 0)
    {
        SeriousErrorIn
        (
            "bool primitiveMesh::checkFaceVertices("
            "const bool, labelHashSet*) const"
        )   << "const bool, labelHashSet*) const: "
            << nErrorFaces << " faces with invalid vertex labels found"
            << endl;

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Face vertices OK.\n" << endl;
        }

        return false;
    }
}


// Warn if more than 2 points shared between faces.
void primitiveMesh::warnCommonPoints
(
    const label faceI,
    const Map<label>& nCommonPoints,
    bool& hasWarned
) const
{
    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();

        label nCommon = iter();

        if (nCommon >= 3)
        {
            // Warning for if more than 2 shared points
            if (!hasWarned)
            {
                WarningIn
                (
                    "bool primitiveMesh::warnCommonPoints(const label"
                    ", const Map<label>&, bool&) const"
                )   << "Face " << faceI << " vertex labels "
                    << faces()[faceI]
                    << " and face " << nbFaceI << " vertex labels "
                    << faces()[nbFaceI]
                    << " share " << nCommon << " vertices."
                    << endl
                    << "This mesh should only be used with "
                    << "faceDecomp finiteElement decomposition."
                    << endl
                    << "Further face-face warnings will be"
                    << " suppressed." << endl;

                hasWarned = true;
            }
        }
    }
}


// Check if all points on face are shared between faces.
bool primitiveMesh::checkDuplicateFaces
(
    const label faceI,
    const Map<label>& nCommonPoints,
    labelHashSet* setPtr
) const
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();

        label nCommon = iter();

        const face& curFace = faces()[faceI];
        const face& nbFace = faces()[nbFaceI];

        if (nCommon == nbFace.size() || nCommon == curFace.size())
        {
            Pout<< "bool primitiveMesh::checkDuplicateFaces("
                << "const label, const Map<label>&) const :" << endl
                << "Face " << faceI << " vertex labels " << curFace
                << " and face " << nbFaceI << " vertex labels " << nbFace
                << " share too many vertices" << endl;

            if (setPtr)
            {
                setPtr->insert(faceI);
                setPtr->insert(nbFaceI);
            }

            error = true;
        }
    }
    return error;
}


// Check that shared points are in consecutive order.
bool primitiveMesh::checkCommonOrder
(
    const label faceI,
    const Map<label>& nCommonPoints,
    labelHashSet* setPtr
) const
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();

        label nCommon = iter();

        const face& curFace = faces()[faceI];
        const face& nbFace = faces()[nbFaceI];

        if
        (
            nCommon >= 2
         && nCommon != nbFace.size()
         && nCommon != curFace.size()
        )
        {
            forAll(curFace, fp)
            {
                // Get the index in the neighbouring face shared with curFace
                label nb = findIndex(nbFace, curFace[fp]);

                if (nb != -1)
                {

                    // Check the whole face from nb onwards for shared vertices
                    // with neighbouring face. Rule is that any shared vertices
                    // should be consecutive on both faces i.e. if they are
                    // vertices fp,fp+1,fp+2 on one face they should be vertices
                    // nb, nb+1, nb+2 (or nb+2, nb+1, nb) on the other face.


                    // Vertices before and after on curFace
                    label fpPlus1 = (fp+1) % curFace.size();
                    label fpMin1 = (fp == 0 ? curFace.size()-1 : fp-1);

                    // Vertices before and after on nbFace
                    label nbPlus1 = (nb+1) % nbFace.size();
                    label nbMin1 = (nb == 0 ? nbFace.size()-1 : nb-1);

                    // Find order of walking by comparing next points on both
                    // faces.
                    label curInc = labelMax;
                    label nbInc = labelMax;

                    if (nbFace[nbPlus1] == curFace[fpPlus1])
                    {
                        curInc = 1;
                        nbInc = 1;
                    }
                    else if (nbFace[nbPlus1] == curFace[fpMin1])
                    {
                        curInc = -1;
                        nbInc = 1;
                    }
                    else if (nbFace[nbMin1] == curFace[fpMin1])
                    {
                        curInc = -1;
                        nbInc = -1;
                    }
                    else
                    {
                        curInc = 1;
                        nbInc = -1;
                    }


                    // Pass1: loop until start of common vertices found.
                    label curNb = nb;
                    label curFp = fp;

                    do
                    {
                        curFp += curInc;

                        if (curFp >= curFace.size())
                        {
                            curFp = 0;
                        }
                        else if (curFp < 0)
                        {
                            curFp = curFace.size()-1;
                        }

                        curNb += nbInc;

                        if (curNb >= nbFace.size())
                        {
                            curNb = 0;
                        }
                        else if (curNb < 0)
                        {
                            curNb = nbFace.size()-1;
                        }
                    } while (curFace[curFp] == nbFace[curNb]);


                    // Pass2: check equality walking from curFp, curNb
                    // in opposite order.

                    curInc = -curInc;
                    nbInc = -nbInc;

                    for (label commonI = 0; commonI < nCommon; commonI++)
                    {
                        curFp += curInc;

                        if (curFp >= curFace.size())
                        {
                            curFp = 0;
                        }
                        else if (curFp < 0)
                        {
                            curFp = curFace.size()-1;
                        }

                        curNb += nbInc;

                        if (curNb >= nbFace.size())
                        {
                            curNb = 0;
                        }
                        else if (curNb < 0)
                        {
                            curNb = nbFace.size()-1;
                        }

                        if (curFace[curFp] != nbFace[curNb])
                        {
                            Pout<< "bool primitiveMesh::checkCommonOrder("
                                << "const label, const Map<label>&) const: "
                                << endl
                                << "Face " << faceI << " vertex labels "
                                << curFace << " and face " << nbFaceI
                                << " with vertex labels " << nbFace <<" share "
                                << nCommon
                                << " vertices which are not in consecutive"
                                << " order" << endl;

                            if (setPtr)
                            {
                                setPtr->insert(faceI);
                                setPtr->insert(nbFaceI);
                            }

                            error = true;

                            break;
                        }
                    }


                    // Done the curFace - nbFace combination.
                    break;
                }
            }
        }
    }

    return error;
}


// Checks common vertices between faces. If more than 2 they should be
// consecutive on both faces.
bool primitiveMesh::checkFaceFaces
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceFaces(const bool, labelHashSet*)"
            << " const: " << "checking face-face connectivity" << endl;
    }

    const labelListList& pf = pointFaces();

    label nErrorDuplicate = 0;
    label nErrorOrder = 0;

    // Whether to print warnings.
    bool hasWarned = false;

    for (label faceI = 0; faceI < nFaces(); faceI++)
    {
        const face& curFace = faces()[faceI];

        // Calculate number of common points between current faceI and
        // neighbouring face. Store on map.

        Map<label> nCommonPoints(6*curFace.size());

        forAll(curFace, fp)
        {
            label pointI = curFace[fp];

            const labelList& nbs = pf[pointI];

            forAll(nbs, nbI)
            {
                label nbFaceI = nbs[nbI];

                if (faceI < nbFaceI)
                {
                    // Only check once for each combination of two faces.

                    Map<label>::iterator fnd = nCommonPoints.find(nbFaceI);

                    if (fnd == nCommonPoints.end())
                    {
                        // First common vertex found.
                        nCommonPoints.insert(nbFaceI, 1);
                    }
                    else
                    {
                        fnd()++;
                    }
                }
            }
        }

        // Perform various checks on common points

        // Warn if more than 2 points shared
        warnCommonPoints(faceI, nCommonPoints, hasWarned);

        // Check all vertices shared (duplicate point)
        if (checkDuplicateFaces(faceI, nCommonPoints, setPtr))
        {
            nErrorDuplicate++;
        }

        // Check common vertices are consecutive on both faces
        if (checkCommonOrder(faceI, nCommonPoints, setPtr))
        {
            nErrorOrder++;
        }
    }

    reduce(nErrorDuplicate, sumOp<label>());
    reduce(nErrorOrder, sumOp<label>());

    if (nErrorDuplicate > 0 || nErrorOrder > 0)
    {
        if (nErrorDuplicate > 0)
        {
            SeriousErrorIn
            (
                "bool primitiveMesh::checkFaceFaces(const bool"
                ", labelHashSet*) const"
            )   << nErrorDuplicate << " duplicate faces found" << endl;
        }
        if (nErrorOrder > 0)
        {
            SeriousErrorIn
            (
                "bool primitiveMesh::checkFaceFaces(const bool"
                ", labelHashSet*) const"
            )   << nErrorOrder
                << " faces with non-consecutive shared points found" << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Face-face connectivity OK.\n" << endl;
        }

        return false;
    }
}


// Checks cells with 1 or less internal faces. Give numerical problems.
bool primitiveMesh::checkFloatingCells
(
    const bool report,    // report,
    labelHashSet* setPtr  // setPtr
) const
{
    /*
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFloatingCells(const bool"
            << ", labelHashSet*) const: "
            << "checking floating cells" << endl;
    }

    const cellList& c = cells();
    
    label nErrorCells = 0;

    scalar minDet = GREAT;
    scalar sumDet = 0;
    label nSummed = 0;

    forAll (c, cellI)
    {
        const labelList& curFaces = c[cellI];

        // Calculate local normalization factor
        scalar avgArea = 0;

        label nInternalFaces = 0;

        forAll(curFaces, i)
        {
            if (isInternalFace(curFaces[i]))
            {
                avgArea += mag(faceAreas()[curFaces[i]]);

                nInternalFaces++;
            }
        }

        if (nInternalFaces == 0)
        {
            if (debug || report)
            {
                Pout<< " Cell " << cellI << " has only " << nInternalFaces
                    << " cells connected other cells. This will probably lead"
                    << " to numerical problems" << endl;
            }

            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nErrorCells++;
        }
        else
        {
            avgArea /= nInternalFaces;

            tensor areaTensor(tensor::zero);

            forAll(curFaces, i)
            {
                if (isInternalFace(curFaces[i]))
                {
                    areaTensor += sqr(faceAreas()[curFaces[i]]/avgArea);
                }
            }

            scalar determinant = mag(det(areaTensor));

            minDet = min(determinant, minDet);
            sumDet += determinant;
            nSummed++;

            if (determinant < 1E-3)
            {
                Pout<< "cell:" << cellI << " det:" << determinant
                    << " avgArea:" << avgArea << endl;

                forAll(curFaces, i)
                {
                    Pout<< "    face:" << curFaces[i]
                        << " norm.area:" << faceAreas()[curFaces[i]]/avgArea
                        << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(cellI);
                }

                nErrorCells++;
            }
        }
    }

    reduce(nErrorCells, sumOp<label>());
    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    reduce(nSummed, sumOp<label>());

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "Cell determinant (wellposedness) : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorCells > 0)
    {
        WarningIn
        (
            "bool primitiveMesh::checkFloatingCells(const bool"
            ", labelHashSet*) const"
        )   << nErrorCells << " cells which have a small cell determinant"
            << " found." << endl
            << "    This cell determinant is a measure for how well a gradient"
            << " can be determined." << endl
            << "    The cells with a low value might"
            << " give numerical problems in a 3D geometry if" << endl
            << "    they are inside"
            << " the domain or if they are used with non-fixed value boundary"
            << endl << "    conditions"
            << endl;

        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "Floating cell check OK.\n" << endl;
        }

        return false;
    }
    */
    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool primitiveMesh::checkGeometry(const bool report) const
{
    label noFailedCheckes = 0;

    if (checkClosedBoundary(report)) noFailedCheckes++;
    if (checkClosedCells(report)) noFailedCheckes++;
    if (checkFaceAreas(report)) noFailedCheckes++;
    if (checkCellVolumes(report)) noFailedCheckes++;
    if (checkFaceDotProduct(report)) noFailedCheckes++;
    if (checkFacePyramids(report)) noFailedCheckes++;
    if (checkFaceSkewness(report)) noFailedCheckes++;

    if (noFailedCheckes == 0)
    {
        if (debug || report)
        {
            Info << "Mesh geometry OK." << endl;
        }
        return false;
    }
    else
    {
        Info<< "Failed " << noFailedCheckes << " mesh geometry checks." << endl;

        return true;
    }
}


bool primitiveMesh::checkTopology(const bool report) const
{
    label noFailedCheckes = 0;

    if (checkPoints(report)) noFailedCheckes++;
    if (checkUpperTriangular(report)) noFailedCheckes++;
    if (checkCellsZipUp(report)) noFailedCheckes++;
    if (checkFaceVertices(report)) noFailedCheckes++;
    if (checkFaceFaces(report)) noFailedCheckes++;
    if (checkFloatingCells(report)) noFailedCheckes++;

    if (noFailedCheckes == 0)
    {
        if (debug || report)
        {
            Info << "Mesh topology OK." << endl;
        }

        return false;
    }
    else
    {
        Info<< "Failed " << noFailedCheckes << " mesh topology checks." << endl;

        return true;
    }
}


bool primitiveMesh::checkMesh(const bool report) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkMesh(const bool report) const: "
            << "checking primitiveMesh" << endl;
    }

    bool failedCheckes = checkTopology(report) | checkGeometry(report);

    if (!failedCheckes)
    {
        if (debug || report)
        {
            Info << "Mesh OK." << endl;
        }

        return false;
    }
    else
    {
        Info << "Failed some mesh checks." << endl;

        return true;
    }
}


scalar primitiveMesh::setOrthWarn(const scalar ow)
{
    scalar oldOrthWarn = orthWarn_;
    orthWarn_ = ow;

    return oldOrthWarn;
}


scalar primitiveMesh::setSkewWarn(const scalar sw)
{
    scalar oldSkewWarn = skewWarn_;
    skewWarn_ = sw;

    return oldSkewWarn;
}


scalar primitiveMesh::setAspectWarn(const scalar aw)
{
    scalar oldAspectWarn = aspectWarn_;
    aspectWarn_ = aw;

    return oldAspectWarn;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
