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

\*---------------------------------------------------------------------------*/

#include "processorFaPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "IPstream.H"
#include "OPstream.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorFaPatch, 0);
addToRunTimeSelectionTable(faPatch, processorFaPatch, dictionary);


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void processorFaPatch::initGeometry()
{
    if (Pstream::parRun())
    {
        OPstream toNeighbProc(neighbProcNo(), 3*size()*sizeof(vector));

        toNeighbProc
            << edgeCentres()
            << edgeLengths()
            << edgeFaceCentres();
    }
}


void processorFaPatch::calcGeometry()
{
    if (Pstream::parRun())
    {
        {
            IPstream fromNeighbProc(neighbProcNo(), 3*size()*sizeof(vector));
            fromNeighbProc
                >> neighbEdgeCentres_
                >> neighbEdgeLengths_
                >> neighbEdgeFaceCentres_;
        }

        const scalarField& magEl = magEdgeLengths();

        forAll(magEl, edgei)
        {
            scalar nmagEl = mag(neighbEdgeLengths_[edgei]);
            scalar avEl = (magEl[edgei] + nmagEl)/2.0;

            if (mag(magEl[edgei] - nmagEl)/avEl > 1e-6)
            {
                FatalErrorIn
                (
                    "processorFvPatch::makeWeights(scalarField& w) const"
                )   << "edge " << edgei
                    << " length does not match neighbour by "
                    << 100*mag(magEl[edgei] - nmagEl)/avEl
                    << "% -- possible edge ordering problem"
                    << exit(FatalError);
            }
        }

        calcTransformTensors
        (
            edgeCentres(),
            neighbEdgeCentres_,
            edgeNormals(),
            neighbEdgeLengths_/mag(neighbEdgeLengths_)
        );
    }
}


void processorFaPatch::initMovePoints(const pointField& p)
{
    faPatch::movePoints(p);
    initGeometry();
}


void processorFaPatch::movePoints(const pointField&)
{
    calcGeometry();
}


// Make patch weighting factors
void processorFaPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        // The face normals point in the opposite direction on the other side
        scalarField neighbEdgeCentresCn
        (
            (
                neighbEdgeLengths()
               /mag(neighbEdgeLengths())
            )
          & (
              neighbEdgeCentres()
            - neighbEdgeFaceCentres())
        );

        w = neighbEdgeCentresCn/
            (
                (edgeNormals() & faPatch::delta())
              + neighbEdgeCentresCn
            );
    }
    else
    {
        w = 1.0;
    }
}


// Make patch edge - neighbour face distances
void processorFaPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (Pstream::parRun())
    {
        dc = (1.0 - weights())/(edgeNormals() & faPatch::delta());
    }
    else
    {
        dc = 1.0/(edgeNormals() & faPatch::delta());
    }
}


// Return delta (P to N) vectors across coupled patch
tmp<vectorField> processorFaPatch::delta() const
{
    if (Pstream::parRun())
    {
        // To the transformation if necessary
        if (parallel())
        {
            return
                faPatch::delta()
              - (
                    neighbEdgeCentres()
                  - neighbEdgeFaceCentres()
                );
        }
        else
        {
            return
                faPatch::delta()
              - transform
                (
                    forwardT(),
                    (
                        neighbEdgeCentres()
                      - neighbEdgeFaceCentres()
                    )
                );
        }
    }
    else
    {
        return faPatch::delta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
