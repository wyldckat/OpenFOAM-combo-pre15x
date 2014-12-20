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
    Calculate cut edge globalProcessor addressing, needed for vector-matrix
    multiply on globalProcessor boundaries.

\*---------------------------------------------------------------------------*/

#include "globalProcessorTetPolyPatchCellDecomp.H"
#include "tetPolyBoundaryMeshCellDecomp.H"
#include "tetPolyMeshCellDecomp.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void globalProcessorTetPolyPatchCellDecomp::calcLocalEdgesIndices() const
{
    if (debug)
    {
        Info<< "labelList globalProcessorTetPolyPatchCellDecomp::"
            << "calcLocalEdgesIndices() const : "
            << "calculating local edge indices"
            << endl;
    }

    // Get reference to the mesh
    const tetPolyMeshCellDecomp& mesh = boundaryMesh().mesh();

    // Get reference to edges
    const edgeList& patchEdges = meshEdges();

    localEdgeIndicesPtr_ = new labelList(patchEdges.size(), -1);
    labelList& localEdgeInd = *localEdgeIndicesPtr_;

    const lduAddressing& lduAddr = mesh.ldu();

    forAll (patchEdges, edgeI)
    {
        localEdgeInd[edgeI] =
            lduAddr.triIndex
            (
                patchEdges[edgeI].start(),
                patchEdges[edgeI].end()
            );
    }

#   ifdef DEBUGtetFemMatrix
    if (localEdgeInd.size() > 0 && min(localEdgeInd) < 0)
    {
        FatalErrorIn
        (
            "void globalProcessorTetPolyPatchCellDecomp::"
            "calcLocalEdgesIndices() const"
        )   << "Problem in local edge addressing"
            << abort(FatalError);
    }
#   endif

    if (debug)
    {
        Info<< "void globalProcessorTetPolyPatchCellDecomp::"
            << "calcLocalEdgesIndices() const : "
            << "finished calculating local edge indices"
            << endl;
    }
}


void globalProcessorTetPolyPatchCellDecomp::calcCutEdgeIndices() const
{
    if (debug)
    {
        Info<< "void globalProcessorTetPolyPatchCellDecomp::"
            << "calcCutEdgeIndices() const : "
            << "calculating cut edge indices"
            << endl;
    }

    if (cutEdgeIndicesPtr_)
    {
        FatalErrorIn
        (
            "void globalProcessorTetPolyPatchCellDecomp::"
            "calcCutEdgesIndices() const"
        )   << "addressing already allocated"
            << abort(FatalError);
    }

    // Get reference to the mesh
    const tetPolyMeshCellDecomp& mesh = boundaryMesh().mesh();

    // Get reference to edges
    const edgeList& patchCutEdges = meshCutEdges();

    cutEdgeIndicesPtr_ = new labelList(patchCutEdges.size(), -1);
    labelList& cutEdgeInd = *cutEdgeIndicesPtr_;

    const lduAddressing& lduAddr = mesh.ldu();

    forAll (patchCutEdges, edgeI)
    {
        cutEdgeInd[edgeI] =
            lduAddr.triIndex
            (
                patchCutEdges[edgeI].start(),
                patchCutEdges[edgeI].end()
            );
    }

    if (debug)
    {
        Info<< "void globalProcessorTetPolyPatchCellDecomp::"
            << "calcCutEdgeIndices() const : "
            << "finished calculating cut edge indices"
            << endl;
    }
}


void globalProcessorTetPolyPatchCellDecomp::calcCutEdgeAddressing() const
{
    if
    (
        cutEdgeOwnerIndicesPtr_
     || cutEdgeOwnerStartPtr_
     || cutEdgeNeighbourIndicesPtr_
     || cutEdgeNeighbourStartPtr_
     || doubleCutEdgeIndicesPtr_
     || doubleCutOwnerPtr_
     || doubleCutNeighbourPtr_
     || ownNeiDoubleMaskPtr_
    )
    {
        FatalErrorIn
        (
            "void globalProcessorTetPolyPatchCellDecomp::"
            "calcCutEdgeAddressing() "
            "const"
        )   << "addressing already allocated"
            << abort(FatalError);
    }


    // Make a list over all edges in the mesh.  Mark the ones that are local
    // to the patch and then collect the rest.
    // For doubly cut edges, mark up the local points

    const tetPolyMeshCellDecomp& mesh = boundaryMesh().mesh();

    // Get reference to mesh points
    const labelList& mp = meshPoints();

    // Get reference to local edge indices
    const labelList& cutEdges = cutEdgeIndices();

    // Get the cut edge mask.  It will be used to create the new mask
    const scalarField& cutMask = meshCutEdgeMask();

    // Mark up the cut edges

    labelList cutIndexLookup(mesh.nEdges(), -1);

    forAll (cutEdges, edgeI)
    {
        cutIndexLookup[cutEdges[edgeI]] = edgeI;
    }

    labelList localPointLabel(mesh.nPoints(), -1);

    forAll (mp, pointI)
    {
        localPointLabel[mp[pointI]] = pointI;
    }

    // Get reference to addressing
    const lduAddressing& ldu = mesh.ldu();

    const labelList& globalOwner = ldu.lowerAddr();
    const labelList& globalNeighbour = ldu.upperAddr();

    // Allocate doubly cut arrays
    doubleCutEdgeIndicesPtr_ = new labelList(cutEdges.size(), -1);
    labelList& doubleCutEdges = *doubleCutEdgeIndicesPtr_;

    doubleCutOwnerPtr_ = new labelList(cutEdges.size(), -1);
    labelList& doubleCutOwn = *doubleCutOwnerPtr_;

    doubleCutNeighbourPtr_ = new labelList(cutEdges.size(), -1);
    labelList& doubleCutNei = *doubleCutNeighbourPtr_;

    label nDoubleCut = 0;

    // Owner side
    //~~~~~~~~~~~

    // Allocate the array
    cutEdgeOwnerIndicesPtr_ = new labelList(cutEdges.size(), -1);
    labelList& own = *cutEdgeOwnerIndicesPtr_;
    label nOwn = 0;

    cutEdgeOwnerStartPtr_ = new labelList(mp.size() + 1, -1);
    labelList& ownStart = *cutEdgeOwnerStartPtr_;

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as cut;
    // add it to the list of cut edges and grab the mask value that goes
    // with it.
    forAll (mp, pointI)
    {
        ownStart[pointI] = nOwn;

        const label curPointID = mp[pointI];

        // get owner edges indices
        const label startFaceOwn = ldu.ownerStartAddr()[curPointID];
        const label endFaceOwn = ldu.ownerStartAddr()[curPointID + 1];

        for
        (
            label edgeLabel = startFaceOwn;
            edgeLabel < endFaceOwn;
            edgeLabel++
        )
        {
            if (cutIndexLookup[edgeLabel] >= 0)
            {
                if (localPointLabel[globalNeighbour[edgeLabel]] == -1)
                {
                    // Singly cut edge
                    own[nOwn] = edgeLabel;
                    nOwn++;
                }
                else
                {
                    // Doubly cut edge
                    doubleCutEdges[nDoubleCut] = edgeLabel;
                    doubleCutOwn[nDoubleCut] = pointI;
                    doubleCutNei[nDoubleCut] =
                        localPointLabel[globalNeighbour[edgeLabel]];

                    nDoubleCut++;
                }
            }
        }
    }

    // reset the size of owner edges
    own.setSize(nOwn);

    // set the last start label by hand
    ownStart[mp.size()] = nOwn;

    // Neighbour side
    //~~~~~~~~~~~~~~~

    // Allocate the array
    cutEdgeNeighbourIndicesPtr_ = new labelList(cutEdges.size(), -1);
    labelList& nei = *cutEdgeNeighbourIndicesPtr_;
    label nNei = 0;

    cutEdgeNeighbourStartPtr_ = new labelList(mp.size() + 1, -1);
    labelList& neiStart = *cutEdgeNeighbourStartPtr_;

    const unallocLabelList& losort = ldu.losortAddr();

    // Go through all the local points and get all the edges coming
    // from that point.  Check if the edge has been marked as local;
    // if not, add it to the list of cut edges.
    forAll (mp, pointI)
    {
        neiStart[pointI] = nNei;

        const label curPointID = mp[pointI];

        // get neighbour edges indices
        const label startFaceNei = ldu.losortStartAddr()[curPointID];
        const label endFaceNei = ldu.losortStartAddr()[curPointID + 1];

        for
        (
            label edgeLabel = startFaceNei;
            edgeLabel < endFaceNei;
            edgeLabel++
        )
        {
            if (cutIndexLookup[losort[edgeLabel]] >= 0)
            {
                if (localPointLabel[globalOwner[losort[edgeLabel]]] == -1)
                {
                    // Singly cut edge
                    nei[nNei] = losort[edgeLabel];
                    nNei++;
                }
                else
                {
                    // Doubly cut edge
                    doubleCutEdges[nDoubleCut] = losort[edgeLabel];
                    doubleCutOwn[nDoubleCut] =
                        localPointLabel[globalOwner[losort[edgeLabel]]];
                    doubleCutNei[nDoubleCut] = pointI;

                    nDoubleCut++;
                }
            }
        }
    }

    // reset the size of neighbour edges
    nei.setSize(nNei);

    // set the last start label by hand
    neiStart[mp.size()] = nNei;

    // Reset the size of double cut edge data
    doubleCutEdges.setSize(nDoubleCut);
    doubleCutOwn.setSize(nDoubleCut);
    doubleCutNei.setSize(nDoubleCut);

    // Allocate the reordered mask that corresponds to the owner-neighbour
    // coefficient ordering
    ownNeiDoubleMaskPtr_ =
        new scalarField
        (
            own.size() + nei.size()
          + doubleCutOwn.size() + doubleCutNei.size()
        );
    scalarField& ownNeiDoubleMask = *ownNeiDoubleMaskPtr_;

    label nMask = 0;

    forAll (own, i)
    {
        ownNeiDoubleMask[nMask] = cutMask[cutIndexLookup[own[i]]];
        nMask++;
    }

    forAll (nei, i)
    {
        ownNeiDoubleMask[nMask] = cutMask[cutIndexLookup[nei[i]]];
        nMask++;
    }

    forAll (doubleCutOwn, i)
    {
        ownNeiDoubleMask[nMask] = cutMask[cutIndexLookup[doubleCutOwn[i]]];
        nMask++;

        ownNeiDoubleMask[nMask] = cutMask[cutIndexLookup[doubleCutNei[i]]];
        nMask++;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& globalProcessorTetPolyPatchCellDecomp::localEdgeIndices() const
{
    if (!localEdgeIndicesPtr_)
    {
        calcLocalEdgesIndices();
    }

    return *localEdgeIndicesPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::cutEdgeIndices() const
{
    if (!cutEdgeIndicesPtr_)
    {
        calcCutEdgeIndices();
    }

    return *cutEdgeIndicesPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::cutEdgeOwnerIndices() const
{
    if (!cutEdgeOwnerIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeOwnerIndicesPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::cutEdgeOwnerStart() const
{
    if (!cutEdgeOwnerStartPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeOwnerStartPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::cutEdgeNeighbourIndices() const
{
    if (!cutEdgeNeighbourIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeNeighbourIndicesPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::cutEdgeNeighbourStart() const
{
    if (!cutEdgeNeighbourStartPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *cutEdgeNeighbourStartPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::doubleCutEdgeIndices() const
{
    if (!doubleCutEdgeIndicesPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutEdgeIndicesPtr_;
}


const labelList& globalProcessorTetPolyPatchCellDecomp::doubleCutOwner() const
{
    if (!doubleCutOwnerPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutOwnerPtr_;
}


const labelList&
globalProcessorTetPolyPatchCellDecomp::doubleCutNeighbour() const
{
    if (!doubleCutNeighbourPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *doubleCutNeighbourPtr_;
}


const scalarField&
globalProcessorTetPolyPatchCellDecomp::ownNeiDoubleMask() const
{
    if (!ownNeiDoubleMaskPtr_)
    {
        calcCutEdgeAddressing();
    }

    return *ownNeiDoubleMaskPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
