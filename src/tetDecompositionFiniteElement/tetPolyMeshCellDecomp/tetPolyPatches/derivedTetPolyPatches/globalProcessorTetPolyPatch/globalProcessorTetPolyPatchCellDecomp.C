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

\*---------------------------------------------------------------------------*/

#include "globalProcessorTetPolyPatchCellDecomp.H"
#include "tetPolyMeshCellDecomp.H"
#include "demandDrivenData.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(globalProcessorTetPolyPatchCellDecomp, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
globalProcessorTetPolyPatchCellDecomp::globalProcessorTetPolyPatchCellDecomp
(
    const label globalPointSize,
    const labelList& meshPoints,
    const labelList& sharedPointAddr,
    const label globalEdgeSize,
    const edgeList& meshEdges,
    const labelList& sharedEdgeAddr,
    const edgeList& meshCutEdges,
    const scalarField& meshCutEdgeMask,
    const tetPolyBoundaryMeshCellDecomp& bm,
    const label index
)
:
    tetPolyPatchCellDecomp(bm),
    globalPointSize_(globalPointSize),
    meshPoints_(meshPoints),
    sharedPointAddr_(sharedPointAddr),
    globalEdgeSize_(globalEdgeSize),
    meshEdges_(meshEdges),
    sharedEdgeAddr_(sharedEdgeAddr),
    meshCutEdges_(meshCutEdges),
    meshCutEdgeMask_(meshCutEdgeMask),
    boundaryIndex_(index),
    localEdgeIndicesPtr_(NULL),
    cutEdgeIndicesPtr_(NULL),
    cutEdgeOwnerIndicesPtr_(NULL),
    cutEdgeOwnerStartPtr_(NULL),
    cutEdgeNeighbourIndicesPtr_(NULL),
    cutEdgeNeighbourStartPtr_(NULL),
    doubleCutEdgeIndicesPtr_(NULL),
    doubleCutOwnerPtr_(NULL),
    doubleCutNeighbourPtr_(NULL),
    ownNeiDoubleMaskPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

globalProcessorTetPolyPatchCellDecomp::~globalProcessorTetPolyPatchCellDecomp()
{
    deleteDemandDrivenData(localEdgeIndicesPtr_);

    clearCutEdgeAddressing();
}


void globalProcessorTetPolyPatchCellDecomp::clearCutEdgeAddressing() const
{
    deleteDemandDrivenData(cutEdgeIndicesPtr_);
    deleteDemandDrivenData(cutEdgeOwnerIndicesPtr_);
    deleteDemandDrivenData(cutEdgeOwnerStartPtr_);
    deleteDemandDrivenData(cutEdgeNeighbourIndicesPtr_);
    deleteDemandDrivenData(cutEdgeNeighbourStartPtr_);

    deleteDemandDrivenData(doubleCutEdgeIndicesPtr_);
    deleteDemandDrivenData(doubleCutOwnerPtr_);
    deleteDemandDrivenData(doubleCutNeighbourPtr_);

    deleteDemandDrivenData(ownNeiDoubleMaskPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointField& globalProcessorTetPolyPatchCellDecomp::localPoints() const
{
    notImplemented("processorTetPolyPatchCellDecomp::localPoints() const");
    return Field<point>::null();
}


const vectorField& globalProcessorTetPolyPatchCellDecomp::pointNormals() const
{
    notImplemented("processorTetPolyPatchCellDecomp::pointNormals() const");
    return Field<point>::null();
}


triFaceList globalProcessorTetPolyPatchCellDecomp::faceTriangles
(
    const label faceID
) const
{
    notImplemented
    (
        "processorTetPolyPatchCellDecomp::faceTriangles(label faceID) const"
    );

    return List<triFace>::null();
}


faceList globalProcessorTetPolyPatchCellDecomp::triFaces() const
{
    notImplemented
    (
        "faceList processorTetPolyPatchCellDecomp::triFaces() const"
    );

    return faceList::null();
}


void globalProcessorTetPolyPatchCellDecomp::updateMesh()
{
    notImplemented("processorTetPolyPatchFaceDecomp::updateMesh()");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
