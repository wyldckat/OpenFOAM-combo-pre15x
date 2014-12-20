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

#include "globalProcessorPointPatch.H"
#include "pointMesh.H"
#include "demandDrivenData.H"
#include "parallelInfo.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(globalProcessorPointPatch, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
globalProcessorPointPatch::globalProcessorPointPatch
(
    const parallelInfo& pi,
    const pointBoundaryMesh& bm,
    const label
)
:
    pointPatch(bm),
    globalPointSize_(pi.nGlobalPoints()),
    meshPoints_(pi.sharedPointLabels()),
    sharedPointAddr_(pi.sharedPointAddr())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

globalProcessorPointPatch::~globalProcessorPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointField& globalProcessorPointPatch::localPoints() const
{
    notImplemented("globalProcessorPointPatch::localPoints() const");
    return pointField::null();
}


const vectorField& globalProcessorPointPatch::pointNormals() const
{
    notImplemented("globalProcessorPointPatch::pointNormals() const");
    return vectorField::null();
}


triFaceList globalProcessorPointPatch::faceTriangles
(
    const label
) const
{
    notImplemented
    (
        "processorPointPatch::faceTriangles(label faceID) const"
    );

    return List<triFace>::null();
}


const edgeList& globalProcessorPointPatch::meshEdges() const
{
    notImplemented("globalProcessorPointPatch::meshEdges() const");
    return edgeList::null();
}


const labelList& globalProcessorPointPatch::sharedEdgeAddr() const
{
    notImplemented("globalProcessorPointPatch::sharedEdgeAddr() const");
    return labelList::null();
}


const edgeList& globalProcessorPointPatch::meshCutEdges() const
{
    notImplemented("globalProcessorPointPatch::meshCutEdges() const");
    return edgeList::null();
}


const scalarField& globalProcessorPointPatch::meshCutEdgeMask() const
{
    notImplemented("globalProcessorPointPatch::meshCutEdgeMask() const");
    return scalarField::null();
}


const labelList& globalProcessorPointPatch::localEdgeIndices() const
{
    notImplemented("globalProcessorPointPatch::localEdgeIndices() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::cutEdgeIndices() const
{
    notImplemented("globalProcessorPointPatch::cutEdgeIndices() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::cutEdgeOwnerIndices() const
{
    notImplemented("globalProcessorPointPatch::cutEdgeOwnerIndices() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::cutEdgeOwnerStart() const
{
    notImplemented("globalProcessorPointPatch::cutEdgeOwnerStart() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::cutEdgeNeighbourIndices() const
{
    notImplemented
    (
        "globalProcessorPointPatch::cutEdgeNeighbourIndices() const"
    );
    return labelList::null();
}


const labelList& globalProcessorPointPatch::cutEdgeNeighbourStart() const
{
    notImplemented("globalProcessorPointPatch::cutEdgeNeighbourStart() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::doubleCutEdgeIndices() const
{
    notImplemented("globalProcessorPointPatch::doubleCutEdgeIndices() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::doubleCutOwner() const
{
    notImplemented("globalProcessorPointPatch::doubleCutOwner() const");
    return labelList::null();
}


const labelList& globalProcessorPointPatch::doubleCutNeighbour() const
{
    notImplemented("globalProcessorPointPatch::doubleCutNeighbour() const");
    return labelList::null();
}


const scalarField& globalProcessorPointPatch::ownNeiDoubleMask() const
{
    notImplemented("globalProcessorPointPatch::ownNeiDoubleMask() const");
    return scalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
