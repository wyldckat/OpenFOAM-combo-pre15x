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
    Tetrahedral point field decomposer.

\*---------------------------------------------------------------------------*/

#include "geometricTetPointFieldDecomposer.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate point addressing
void geometricTetPointFieldDecomposer::calcAddressing() const
{
    if (directAddressingPtr_)
    {
        FatalErrorIn
        (
            "void geometricTetPointFieldDecomposer::calcAddressing() const"
        )   << "addressing already calculated"
            << abort(FatalError);
    }

    // Allocate the addressing
    directAddressingPtr_ = new labelList(processorMesh_.nPoints(), -1);
    labelList& addr = *directAddressingPtr_;

    label nAddr = 0;

    // Insert point addressing
    forAll (pointAddressing_, pointI)
    {
        addr[nAddr] = pointAddressing_[pointI];
        nAddr++;
    }

#   ifdef FACE_DECOMP
    // Insert face addressing.  Only for face decomposition
    const label faceOffset = originalMesh_.faceOffset();

    forAll (faceAddressing_, faceI)
    {
        // Remember to decrement the index by one (turning index)
        addr[nAddr] = faceOffset + mag(faceAddressing_[faceI]) - 1;
        nAddr++;
    }
#   endif

    // Insert cell addressing
    const label cellOffset = originalMesh_.cellOffset();

    forAll (cellAddressing_, cellI)
    {
        addr[nAddr] = cellOffset + cellAddressing_[cellI];
        nAddr++;
    }
}


const labelList& geometricTetPointFieldDecomposer::directAddressing() const
{
    if (!directAddressingPtr_)
    {
        calcAddressing();
    }

    return *directAddressingPtr_;
}


// Calculate patch addressing
void geometricTetPointFieldDecomposer::
tetPolyPatchFieldDecomposer::calcPatchAddressing() const
{
    if (directPatchAddressingPtr_)
    {
        FatalErrorIn
        (
            "void geometricTetPointFieldDecomposer::"
            "tetPolyPatchFieldDecomposer::calcPatchAddressing() const"
        )   << "addressing already calculated"
            << abort(FatalError);
    }

    // Allocate the addressing
    directPatchAddressingPtr_ = new labelList(targetPatch().size(), -1);
    labelList& addr = *directPatchAddressingPtr_;

    // Algorithm:
    // Go to the source patch, create a lookup list the size of all
    // points in the mesh and then gather the points for the current
    // patch.
    labelList pointLookup
        (sourcePatch().boundaryMesh().mesh().nPoints(), -1);

    const labelList& sourcePatchPoints = sourcePatch().meshPoints();

    forAll (sourcePatchPoints, pointI)
    {
        pointLookup[sourcePatchPoints[pointI]] = pointI;
    }

    // Gather the information
    
    const labelList& targetPatchPoints = targetPatch().meshPoints();

    forAll (targetPatchPoints, pointI)
    {
        addr[pointI] =
            pointLookup[directAddressing_[targetPatchPoints[pointI]]];
    }

    if (addr.size() && min(addr) < 0)
    {
        FatalErrorIn
        (
            "void geometricTetPointFieldDecomposer::"
            "tetPolyPatchFieldDecomposer::calcPatchAddressing() const"
        )   << "error in addressing"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
geometricTetPointFieldDecomposer::geometricTetPointFieldDecomposer
(
    const tetPolyMesh& originalMesh,
    const tetPolyMesh& processorMesh,
    const labelList& pointAddressing,
    const labelList& faceAddressing,
    const labelList& cellAddressing,
    const labelList& boundaryAddressing
)
:
    originalMesh_(originalMesh),
    processorMesh_(processorMesh),
    pointAddressing_(pointAddressing),
    faceAddressing_(faceAddressing),
    cellAddressing_(cellAddressing),
    boundaryAddressing_(boundaryAddressing),
    patchFieldDecompPtrs_
    (
        processorMesh_.boundary().size(),
        reinterpret_cast<tetPolyPatchFieldDecomposer*>(NULL)
    ),
    directAddressingPtr_(NULL)
{
    // Hook the patch field decomposers for all non-processor patch fields
    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFieldDecompPtrs_[patchI] =
                new tetPolyPatchFieldDecomposer
                (
                    originalMesh_.boundary()[boundaryAddressing_[patchI]],
                    processorMesh_.boundary()[patchI],
                    directAddressing()
                );
        }
        // for processor patches the pointer stays null
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

geometricTetPointFieldDecomposer::~geometricTetPointFieldDecomposer()
{
    forAll (patchFieldDecompPtrs_, patchI)
    {
        if (patchFieldDecompPtrs_[patchI] != NULL)
        {
            delete(patchFieldDecompPtrs_[patchI]);
            patchFieldDecompPtrs_[patchI] = NULL;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
