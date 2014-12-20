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

#include "processorTopology.H"
#include "processorPolyPatch.H"
#include "ListOps.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::processorTopology::procNeighbours
(
    const polyPatchList& patches
)
{
    // Determine number of processor neighbours and max neighbour id.

    label nNeighbours = 0;

    label maxNb = labelMin;

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (typeid(patch) == typeid(processorPolyPatch))
        {
            const processorPolyPatch& procPatch = 
                refCast<const processorPolyPatch>(patch);

            nNeighbours++;

            maxNb = max(maxNb, procPatch.neighbProcNo());
        }
    }

    labelList neighbours(nNeighbours);

    procPatchMap_.setSize(maxNb + 1);
    procPatchMap_ = -1;

    nNeighbours = 0;

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (typeid(patch) == typeid(processorPolyPatch))
        {
            const processorPolyPatch& procPatch = 
                refCast<const processorPolyPatch>(patch);

            neighbours[nNeighbours++] = procPatch.neighbProcNo();

            // Construct reverse map
            procPatchMap_[procPatch.neighbProcNo()] = patchI;
        }
    }

    return neighbours;
}


void Foam::processorTopology::calcAddressing() const
{
    // Construct cellFaces

    cellFacesPtr_ = new labelListList(size());

    labelListList& cellFaces = *cellFacesPtr_;


    label faceI = 0;

    forAll(*this, procI)
    {
        const labelList& myNbs = operator[](procI);

        cellFaces[procI].setSize(myNbs.size());

        cellFaces[procI] = -1;

        // Introduce face to lower numbered cell (so cellFaces already sized)
        forAll(myNbs, nbI)
        {
            label nb = myNbs[nbI];

            if (nb < procI)
            {
                // Insert into first free spot in cellFaces[procI]
                labelList& procFaces = cellFaces[procI];

                label pos = findIndex(procFaces, -1);

                procFaces[pos] = faceI;

                // Insert into first free spot in cellFaces[nb]
                labelList& nbFaces = cellFaces[nb];

                pos = findIndex(nbFaces, -1);

                nbFaces[pos] = faceI;

                faceI++;
            }
        }

    }

    // Do faceCells

    faceCellsPtr_ = new labelListList(faceI);

    labelListList& faceCells = *faceCellsPtr_;


    forAll(faceCells, faceI)
    {
        faceCells[faceI].setSize(2);
    }

    forAll(*this, procI)
    {
        const labelList& myNbs = operator[](procI);

        const labelList& myFaces = cellFaces[procI];
    
        forAll(myNbs, nbI)
        {
            label nb = myNbs[nbI];

            label faceI = myFaces[nbI];

            faceCells[faceI][0] = procI;
            faceCells[faceI][1] = nb;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::processorTopology::processorTopology(const polyPatchList& patches)
:
    labelListList(Pstream::nProcs()),
    procPatchMap_(),
    cellFacesPtr_(NULL),
    faceCellsPtr_(NULL)
{
    // Fill my 'slot' with my neighbours
    operator[](Pstream::myProcNo()) = procNeighbours(patches);

    // Distribute to all processors
    Pstream::gatherList(*this);
    Pstream::scatterList(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorTopology::~processorTopology()
{
    deleteDemandDrivenData(cellFacesPtr_);
    deleteDemandDrivenData(faceCellsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::processorTopology::getFace
(
    const label myProcID,
    const label nbrProcID
) const
{
    const labelList& myFaces = cellFaces()[myProcID];

    forAll(myFaces, myFaceI)
    {
        label faceI = myFaces[myFaceI];

        const labelList& fCells = faceCells()[faceI];

        if ((fCells[0] == nbrProcID) || (fCells[1] == nbrProcID))
        {
            return faceI;
        }
    }

    FatalErrorIn("processorTopology::getFace(const label, const label)")
        << "Cannot find connection (face) between processor " << myProcID
        << " and " << nbrProcID << abort(FatalError);

    return -1;
}


Foam::label Foam::processorTopology::getFace(const label nbrProcID) const
{
    return getFace(Pstream::myProcNo(), nbrProcID);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
