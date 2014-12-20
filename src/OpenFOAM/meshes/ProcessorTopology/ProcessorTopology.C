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

#include "ProcessorTopology.H"
#include "ListSearch.H"
#include "commSchedule.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Patch, class ProcPatch>
Foam::labelList Foam::ProcessorTopology<Patch, ProcPatch>::procNeighbours
(
    const ptrList<Patch>& patches
)
{
    // Determine number of processor neighbours and max neighbour id.

    label nNeighbours = 0;

    label maxNb = labelMin;

    forAll(patches, patchi)
    {
        const Patch& patch = patches[patchi];

        if (typeid(patch) == typeid(ProcPatch))
        {
            const ProcPatch& procPatch = 
                refCast<const ProcPatch>(patch);

            nNeighbours++;

            maxNb = max(maxNb, procPatch.neighbProcNo());
        }
    }

    labelList neighbours(nNeighbours);

    procPatchMap_.setSize(maxNb + 1);
    procPatchMap_ = -1;

    nNeighbours = 0;

    forAll(patches, patchi)
    {
        const Patch& patch = patches[patchi];

        if (typeid(patch) == typeid(ProcPatch))
        {
            const ProcPatch& procPatch = 
                refCast<const ProcPatch>(patch);

            neighbours[nNeighbours++] = procPatch.neighbProcNo();

            // Construct reverse map
            procPatchMap_[procPatch.neighbProcNo()] = patchi;
        }
    }

    return neighbours;
}


template<class Patch, class ProcPatch>
void Foam::ProcessorTopology<Patch, ProcPatch>::calcAddressing()
{
    cellFaces_.setSize(size());

    // Size and set to -1.
    forAll(*this, procI)
    {
        cellFaces_[procI].setSize(operator[](procI).size());
        cellFaces_[procI] = -1;
    }

    label nFaces = 0;

    forAll(*this, procI)
    {
        // Sort my list of neighbours
        labelList myNbs(operator[](procI));
        sort(myNbs);

        // Introduce face to lower numbered cell. (note: since myNbs already
        // sorted this is equivalent to upper-triangular order)
        forAll(myNbs, nbI)
        {
            label nb = myNbs[nbI];

            if (nb < procI)
            {
                // Insert into first free spot in cellFaces_[procI]
                labelList& procFaces = cellFaces_[procI];

                label pos = findIndex(procFaces, -1);

                if (pos == -1)
                {
                    FatalErrorIn
                    (
                        "ProcessorTopology<Patch, ProcPatch>::"
                        "calcAddressing()"
                    )   << "Cannot find empty slot (-1) in faces "
                        << procFaces << " of processor " << procI << endl
                        << "when trying to insert new face for connection to"
                        << " processor " << nb << endl
                        << abort(FatalError);
                }

                procFaces[pos] = nFaces;

                // Insert into first free spot in cellFaces_[nb]
                labelList& nbFaces = cellFaces_[nb];

                pos = findIndex(nbFaces, -1);

                if (pos == -1)
                {
                    FatalErrorIn
                    (
                        "ProcessorTopology<Patch, ProcPatch>::"
                        "calcAddressing()"
                    )   << "Cannot find empty slot (-1) in faces "
                        << nbFaces << " of processor " << nb << endl
                        << "when trying to insert new face for connection to"
                        << " processor " << procI << endl
                        << abort(FatalError);
                }

                nbFaces[pos] = nFaces;

                nFaces++;
            }
        }

    }

    // Do faceCells

    faceCells_.setSize(nFaces);

    // Size and set to -1 (needed to recognize unfilled slots)
    forAll(faceCells_, faceI)
    {
        faceCells_[faceI].setSize(2);
        faceCells_[faceI] = -1;
    }

    forAll(*this, procI)
    {
        const labelList& myFaces = cellFaces_[procI];
    
        forAll(myFaces, i)
        {
            label faceI = myFaces[i];

            if (faceCells_[faceI][0] == -1)
            {
                faceCells_[faceI][0] = procI;
            }
            else if (faceCells_[faceI][1] == -1)
            {
                faceCells_[faceI][1] = procI;
            }
            else
            {
                FatalErrorIn
                (
                    "ProcessorTopology<Patch, ProcPatch>::calcAddressing()"
                )   << "More than two processors using face " << faceI << endl
                    << "Processor1:" << faceCells_[faceI][0] << endl
                    << "Processor3:" << faceCells_[faceI][1] << endl
                    << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Patch, class ProcPatch>
Foam::ProcessorTopology<Patch, ProcPatch>::ProcessorTopology
(
    const ptrList<Patch>& patches
)
:
    labelListList(Pstream::nProcs()),
    patchSchedule_(2*patches.size())
{
    if (Pstream::parRun())
    {
        // Fill my 'slot' with my neighbours
        operator[](Pstream::myProcNo()) = procNeighbours(patches);

        // Distribute to all processors
        Pstream::gatherList(*this);
        Pstream::scatterList(*this);

        calcAddressing();
    }

    if 
    (
        Pstream::parRun()
     && debug::optimisationSwitch("scheduledTransfer", false)
    )
    {
        label patchEvali = 0;

        forAll(patches, patchi)
        {
            if (typeid(patches[patchi]) != typeid(ProcPatch))
            {
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali].init = true;
                patchSchedule_[patchEvali++].bufferedTransfer = false;
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali].init = false;
                patchSchedule_[patchEvali++].bufferedTransfer = false;
            }
        }

        labelList schedule
        (
            commSchedule
            (
                this->cellFaces(),
                this->faceCells()
            )[Pstream::myProcNo()]
        );

        forAll(schedule, iter)
        {
            label nb = schedule[iter];
            label patchi = procPatchMap_[nb];

            if (Pstream::myProcNo() > nb)
            {
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali].init = true;
                patchSchedule_[patchEvali++].bufferedTransfer = false;
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali].init = false;
                patchSchedule_[patchEvali++].bufferedTransfer = false;
            }
            else
            {
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali].init = false;
                patchSchedule_[patchEvali++].bufferedTransfer = false;
                patchSchedule_[patchEvali].patch = patchi;
                patchSchedule_[patchEvali].init = true;
                patchSchedule_[patchEvali++].bufferedTransfer = false;
            }
        }
    }
    else
    {
        label patchEvali = 0;

        forAll(patches, patchi)
        {
            patchSchedule_[patchEvali].patch = patchi;
            patchSchedule_[patchEvali].init = true;
            patchSchedule_[patchEvali++].bufferedTransfer = true;
        }

        forAll(patches, patchi)
        {
            patchSchedule_[patchEvali].patch = patchi;
            patchSchedule_[patchEvali].init = false;
            patchSchedule_[patchEvali++].bufferedTransfer = true;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Patch, class ProcPatch>
Foam::label Foam::ProcessorTopology<Patch, ProcPatch>::getFace
(
    const label myProcID,
    const label nbrProcID
) const
{
    const labelList& myFaces = cellFaces_[myProcID];

    forAll(myFaces, myFaceI)
    {
        label faceI = myFaces[myFaceI];

        const labelList& fCells = faceCells_[faceI];

        if ((fCells[0] == nbrProcID) || (fCells[1] == nbrProcID))
        {
            return faceI;
        }
    }

    FatalErrorIn
    (
        "ProcessorTopology<Patch, ProcPatch>::getFace(const label, "
        "const label)"
    )   << "Cannot find connection (face) between processor " << myProcID
        << " and " << nbrProcID
        << abort(FatalError);

    return -1;
}


template<class Patch, class ProcPatch>
Foam::label
Foam::ProcessorTopology<Patch, ProcPatch>::getFace
(
    const label nbrProcID
) const
{
    return getFace(Pstream::myProcNo(), nbrProcID);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
