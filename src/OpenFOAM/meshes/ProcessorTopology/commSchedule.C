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

#include "commSchedule.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::commSchedule::nBusyNbs
(
    const labelListList& cellFaces,
    const labelListList& faceCells,
    const boolList& busy,
    const label cellI
) const
{
    label nBusy = 0;

    const labelList& myFaces = cellFaces[cellI];

    forAll(myFaces, myFaceI)
    {
        label faceI = myFaces[myFaceI];

        const labelList& myCells = faceCells[faceI];

        // Check if the other cell on the face is busy.
        if
        (
            (
                (myCells[0] == cellI)
             && busy[myCells[1]]
            )
         || (
                (myCells[1] == cellI)
             && busy[myCells[0]]
            )
        )
        {
            nBusy++;
        }
    }
    return nBusy;
}



// Find faces which are not yet used in a schedule and
// assign them to the current schedule such that a cell is only used once.
bool Foam::commSchedule::scheduleIteration
(
    const labelListList& cellFaces,
    const labelListList& faceCells,
    const label commIter
)
{
    // Whether cell is busy (i.e. is involved in a communication with
    // any of its neighbours)
    boolList busy(cellFaces.size(), false);

    // To detect whether all faces have been scheduled
    bool seenUnscheduledFace = false;

    do
    {
        // Get an unscheduled face. Alg is:
        // - for all candidate faces
        // - choose one with cells with max number of busy neighbours.
        //   (should cluster next to previously scheduled cellpair)
        // - if multiple ones with same number choose out of these one
        //   with min number of faces
        //   (if no scheduled cellpairs chooses corner cell)
        label routeFaceI = -1;

        label maxBusy = labelMin;

        label minNFaces = labelMax;

        forAll(faceCells, faceI)
        {
            if (faceSchedule_[faceI] != -1)
            {
                // Already scheduled. Skip.
                continue;
            }

            // Remember whether we have done something in this loop.
            seenUnscheduledFace = true;

            const labelList& myCells = faceCells[faceI];

            label cell0 = myCells[0];
            label cell1 = myCells[1];

            if (!busy[cell0] && !busy[cell1])
            {
                // Face is ok candidate: unscheduled and cells not busy

                label busy0 = nBusyNbs(cellFaces, faceCells, busy, cell0);
                label busy1 = nBusyNbs(cellFaces, faceCells, busy, cell1);

                //label totBusy = max(busy0, busy1);
                label totBusy = busy0 + busy1;

                if (totBusy > maxBusy)
                {
                    maxBusy = totBusy;

                    minNFaces = labelMax;

                    routeFaceI = faceI;
                }
                else if (totBusy == maxBusy)
                {
                    label nFaces =
                        min
                        (
                            cellFaces[cell0].size(),
                            cellFaces[cell1].size()
                        );

                    if (nFaces < minNFaces)
                    {
                        minNFaces = nFaces;

                        routeFaceI = faceI;
                    }
                }
            }
        }


        if (routeFaceI == -1)
        {
            break;
        }

        // Make connection across routeFaceI
        const labelList& myCells = faceCells[routeFaceI];

        faceSchedule_[routeFaceI] = commIter;

        busy[myCells[0]] = true;

        busy[myCells[1]] = true;
    }
    while (true);

    if (seenUnscheduledFace)
    {
        // At least one unscheduled face has been encountered.(might be
        // scheduled now)
        return true;
    }
    else
    {
        return false;
    }
}


// Schedule all. Fill faceSchedule_ and *this
void Foam::commSchedule::scheduleAll
(
    const labelListList& cellFaces,
    const labelListList& faceCells
)
{
    // Current communication iteration.
    label commIter = 0;

    do
    {
        // Get all proc-proc communication across unscheduled face.

        if (!scheduleIteration(cellFaces, faceCells, commIter))
        {
            // Done
            break;
        }

        commIter++;

    } while (true);


    // Rework face wise schedule into cell wise. Gives processor to communicate
    // with. (i.e. cellCells with cells ordered according to communication)
    // Does not determine yet whether to send or receive.

    forAll(cellFaces, cellI)
    {
        const labelList& myFaces = cellFaces[cellI];

        // Sort faces according to schedule.
        SortableList<label> sortedSchedule(myFaces.size());
        forAll(myFaces, myFaceI)
        {
            sortedSchedule[myFaceI] = faceSchedule_[myFaces[myFaceI]];
        }
        sortedSchedule.sort();


        labelList& cellSchedule = operator[](cellI);

        cellSchedule.setSize(myFaces.size());

        const labelList& indcs = sortedSchedule.indices();

        forAll(indcs, indexI)
        {
            label faceI = myFaces[indcs[indexI]];

            const labelList& myNbs = faceCells[faceI];

            if (myNbs[0] == cellI)
            {
                cellSchedule[indexI] = myNbs[1];
            }
            else
            {
                cellSchedule[indexI] = myNbs[0];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from separate addressing
Foam::commSchedule::commSchedule
(
    const labelListList& cellFaces,
    const labelListList& faceCells
)
:
    labelListList(cellFaces.size()),
    faceSchedule_(faceCells.size(), -1)
{
    scheduleAll(cellFaces, faceCells);
}


// ************************************************************************* //
