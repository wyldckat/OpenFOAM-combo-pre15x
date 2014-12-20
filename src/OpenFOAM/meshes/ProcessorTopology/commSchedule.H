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

Class
    commSchedule

Description
    Given the processor topology works out a non-blocking communication
    schedule. Instead of doing all sends first and having a lot of buffering
    this works out how to have each processor do a sequence of sends and
    receives such that there can never be the case that two processors are
    waiting for each other to send something.

    *this[procI] is for processor procI a list with for every iteration
    the neighbour it has to communicate with. It does not prespecify which
    half should start sending or receiving. You might want to group the two
    close together or separate them in time as long as one is sending and the
    other is receiving.


    A typical use with matched send and receives is something like:

    const labelList& mySchedule = operator[](Pstream::myProcNo()];

    forAll(mySchedule, iterI)
    {
        label nb = mySchedule[iterI];

        if (nb < Pstream::myProcNo())
        {
            // Send to neighbour
            ..
            // Receive from neighbour
            ..
        }
        else
        {
            // Receive from neighbour
            ..
            // Send to neighbour
            ..
        }
    }
    


SourceFiles
    commSchedule.C

\*---------------------------------------------------------------------------*/

#ifndef commSchedule_H
#define commSchedule_H

#include "labelList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class commSchedule Declaration
\*---------------------------------------------------------------------------*/

class commSchedule
:
    public labelListList
{
    // Private data

        //- Schedule per face. Gives for each face whether communication
        //  takes place. Is -1 or the iteration in which communication should
        //  take place. *this is the same information but sorted per iteration.
        labelList faceSchedule_;


    // Private Member Functions

        //- Count number of 'busy' neighbouring processors
        label nBusyNbs
        (
            const labelListList& cellFaces,
            const labelListList& faceCells,
            const boolList& busy,
            const label cellI
        ) const;

        //- Find faces which are not yet used in a schedule and
        //  assign them to the current schedule such that a cell is only
        //  used once. Returns true if any faces have been 'scheduled', false
        //  otherwise
        bool scheduleIteration
        (
            const labelListList& cellFaces,
            const labelListList& faceCells,
            const label commIter
        );

        //- Schedule all. Fill faceSchedule_ and *this
        void scheduleAll
        (
            const labelListList& cellFaces,
            const labelListList& faceCells
        );


public:

    // Constructors

        //- Construct from addressing
        commSchedule
        (
            const labelListList& cellFaces,
            const labelListList& faceCells
        );


    // Member Functions

        //- Schedule on a face by face basis.
        const labelList& faceSchedule() const
        {
            return faceSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
