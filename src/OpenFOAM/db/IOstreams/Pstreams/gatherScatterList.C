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
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    Values[Pstream::myProcNo()].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T>
void Pstream::gatherList
(
    const List<Pstream::commsStruct>& comms,
    List<T>& Values
)
{
    if (Pstream::parRun())
    {
        if (Values.size() != Pstream::nProcs())
        {
            FatalErrorIn
            (
                "Pstream::gatherList(const List<Pstream::commsStruct>&"
                ", List<T>)"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << Pstream::nProcs()
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            //IPstream fromBelow(belowID, Pstream::nProcs()*sizeof(T));
            IPstream fromBelow(belowID);

            // Receive from belowID first and put in correct slot in Values
            fromBelow >> Values[belowID];

            if (debug)
            {
                Pout<< " received through "
                    << belowID << " data from:" << belowID
                    << " data:" << Values[belowID] << endl;
            }

            // Receive from all other processors below belowID
            const labelList& belowLeaves = comms[belowID].allBelow();

            forAll(belowLeaves, leafI)
            {
                label leafID = belowLeaves[leafI];

                fromBelow >> Values[leafID];

                if (debug)
                {
                    Pout<< " received through "
                        << belowID << " data from:" << leafID
                        << " data:" << Values[leafID] << endl;
                }
            }
        }

        // Send up from Values:
        // - my own value first
        // - all belowLeaves next
        if (myComm.above() != -1)
        {
            //OPstream toAbove(myComm.above(), Pstream::nProcs()*sizeof(T), false);
            OPstream toAbove(myComm.above(), 0, false);

            if (debug)
            {
                Pout<< " sending to " << myComm.above()
                    << " data from me:" << Pstream::myProcNo()
                    << " data:" << Values[Pstream::myProcNo()] << endl;
            }
            toAbove << Values[Pstream::myProcNo()];

            forAll(myComm.allBelow(), leafI)
            {
                label leafID = myComm.allBelow()[leafI];

                if (debug)
                {
                    Pout<< " sending to "
                        << myComm.above() << " data from:" << leafID
                        << " data:" << Values[leafID] << endl;
                }
                toAbove << Values[leafID];
            }
        }
    }
}


template <class T>
void Pstream::gatherList(List<T>& Values)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        gatherList(Pstream::linearCommunication(), Values);
    }
    else
    {
        gatherList(Pstream::treeCommunication(), Values);
    }
}


template <class T>
void Pstream::scatterList
(
    const List<Pstream::commsStruct>& comms,
    List<T>& Values
)
{
    if (Pstream::parRun())
    {
        if (Values.size() != Pstream::nProcs())
        {
            FatalErrorIn
            (
                "Pstream::scatterList(const List<Pstream::commsStruct>&"
                ", List<T>)"
            )   << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << Pstream::nProcs()
                << Foam::abort(FatalError);
        }

        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            //IPstream fromAbove(myComm.above(), Pstream::nProcs()*sizeof(T));
            IPstream fromAbove(myComm.above());

            forAll(myComm.allNotBelow(), leafI)
            {
                label leafID = myComm.allNotBelow()[leafI];

                fromAbove >> Values[leafID];

                if (debug)
                {
                    Pout<< " received through "
                        << myComm.above() << " data for:" << leafID
                        << " data:" << Values[leafID] << endl;
                }
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            //OPstream toBelow(belowID, Pstream::nProcs()*sizeof(T), false);
            OPstream toBelow(belowID, 0, false);

            // Send data destined for all other processors below belowID
            const labelList& belowLeaves = comms[belowID].allNotBelow();

            forAll(belowLeaves, leafI)
            {
                label leafID = belowLeaves[leafI];

                toBelow << Values[leafID];

                if (debug)
                {
                    Pout<< " sent through "
                        << belowID << " data for:" << leafID
                        << " data:" << Values[leafID] << endl;
                }
            }
        }
    }
}


template <class T>
void Pstream::scatterList(List<T>& Values)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        scatterList(Pstream::linearCommunication(), Values);
    }
    else
    {
        scatterList(Pstream::treeCommunication(), Values);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
