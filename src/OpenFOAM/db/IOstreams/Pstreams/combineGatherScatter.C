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
    Variant of gather, scatter.
    Normal gather uses:
    - construct null and read (>>) from Istream
    - binary operator and assignment operator to combine values

    combineGather uses:
    - construct from Istream
    - modify operator which modifies its lhs

\*---------------------------------------------------------------------------*/

#include "OPstream.H"
#include "IPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class T, class CombineOp>
void Pstream::combineGather
(
    const List<Pstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop
)
{
    if (Pstream::parRun())
    {
        // Get my communication order
        const commsStruct& myComm = comms[Pstream::myProcNo()];

        // Receive from my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            //IPstream fromBelow(belowID, sizeof(T));
            IPstream fromBelow(belowID);
            T value(fromBelow);

            if (debug)
            {
                Pout<< " received from "
                    << belowID << " data:" << value << endl;
            }
            cop(Value, value);
        }

        // Send up Value
        if (myComm.above() != -1)
        {
            //OPstream toAbove(myComm.above(), sizeof(T), false);
            OPstream toAbove(myComm.above(), 0, false);

            if (debug)
            {
                Pout<< " sending to " << myComm.above()
                    << " data:" << Value << endl;
            }
            toAbove << Value;
        }
    }
}


template <class T, class CombineOp>
void Pstream::combineGather(T& Value, const CombineOp& cop)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        combineGather(Pstream::linearCommunication(), Value, cop);
    }
    else
    {
        combineGather(Pstream::treeCommunication(), Value, cop);
    }
}


template <class T>
void Pstream::combineScatter(const List<Pstream::commsStruct>& comms, T& Value)
{
    if (Pstream::parRun())
    {
        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            //IPstream fromAbove(myComm.above(), sizeof(T));
            IPstream fromAbove(myComm.above());
            Value = T(fromAbove);

            if (debug)
            {
                Pout<< " received from "
                    << myComm.above() << " data:" << Value << endl;
            }
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            label belowID = myComm.below()[belowI];

            //OPstream toBelow(belowID, sizeof(T), false);
            OPstream toBelow(belowID, 0, false);

            toBelow << Value;

            if (debug)
            {
                Pout<< " sent to " << belowID << " data:" << Value << endl;
            }
        }
    }
}


template <class T>
void Pstream::combineScatter(T& Value)
{
    if (Pstream::nProcs() < Pstream::nProcsSimpleSum)
    {
        combineScatter(Pstream::linearCommunication(), Value);
    }
    else
    {
        combineScatter(Pstream::treeCommunication(), Value);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
