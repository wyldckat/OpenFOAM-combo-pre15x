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
    Read token and binary block from IPstream

\*---------------------------------------------------------------------------*/

#include "IPstream.H"

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

IPstream::IPstream
(
    const int fromProcNo,
    const label bufSize,
    streamFormat format,
    versionNumber version
)
:
    Pstream(bufSize),
    Istream(format, version),
    fromProcNo_(fromProcNo)
{
    setOpened();
    setGood();

    MPI_Status status;
    int messageSize;

    // If the buffer size is not specified then probe the incomming message

    if (!bufSize)
    {
        MPI_Probe(procID(fromProcNo_), msgType(), MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_BYTE, &messageSize);

        buf_.setSize(messageSize);
    }

    // Read message into buffer

    if
    (
        MPI_Recv
        (
            buf_.begin(),
            buf_.size(),
            MPI_PACKED,
            procID(fromProcNo_),
            msgType(),
            MPI_COMM_WORLD,
            &status
        )
    )
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "MPI_Recv cannot receive incomming message";
        abort();
    }


    // Check size of message read

    MPI_Get_count(&status, MPI_BYTE, &messageSize);

    if (messageSize > buf_.size())
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "buffer (" << buf_.size()
            << ") not large enough for incomming message ("
            << messageSize<< ')';
        abort();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
