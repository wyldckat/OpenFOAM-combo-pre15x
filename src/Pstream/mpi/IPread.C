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

    if (!read(fromProcNo_, buf_.begin(), buf_.size()))
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "read failed";
        abort();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool IPstream::read
(
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize
)
{
    MPI_Status status;
    int messageSize;

    if
    (
        MPI_Recv
        (
            buf,
            bufSize,
            MPI_PACKED,
            procID(fromProcNo),
            msgType(),
            MPI_COMM_WORLD,
            &status
        )
    )
    {
        FatalErrorIn
        (
            "IPstream::read"
            "(const int fromProcNo, char* buf, std::streamsize bufSize)"
        )   << "MPI_Recv cannot receive incomming message" << endl;

        return false;
    }


    // Check size of message read

    MPI_Get_count(&status, MPI_BYTE, &messageSize);

    if (messageSize > bufSize)
    {
        FatalErrorIn
        (
            "IPstream::read"
            "(const int fromProcNo, char* buf, std::streamsize bufSize)"
        )   << "buffer (" << int(bufSize)
            << ") not large enough for incomming message ("
            << messageSize << ')' << endl;

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
