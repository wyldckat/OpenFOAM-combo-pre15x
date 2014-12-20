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
    Read token and binary block from IPstream using pvm.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "IPstream.H"

#include <pvm3.h>

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

    int bufid, len, tag, tid;

    // If the buffer size is not specified then probe the incomming message

    if (!bufSize)
    {
        // Probe read buffer until message arrives.
        while (!(bufid = pvm_probe(procID(fromProcNo_), msgType())));

        // When the message arrives find it's size
        pvm_bufinfo(bufid, &len, &tag, &tid);

        // Resize buffer to message size
        buf_.setSize(len);
    }


    // Read message into buffer

    if
    (
        pvm_precv
        (
            procID(fromProcNo_),
            msgType(),
            buf_.begin(),
            buf_.size(),
            PVM_BYTE,
            &tid, &tag, &len
        ) != PvmOk
    )
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "pvm_precv cannot receive incomming message"
            << ::abort;
    }


    // Check size of message read

    if (len > buf_.size())
    {
        FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
            << "buffer (" << buf_.size()
            << ") not large enough for incomming message (" << len << ')'
            << ::abort;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
