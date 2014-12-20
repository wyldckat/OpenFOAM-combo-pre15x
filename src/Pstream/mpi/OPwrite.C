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
    Write primitive and binary block from OPstream

\*---------------------------------------------------------------------------*/

#include "OPstream.H"

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

OPstream::~OPstream()
{
    bool transferFailed = true;

    if (bufferedTransfer_)
    {
        transferFailed = MPI_Bsend
        (
            buf_.begin(),
            bufPosition_,
            MPI_PACKED,
            procID(toProcNo_),
            msgType(),
            MPI_COMM_WORLD
        );
    }
    else
    {
        transferFailed = MPI_Send
        (
            buf_.begin(),
            bufPosition_,
            MPI_PACKED,
            procID(toProcNo_),
            msgType(),
            MPI_COMM_WORLD
        );
    }

    if (transferFailed)
    {
        FatalErrorIn("OPstream::~OPstream()")
            << "MPI_Bsend cannot send outgoing message";
        abort();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
