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
    Pstream

Description

\*---------------------------------------------------------------------------*/

#ifndef PstreamReduceOpsMPINew_H
#define PstreamReduceOpsMPINew_H

#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void reduce(scalar& Value, sumOp<scalar> bop)
{
    if (Pstream::parRun())
    {
        MPI_Status status;

        int myProcNo = Pstream::myProcNo();
        int nProcs = Pstream::nProcs();

        //
        // receive from children
        //
        int level = 1;
        int thisLevelOffset = 2;
        int childLevelOffset = thisLevelOffset/2;
        int childProcId = 0;

        while
        (
            (childLevelOffset < nProcs)
         && (myProcNo % thisLevelOffset) == 0
        )
        {
            childProcId = myProcNo + childLevelOffset;

            scalar value;

            if (childProcId < nProcs)
            {
                if
                (
                    MPI_Recv
                    (
                        &value,
                        1,
                        MPI_DOUBLE,
                        Pstream::procID(childProcId),
                        Pstream::msgType(),
                        MPI_COMM_WORLD,
                        &status
                    )
                )
                {
                    FatalErrorIn("reduce(scalar& Value, sumOp<scalar> sumOp)")
                        << "MPI_Recv failed"
                        << abort(FatalError);
                }

	            Value = bop(Value, value);
	        }

            level++;
            thisLevelOffset <<= 1;
            childLevelOffset = thisLevelOffset/2;
        }

        //
        // send and receive from parent
        //
        if (!Pstream::master())
        {
            int parentId = myProcNo - (myProcNo % thisLevelOffset);

            if
            (
                MPI_Send
                (
                    &Value,
                    1,
                    MPI_DOUBLE,
                    Pstream::procID(parentId),
                    Pstream::msgType(),
                    MPI_COMM_WORLD
                )
            )
            {
                FatalErrorIn("reduce(scalar& Value, sumOp<scalar> sumOp)")
                    << "MPI_Send failed"
                    << abort(FatalError);
            }

            if
            (
                MPI_Recv
                (
                    &Value,
                    1,
                    MPI_DOUBLE,
                    Pstream::procID(parentId),
                    Pstream::msgType(),
                    MPI_COMM_WORLD,
                    &status
                )
            )
            {
                FatalErrorIn("reduce(scalar& Value, sumOp<scalar> sumOp)")
                    << "MPI_Recv failed"
                    << abort(FatalError);
            }
        }


        //
        // distribute to my children
        //
        level--;
        thisLevelOffset >>= 1;
        childLevelOffset = thisLevelOffset/2;

        while (level > 0)
        {
            childProcId = myProcNo + childLevelOffset;

            if (childProcId < nProcs)
            {
                if
                (
                    MPI_Send
                    (
                        &Value,
                        1,
                        MPI_DOUBLE,
                        Pstream::procID(childProcId),
                        Pstream::msgType(),
                        MPI_COMM_WORLD
                    )
                )
                {
                    FatalErrorIn("reduce(scalar& Value, sumOp<scalar> sumOp)")
                        << "MPI_Send failed"
                        << abort(FatalError);
                }
            }

            level--;
            thisLevelOffset >>= 1;
            childLevelOffset = thisLevelOffset/2;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
