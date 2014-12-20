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

\*---------------------------------------------------------------------------*/

#include "processorTopology.H"
#include "commSchedule.H"
#include "clock.H"
#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "polyPatchList.H"
#include "processorPolyPatch.H"
#include "labelList.H"
#include "boolList.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "OFstream.H"
#include "IFstream.H"
#include "OPstream.H"
#include "IPstream.H"
#include <mpi.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//void write(const word& name, const scalarField& elems)
//{
//    Pout<< name << endl;
//
//    forAll(elems, elemI)
//    {
//        Pout<< "    " << elemI << "  " << elems[elemI] << endl;
//    }
//}
//
//    forAll(patches, patchI)
//    {
//        write
//        (
//            word("p") + "_" + patches[patchI].name(),
//            p.boundaryField()[patchI]
//        );
//    }
//    Pout<< endl << endl;


void send
(
    const bool blockingSend,
    const label procPatchI, 
    const processorPolyPatch& procPatch,
    scalarField& sendBuf,
    MPI_Request& request
)
{
    bool transferFailed;

//    Pout<< "ProcPatch:" << procPatchI
//        << "  sending:" << sendBuf.size() * sizeof(scalar)
//        << " bytes" << endl;

    if (blockingSend)
    {
        transferFailed = MPI_Bsend
        (
            sendBuf.begin(),
            sendBuf.size() * sizeof(scalar),
            MPI_PACKED,
            Pstream::procID(procPatch.neighbProcNo()),
            Pstream::msgType(),
            MPI_COMM_WORLD
        );
    }
    else
    {
        transferFailed = MPI_Isend
        (
            sendBuf.begin(),
            sendBuf.size() * sizeof(scalar),
            MPI_PACKED,
            Pstream::procID(procPatch.neighbProcNo()),
            Pstream::msgType(),
            MPI_COMM_WORLD,
            &request
        );
    }

    if (transferFailed)
    {
        FatalErrorIn("OPstream::~OPstream()")
            << "MPI_Send cannot send outgoing message";
        Pstream::abort();
    }
//    Pout<< "ProcPatch:" << procPatchI
//        << "  finished sending:" << sendBuf.size() * sizeof(scalar)
//        << " bytes" << endl;
}



void recv
(
    const bool blockingRecv,
    const label procPatchI,
    const processorPolyPatch& procPatch,
    scalarField& recvBuf,
    MPI_Request& request
)
{
    bool transferFailed;

//    Pout<< "ProcPatch:" << procPatchI
//        << "  receiving:" << recvBuf.size() * sizeof(scalar)
//        << " bytes" << endl;

    if (blockingRecv)
    {
        MPI_Status status;

        transferFailed = MPI_Recv
        (
            recvBuf.begin(),
            recvBuf.size() * sizeof(scalar),
            MPI_PACKED,
            Pstream::procID(procPatch.neighbProcNo()),
            Pstream::msgType(),
            MPI_COMM_WORLD,
            &status
        );

        int messageSize;

        MPI_Get_count(&status, MPI_BYTE, &messageSize);

        if (messageSize > (recvBuf.size() * sizeof(scalar)))
        {
            FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
                << "buffer (" << recvBuf.size()
                << ") not large enough for incomming message ("
                << messageSize<< ')';
            Pstream::abort();
        }
    }
    else
    {
        transferFailed = MPI_Irecv
        (
            recvBuf.begin(),
            recvBuf.size() * sizeof(scalar),
            MPI_PACKED,
            Pstream::procID(procPatch.neighbProcNo()),
            Pstream::msgType(),
            MPI_COMM_WORLD,
            &request
        );
    }

    if (transferFailed)
    {
        FatalErrorIn("OPstream::~OPstream()")
            << "MPI_Irecv cannot receive incoming message";
        Pstream::abort();
    }

//    Pout<< "ProcPatch:" << procPatchI
//        << "  finished receiving:"
//        << recvBuf.size() * sizeof(scalar)
//        << " bytes" << endl;
}


// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::validArgs.append("number of swaps");
    Foam::argList::validOptions.insert("blockingSend", "");
    Foam::argList::validOptions.insert("blockingRecv", "");

#   include "setRootCase.H"

    label nSwaps(readLabel(IStringStream(args.args()[3])()));

    bool blockingSend = args.options().found("blockingSend");
    bool blockingRecv = args.options().found("blockingRecv");

    Info<< "Doing " << endl
        << "    swaps   :" << nSwaps << endl
        << "    blockingSend:" << blockingSend << endl
        << "    blockingRecv:" << blockingRecv << endl
        << endl;

    Info<< "To set network buffer size:" << endl
        << "    echo $MPI_BUFFER_SIZE > /proc/sys/net/core/rmem_max" << endl
        << "    echo $MPI_BUFFER_SIZE > /proc/sys/net/core/wmem_max" << endl << endl;



#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

    Info<< "Mesh and fields read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;

    processorTopology topo(mesh.boundaryMesh());

    commSchedule commIters(topo);

    if (Pstream::master())
    {
        fileName fName("schedule");

        Info<< endl << "Writing schedule to " << fName << endl << endl;

        OFstream scheduleStream(fName);

        scheduleStream << commIters;
    }

    const polyPatchList& patches = mesh.boundaryMesh();

    //
    // Count number of processor patches and fill patch addressing arrays.
    //

    // List of all patch numbers that are processor patches
    // (i.e. processorPatch to patch number)
    labelList procPatches(patches.size());

    // -1 or procPatch number (i.e. from patch to processorPatch)
    labelList procPatchNo(patches.size(), -1);

    label nProcPatches = 0;

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (patch.type() == processorPolyPatch::typeName)
        {
            procPatches[nProcPatches] = patchI;

            procPatchNo[patchI] = nProcPatches;

            nProcPatches++;
        }
    }
    procPatches.setSize(nProcPatches);



    //
    // Allocate send and receive buffers.
    //


    List<scalarField> sendBuf(nProcPatches);
    List<scalarField> recvBuf(nProcPatches);
    // send and recv requests stored in same list.
    List<MPI_Request> requests(2*nProcPatches);

    forAll(procPatches, procPatchI)
    {
        label patchI = procPatches[procPatchI];

        sendBuf[procPatchI].setSize(patches[patchI].size());
        recvBuf[procPatchI].setSize(patches[patchI].size());

        Pout<< "Patch:" << patchI << " procPatch:" << procPatchI
            << "  size:" << patches[patchI].size() << endl;
    }


    Info<< "Allocated buffers in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;


    Foam::clock wallTimer;

    for (label iter = 0; iter < nSwaps; iter++)
    {
        const labelList& schedule = commIters[Pstream::myProcNo()];

        forAll(schedule, iterI)
        {
            label nb = schedule[iterI];

            label patchI = topo.procPatchMap()[nb];

            label procPatchI = procPatchNo[patchI];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            if (nb < Pstream::myProcNo())
            {
                sendBuf[procPatchI] = p.boundaryField()[patchI];

                send
                (
                    blockingSend,
                    procPatchI,
                    procPatch,
                    sendBuf[procPatchI],
                    requests[procPatchI]
                );

                recv
                (
                    blockingRecv,
                    procPatchI,
                    procPatch,
                    recvBuf[procPatchI],
                    requests[procPatchI + nProcPatches]
                );
            }
            else
            {
                recv
                (
                    blockingRecv,
                    procPatchI,
                    procPatch,
                    recvBuf[procPatchI],
                    requests[procPatchI + nProcPatches]
                );

                sendBuf[procPatchI] = p.boundaryField()[patchI];

                send
                (
                    blockingSend,
                    procPatchI,
                    procPatch,
                    sendBuf[procPatchI],
                    requests[procPatchI]
                );
            }
        }


        //
        // Wait for all outstanding requests to finish.
        //

        List<MPI_Status> stats(requests.size());

        if (blockingSend)
        {
            if (blockingRecv)
            {
                // All is done blocking. No outstanding requests
            }
            else
            {
                // Wait for outstanding receives.
                if
                (
                    MPI_Waitall
                    (
                        nProcPatches,       // slice of receive requests only
                        &requests[nProcPatches],
                        &stats[nProcPatches]
                    )
                )
                {
                    FatalErrorIn("OPstream::~OPstream()")
                        << "MPI_Waitall error";
                    Pstream::abort();
                }
            }
        }
        else
        {
            if (blockingRecv)
            {
                // Wait for outstanding send requests
                if
                (
                    MPI_Waitall
                    (
                        nProcPatches,       // slice of send requests only
                        requests.begin(),
                        stats.begin()
                    )
                )
                {
                    FatalErrorIn("OPstream::~OPstream()")
                        << "MPI_Waitall error";
                    Pstream::abort();
                }
            }
            else
            {
                // Wait for outstanding send and receives
                if 
                (
                    MPI_Waitall
                    (
                        requests.size(),
                        requests.begin(),
                        stats.begin()
                    )
                )
                {
                    FatalErrorIn("OPstream::~OPstream()")
                        << "MPI_Waitall error";
                    Pstream::abort();
                }
            }
        }

        //
        // Check receive status and copy receive buffers into fieldd
        //

        forAll(procPatches, procPatchI)
        {
            label patchI = procPatches[procPatchI];

            // Check status
            if (!blockingRecv)
            {
                int messageSize;

                label recPatchI = procPatchI + nProcPatches;

                MPI_Get_count(&stats[recPatchI], MPI_BYTE, &messageSize);

                if (messageSize > (recvBuf[procPatchI].size() * sizeof(scalar)))
                {
                    FatalErrorIn("IPstream::IPstream(const int fromProcNo)")
                        << "buffer (" << recvBuf[procPatchI].size()
                        << ") not large enough for incomming message ("
                        << messageSize<< ')' << abort(FatalError);
                }
            }

            // Copy from receive buffer
            p.boundaryField()[patchI] = recvBuf[procPatchI];
        }
    }

    Info<< "Done"
        << " swaps:" << nSwaps
        << " blockingSend:" << blockingSend
        << " blockingRecv:" << blockingRecv
        << " user time:" << runTime.cpuTimeIncrement()
        << " wall time:" << wallTimer.elapsedClockTime()
        << endl << endl;



    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
