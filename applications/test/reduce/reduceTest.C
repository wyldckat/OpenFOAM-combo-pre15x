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

Application
    icoFoam

Description
    Incompressible laminar CFD code.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#define PstreamReduceOpsMPI_H
#include "PstreamReduceOps.H"
#include "PstreamReduceOpsMPINew.H"
#include "IPstream.H"
#include "OPstream.H"
#include "cpuTime.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList args(argc, argv);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting transfers\n" << endl;

    cpuTime execTime;

    scalar sum = 0;
    for (int i=0; i<100000; i++)
    {
        scalar data = i;
        reduce(data, sumOp<scalar>());
        sum+= data;
    }
    Info<< "\nResult of summation " << sum << endl;

    Info<< "CPU time " << execTime.elapsedCpuTime() << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
