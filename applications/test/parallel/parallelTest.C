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
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Serr<< "\nStarting transfers\n" << endl;

    vector data(0, 1, 2);

    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            {
                Serr<< "slave sending to master "
                    << Pstream::masterNo() << endl;
                OPstream toMaster(Pstream::masterNo());
                toMaster << data;
            }

            Serr<< "slave receiving from master " 
                << Pstream::masterNo() << endl;
            IPstream fromMaster(Pstream::masterNo());
            fromMaster >> data;

            Serr<< data << endl;
        }
        else
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Serr << "master receiving from slave " << slave << endl;
                IPstream fromSlave(slave);
                fromSlave >> data;

                Serr<< data << endl;
            }

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Serr << "master sending to slave " << slave << endl;
                OPstream toSlave(slave);
                toSlave << data;
            }
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
