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

Application
    test

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOmanip.H"
#include "scalar.H"
#include "List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(void)
{
/*    dictionary dotguiserc(string(getenv("HOME")) + "/.guiserc");
    int rtn = dotguiserc.read();

    // Look up debug switches in the start-up
    readDebugSwitches(dotguiserc);
*/
//    List<scalar> slist(100000);

//    register label i;
/*    for (i=0; i<100000; i++)
    {
        slist[i] = i;
    }
    for (i=0; i<10; i++)
    {
        Info<< slist[i] << endl;
    }
    Ostream data_out("/usr/tmp/data");
    data_out << slist;
    data_out.close();
*/
/*    Istream data_in("/usr/tmp/data", BINARY);
    data_in >> slist;

    for (i=0; i<10; i++)
    {
        Info<< slist[i] << endl;
    }
*/

    string st("sfdsf  sdfs23df sdf32f .  sdfsdff23/2sf32");
    Info<< word(st) << "END" << endl;

    string st1("1234567");

    Info<< st1.size() << tab << string(word(st1)) << endl;

    Info<< setw(20) << setprecision(3) << 1.234234 << endl;

    Info<< hex << 255 << endl;
}

// ************************************************************************* //

