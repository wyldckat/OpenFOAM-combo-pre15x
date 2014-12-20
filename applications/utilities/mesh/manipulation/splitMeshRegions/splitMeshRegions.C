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
    Splits mesh into multiple regions and writes them to consecutive
    time directories. Each region is defined as a domain
    whose cells can all be reached by cell-face-cell walking.
    Uses meshWave.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "meshWave.H"
#include "regionSplit.H"
#include "meshSubset.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Create mesh\n" << endl;
    meshSubset mesh
    (
        IOobject
        (
            meshSubset::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );
    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;


    // Determine connected regions
    regionSplit rs(mesh);

    Info<< endl << "Number of regions:" << rs.nRegions() << endl;

    if (rs.nRegions() == 1)
    {
        Info<< "Only one region. Doing nothing." << endl;
    }
    else
    {
        const labelList& region = rs.cellToRegion();

        for(label regionI = 0; regionI < rs.nRegions(); regionI++)
        {
            Info<< "Region " << regionI << endl;
            Info<< "-------- " << endl;

            // Estimated size of hash
            labelHashSet cellsToSubset(mesh.nCells()/rs.nRegions() + 1);

            // Insert selected cells
            forAll(region, cellI)
            {
                if (region[cellI] == regionI)
                {
                    cellsToSubset.insert(cellI);
                }
            }
            Info<< "Region " << regionI << " mesh has " << cellsToSubset.size()
                << " cells" << endl;

            mesh.setLargeCellSubset(cellsToSubset);
 
            Info<< "Mesh subset in = "
                << runTime.cpuTimeIncrement() << " s\n" << endl << endl;

            runTime++;

            Info<< "Writing region " << regionI << " polyMesh to time "
                << runTime.timeName() << endl;

            mesh.subMesh().write();

            Info<< endl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
