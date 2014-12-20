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
    Could work in parallel but never tested.

\*---------------------------------------------------------------------------*/

#include "SortableList.H"
#include "argList.H"
#include "polyMesh.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("largestOnly", "");

#   include "setRootCase.H"
#   include "createTime.H"

    bool largestOnly = args.options().found("largestOnly");

    Info<< "Create mesh\n" << endl;

    fvMeshSubset mesh
    (
        IOobject
        (
            fvMeshSubset::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );
    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;


    // Determine connected regions
    regionSplit rs(mesh);

    Info<< endl << "Number of regions:" << rs.nRegions() << nl << endl;

    if (rs.nRegions() == 1)
    {
        Info<< "Only one region. Doing nothing." << endl;
    }
    else
    {
        const labelList& region = rs.cellToRegion();

        SortableList<label> regionSizes(rs.nRegions(), 0);

        forAll(region, cellI)
        {
            regionSizes[region[cellI]]++;
        }
        forAll(regionSizes, regionI)
        {
            reduce(regionSizes[regionI], sumOp<label>());
        }

        // So now all processors have same data. Sort in ascending order
        // so largest is last. Walk through in reverse order so largest is
        // first (cannot reverse list since would not reverse
        // SortableList::indices)

        regionSizes.sort();

        label regionI = 0;

        for (label i = rs.nRegions()-1; i >= 0; i--)
        {
            label oldRegionI = regionSizes.indices()[i];

            Info<< "Region " << regionI << nl
                << "-------- " << endl;

            Info<< "Region " << regionI
                << " mesh has " << regionSizes[i]
                << " cells" << endl;

            // Subset. Create "oldInternalFaces" patch.
            mesh.setLargeCellSubset(region, oldRegionI, -1);
 
            Info<< "Mesh subset in = "
                << runTime.cpuTimeIncrement() << " s\n" << endl;

            runTime++;

            Info<< "Writing region " << regionI << " polyMesh to time "
                << runTime.timeName() << endl;

            mesh.subMesh().write();

            Info<< endl;

            if (largestOnly)
            {
                break;
            }
            regionI++;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
