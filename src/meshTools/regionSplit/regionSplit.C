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
    This class separates the mesh into distinct unconnected regions,
    each of which is then given a label.

\*---------------------------------------------------------------------------*/

#include "regionSplit.H"
#include "meshWave.H"
#include "regionInfo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionSplit::fillSeedMask
(
    labelList& region,
    const label seedCellID,
    const label mark
) const
{
    // Set faces on unset cell
    const labelList& myFaces = mesh_.cells()[seedCellID];

    labelList changedFaces(myFaces.size());
    List<regionInfo> changedFacesInfo(myFaces.size());

    forAll(myFaces, myFaceI)
    {
        changedFaces[myFaceI] = myFaces[myFaceI];

        changedFacesInfo[myFaceI] = regionInfo(mark);
    }

    // Propagate information over whole domain.

    meshWave<regionInfo> regionCalc
    (
        mesh_,
        changedFaces,
        changedFacesInfo,
        mesh_.nCells()  // max iterations
    );

    // Now regionCalc should hold info on cells that are reachable from
    // changedFaces -> these are part of the current region. Merge these
    // into the overall region array.
    const List<regionInfo>& cellInfo = regionCalc.allCellInfo();

    forAll(cellInfo, cellI)
    {
        if (cellInfo[cellI].valid())
        {
            region[cellI] = mark;
        }
    }
}


void Foam::regionSplit::calcRegionSplit() const
{
    // Start with regionI
    nRegions_ = 0;

    // Region per cell.
    cellToRegionPtr_ = new labelList(mesh_.nCells(), -1);
    labelList& region = *cellToRegionPtr_;
    
    do
    {
        // Find first unset cell

        label unsetCellI = -1;

        forAll(region, cellI)
        {
            if (region[cellI] == -1)
            {
                unsetCellI = cellI;

                break; 
            }
        }

        if (unsetCellI == -1)
        {
            break;
        }

        fillSeedMask(region, unsetCellI, nRegions_);

        // Go to next region
        nRegions_++;
    }
    while(true);

    // If the number of regions is one, there is no need to store the
    // decomposition
    if (nRegions_ == 1)
    {
        deleteDemandDrivenData(cellToRegionPtr_);
    }
}


Foam::labelList Foam::regionSplit::seedMask(const label seedCellID) const
{
    labelList mask(mesh_.nCells(), 0);

    fillSeedMask(mask, seedCellID, 1);

    return mask;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::regionSplit::nRegions() const
{
    if (nRegions_ < 0)
    {
        calcRegionSplit();
    }

    return nRegions_;
}


const Foam::labelList& Foam::regionSplit::cellToRegion() const
{
    // The if-statement forces the calculation of the split
    if (nRegions() == 1)
    {
        // Only 1 region present.  Allocate uniform array
        cellToRegionPtr_ = new labelList(mesh_.nCells(), 0);
    }

    return *cellToRegionPtr_;
}

// ************************************************************************* //
