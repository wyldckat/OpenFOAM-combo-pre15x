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

#include "vtkMesh.H"
#include "Time.H"
#include "cellSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkMesh::calcMapping() const
{
    cMapPtr_ = new labelList(nCells());
    labelList& cMap = *cMapPtr_;

    forAll(cMap, i)
    {
        cMap[i] = i;
    }

    pMapPtr_ = new labelList(nPoints());
    labelList& pMap = *pMapPtr_;

    forAll(pMap, i)
    {
        pMap[i] = i;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::vtkMesh::vtkMesh
(
    const IOobject& io,
    const word& setName
)
:
    fvMeshSubset(io),
    setName_(setName),
    cMapPtr_(NULL),
    pMapPtr_(NULL),
    topoPtr_(NULL)
{
    if (setName.size() > 0)
    {
        // Read cellSet using whole mesh
        cellSet currentSet(*this, setName_);

        // Set current subset
        setLargeCellSubset(currentSet);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkMesh::~vtkMesh()
{
    deleteDemandDrivenData(cMapPtr_);
    deleteDemandDrivenData(pMapPtr_);
    deleteDemandDrivenData(topoPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polyMesh::readUpdateState Foam::vtkMesh::readUpdate()
{
    polyMesh::readUpdateState meshState = fvMeshSubset::readUpdate();

    if (meshState != polyMesh::UNCHANGED)
    {
        // Note: since fvMeshSubset has no movePoints() functionality reconstruct
        // the subset even if only movement.

        deleteDemandDrivenData(cMapPtr_);
        deleteDemandDrivenData(pMapPtr_);
        deleteDemandDrivenData(topoPtr_);

        if (setName_.size() > 0)
        {
            Pout<< "Subsetting mesh based on cellSet " << setName_ << endl;

            // Read cellSet using whole mesh
            cellSet currentSet(*this, setName_);

            setLargeCellSubset(currentSet);
        }
    }

    return meshState;
}


// ************************************************************************* //
