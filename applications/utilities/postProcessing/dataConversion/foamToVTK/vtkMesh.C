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

Description

\*---------------------------------------------------------------------------*/

#include "vtkMesh.H"
#include "vtkTopo.H"
#include "meshSubset.H"
#include "fvMesh.H"
#include "Time.H"
#include "cellSet.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


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
    meshSubset(io),
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
        setCellSubset(currentSet);
    }

    // Construct topo for mesh or subset
    topoPtr_ = new vtkTopo(mesh());
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
    polyMesh::readUpdateState meshState = meshSubset::readUpdate();

    if (meshState != polyMesh::UNCHANGED)
    {
        // Note: since meshSubset has no movePoints() functionality reconstruct
        // the subset even if only movement.

        deleteDemandDrivenData(cMapPtr_);
        deleteDemandDrivenData(pMapPtr_);
        deleteDemandDrivenData(topoPtr_);

        if (setName_.size() > 0)
        {
            Info<< "Subsetting mesh based on cellSet " << setName_ << endl;

            // Read cellSet using whole mesh
            cellSet currentSet(*this, setName_);

            setCellSubset(currentSet);
        }

        // Construct topo for subset
        topoPtr_ = new vtkTopo(mesh());
    }

    return meshState;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
