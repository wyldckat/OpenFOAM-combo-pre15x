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

#include "error.H"

#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void primitiveMesh::printAllocated() const
{
    Info<< "primitiveMesh allocated :" << endl;

    // Motion
    if (oldPointsPtr_)
    {
        Info<< "    Old points" << endl;
    }

    // Topology
    if (cellShapesPtr_)
    {
        Info<< "    Cell shapes" << endl;
    }

    if (edgesPtr_)
    {
        Info<< "    Edges" << endl;
    }

    if (ccPtr_)
    {
        Info<< "    Cell-cells" << endl;
    }

    if (ecPtr_)
    {
        Info<< "    Edge-cells" << endl;
    }

    if (pcPtr_)
    {
        Info<< "    Point-cells" << endl;
    }

    if (cfPtr_)
    {
        Info<< "    Cell-faces" << endl;
    }

    if (efPtr_)
    {
        Info<< "    Edge-faces" << endl;
    }

    if (pfPtr_)
    {
        Info<< "    Point-faces" << endl;
    }

    if (cePtr_)
    {
        Info<< "    Cell-edges" << endl;
    }

    if (fePtr_)
    {
        Info<< "    Face-edges" << endl;
    }

    if (pePtr_)
    {
        Info<< "    Point-edges" << endl;
    }

    if (ppPtr_)
    {
        Info<< "    Point-point" << endl;
    }
    
    if (cpPtr_)
    {
        Info<< "    Cell-point" << endl;
    }

    // Geometry
    if (cellCentresPtr_)
    {
        Info<< "    Cell-centres" << endl;
    }

    if (faceCentresPtr_)
    {
        Info<< "    Face-centres" << endl;
    }

    if (edgeVectorsPtr_)
    {
        Info<< "    Edge-vectors" << endl;
    }

    if (cellVolumesPtr_)
    {
        Info<< "    Cell-volumes" << endl;
    }

    if (faceAreasPtr_)
    {
        Info<< "    Face-areas" << endl;
    }

}


void primitiveMesh::clearGeom()
{
    if (debug)
    {
        Info<< "primitiveMesh::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    deleteDemandDrivenData(cellCentresPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(edgeVectorsPtr_);
    deleteDemandDrivenData(cellVolumesPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
}


void primitiveMesh::clearAddressing()
{
    if (debug)
    {
        Info<< "primitiveMesh::clearAddressing() : "
            << "clearing topology"
            << endl;
    }

    deleteDemandDrivenData(cellShapesPtr_);

    clearOutEdges();

    deleteDemandDrivenData(ccPtr_);
    deleteDemandDrivenData(ecPtr_);
    deleteDemandDrivenData(pcPtr_);

    deleteDemandDrivenData(cfPtr_);
    deleteDemandDrivenData(efPtr_);
    deleteDemandDrivenData(pfPtr_);

    deleteDemandDrivenData(cePtr_);
    deleteDemandDrivenData(fePtr_);
    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(ppPtr_);
    deleteDemandDrivenData(cpPtr_);
}


void primitiveMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


void primitiveMesh::clearPrimitives()
{
    deleteDemandDrivenData(oldPointsPtr_);

    clearedPrimitives_ = true;
}


void primitiveMesh::clearAll()
{
    clearGeom();
    clearAddressing();
    clearPrimitives();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
