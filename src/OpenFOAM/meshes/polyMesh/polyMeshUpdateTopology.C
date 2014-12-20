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
    Update topology of a polyMesh using the mesh morphing engine.

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "parallelInfo.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyMesh::setMorphTimeIndex(const label newTimeIndex) const
{
    if (curMorphTimeIndex_ == newTimeIndex)
    {
        return false;
    }
    else
    {
        curMorphTimeIndex_ = newTimeIndex;

        // Delete the morph map as it is out of date
        deleteDemandDrivenData(morphMap_);

        return true;
    }
}


bool Foam::polyMesh::morphing() const
{
    if (setMorphTimeIndex(time().timeIndex()))
    {
        // In a different time step.  Find out if the morphing engine requires
        // a topology change
        morphing_ = morphEngine().changeTopology();
    }

    return morphing_;
}


void Foam::polyMesh::updateTopology()
{
    if (morphing())
    {
        // If the mesh is morphing it cannot be moving simultaneously
        moving_ = false;

        morph(morphEngine().topoChangeRequest());

        // Force the mesh modifiers to auto-write.  This allows us to
        // preserve the current state of modifiers corresponding with
        // the mesh.  
        if (morphEnginePtr_)
        {
            morphEnginePtr_->writeOpt() = IOobject::AUTO_WRITE;
            morphEnginePtr_->instance() = time().timeName();
        }

        // Update zones
        pointZones_.clearAddressing();
        faceZones_.clearAddressing();
        cellZones_.clearAddressing();

        // Update parallel data
        if (parallelDataPtr_)
        {
            parallelDataPtr_->updateTopology(*morphMap_);
        }

        // Update mesh modifiers
        if (morphEnginePtr_)
        {
            morphEnginePtr_->updateTopology(*morphMap_);
        }
    }
}


void Foam::polyMesh::updateTopology(const polyTopoChange& ref)
{
    if (morphEnginePtr_)
    {
        FatalErrorIn("polyMesh::updateTopology(const polyTopoChange&)")
            << "Cannot mix meshModifiers with explicitly provided mesh changes"
            << endl
            << "Please restart without a meshModifiers file"
            << exit(FatalError);
    }

    // If the mesh is morphing it cannot be moving simultaneously
    moving_ = false;

    morph(ref);

    // Update zones
    pointZones_.clearAddressing();
    faceZones_.clearAddressing();
    cellZones_.clearAddressing();

    // Update parallel data
    if (parallelDataPtr_)
    {
        parallelDataPtr_->updateTopology(*morphMap_);
    }
}


const Foam::mapPolyMesh& Foam::polyMesh::morphMap() const
{
    if (!morphMap_)
    {
        FatalErrorIn
        (
            "const mapPolyMesh& polyMesh::morphMap() const"
        )   << "Morph map not available as the topology has not changed "
            << "within the current time-step"
            << abort(FatalError);
    }

    return *morphMap_;
}


void Foam::polyMesh::resetMorph() const
{
    curMorphTimeIndex_ = 0;
    deleteDemandDrivenData(morphMap_);
}


// ************************************************************************* //
