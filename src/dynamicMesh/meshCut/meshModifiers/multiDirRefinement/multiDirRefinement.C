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

#include "multiDirRefinement.H"
#include "Time.H"
#include "undoableMeshCutter.H"
#include "hexCellLooper.H"
#include "geomCellLooper.H"
#include "topoSet.H"
#include "directions.H"
#include "morphMesh.H"
#include "hexRef8.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "ListSearch.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiDirRefinement, 0);
}


// * * * * * * * * * * * * * Private Statc Functions * * * * * * * * * * * * //

// Update refCells pattern for split cells. Note refCells is current
// list of cells to refine (these should all have been refined)
void Foam::multiDirRefinement::addCells
(
    const Map<label>& splitMap,
    List<refineCell>& refCells
)
{
    label newRefI = refCells.size();

    label oldSize = refCells.size();

    refCells.setSize(newRefI + splitMap.size());

    for (label refI = 0; refI < oldSize; refI++)
    {
        const refineCell& refCell = refCells[refI];

        Map<label>::const_iterator iter = splitMap.find(refCell.cellNo());

        if (iter == splitMap.end())
        {
            FatalErrorIn
            (
                "multiDirRefinement::addCells(const Map<label>&"
                ", List<refineCell>&)"
            )   << "Problem : cannot find added cell for cell "
                << refCell.cellNo() << abort(FatalError);
        }

        refCells[newRefI++] = refineCell(iter(), refCell.direction());
    }
}


// Update vectorField for all the new cells. Takes over value of
// original cell.
void Foam::multiDirRefinement::update
(
    const Map<label>& splitMap,
    vectorField& field
)
{
    field.setSize(field.size() + splitMap.size());

    for
    (
        Map<label>::const_iterator iter = splitMap.begin();
        iter != splitMap.end();
        ++iter
    )
    {
        field[iter()] = field[iter.key()];
    }
}


// Add added cells to labelList
void Foam::multiDirRefinement::addCells
(
    const Map<label>& splitMap,
    labelList& labels
)
{
    label newCellI = labels.size();

    labels.setSize(labels.size() + splitMap.size());

    for
    (
        Map<label>::const_iterator iter = splitMap.begin();
        iter != splitMap.end();
        ++iter
    )
    {
        labels[newCellI++] = iter();
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiDirRefinement::addCells
(
    const primitiveMesh& mesh,
    const Map<label>& splitMap
)
{
    // Construct inverse addressing: from new to original cell.
    labelList origCell(mesh.nCells(), -1);

    forAll(addedCells_, cellI)
    {
        const labelList& added = addedCells_[cellI];

        forAll(added, i)
        {
            label slave = added[i];

            if (origCell[slave] == -1)
            {
                origCell[slave] = cellI;
            }
            else if (origCell[slave] != cellI)
            {
                FatalErrorIn
                (
                    "multiDirRefinement::addCells(const primitiveMesh&"
                    ", const Map<label>&"
                )   << "Added cell " << slave << " has two different masters:"
                    << origCell[slave] << " , " << cellI
                    << abort(FatalError);
            }
        }
    }


    for
    (
        Map<label>::const_iterator iter = splitMap.begin();
        iter != splitMap.end();
        ++iter
    )
    {
        label masterI = iter.key();
        label newCellI = iter();

        while (origCell[masterI] != -1 && origCell[masterI] != masterI)
        {
            masterI = origCell[masterI];
        }

        if (masterI >= addedCells_.size())
        {
            FatalErrorIn
            (
                "multiDirRefinement::addCells(const primitiveMesh&"
                ", const Map<label>&"
            )   << "Map of added cells contains master cell " << masterI
                << " which is not a valid cell number" << endl
                << "This means that the mesh is not consistent with the"
                << " done refinement" << endl
                << "newCell:" << newCellI << abort(FatalError);
        }

        labelList& added = addedCells_[masterI];

        if (added.size() == 0)
        {
            added.setSize(2);
            added[0] = masterI;
            added[1] = newCellI;
        }
        else if (findIndex(added, newCellI) == -1)
        {
            label sz = added.size();
            added.setSize(sz + 1);
            added[sz] = newCellI;
        }
    }
}


Foam::labelList Foam::multiDirRefinement::splitOffHex(const primitiveMesh& mesh)
{
    const cellModel& hex = *(cellModeller::lookup("hex"));

    const cellShapeList& cellShapes = mesh.cellShapes();

    // Split cellLabels_ into two lists.

    labelList nonHexLabels(cellLabels_.size());
    label nonHexI = 0;

    labelList hexLabels(cellLabels_.size());
    label hexI = 0;

    forAll(cellLabels_, i)
    {
        label cellI = cellLabels_[i];

        if (cellShapes[cellI].model() == hex)
        {
            hexLabels[hexI++] = cellI;
        }
        else
        {
            nonHexLabels[nonHexI++] = cellI;
        }
    }

    nonHexLabels.setSize(nonHexI);

    cellLabels_.transfer(nonHexLabels);

    hexLabels.setSize(hexI);

    return hexLabels;
}


void Foam::multiDirRefinement::refineHex8
(
    morphMesh& mesh,
    const labelList& hexCells,
    const bool writeMesh
)
{
    if (debug)
    {
        Info<< "multiDirRefinement : Refining hexes " << hexCells.size()
            << endl;
    }

    hexRef8 hexRefiner(mesh);

    polyTopoChange meshMod(mesh);

    hexRefiner.setRefinement(hexCells, meshMod);

    if (writeMesh)
    {
        ((Time&)mesh.time())++;
        mesh.setMorphTimeIndex(mesh.time().timeIndex());
    }
    else
    {
        mesh.resetMorph();
    }

    mesh.updateTopology(meshMod);

    mesh.movePoints(mesh.morphMap().preMotionPoints());

    if (writeMesh)
    {
        mesh.write();
    }

    if (debug)
    {
        Info<< "multiDirRefinement : updated mesh at time "
            << mesh.time().timeName() << endl;
    }

    // Take over split pattern from hex refiner. (should be empty at this
    // point)
    addedCells_ = hexRefiner.addedCells();
}


void Foam::multiDirRefinement::refineAllDirs
(
    morphMesh& mesh,
    List<vectorField>& cellDirections,
    const cellLooper& cellWalker,
    undoableMeshCutter& cutter,
    const bool writeMesh
)
{
    // Iterator
    refinementIterator refiner(mesh, cutter, cellWalker, writeMesh);

    forAll(cellDirections, dirI)
    {
        if (debug)
        {
            Info<< "multiDirRefinement : Refining " << cellLabels_.size()
                << " cells in direction " << dirI << endl
                << endl;
        }

        const vectorField& dirField = cellDirections[dirI];

        // Combine cell to be cut with direction to cut in.
        // If dirField is only one element use this for all cells.

        List<refineCell> refCells(cellLabels_.size());

        if (dirField.size() == 1)
        {
            // Uniform directions.
            if (debug)
            {
                Info<< "multiDirRefinement : Uniform refinement:"
                    << dirField[0] << endl;
            }

            forAll(refCells, refI)
            {
                label cellI = cellLabels_[refI];

                refCells[refI] = refineCell(cellI, dirField[0]);
            }
        }
        else
        {
            // Non uniform directions.
            forAll(refCells, refI)
            {
                label cellI = cellLabels_[refI];

                refCells[refI] = refineCell(cellI, dirField[cellI]);
            }
        }

        // Do refine mesh (multiple iterations). Remember added cells.
        Map<label> splitMap = refiner.setRefinement(refCells);

        // Update overall map of added cells
        addCells(mesh, splitMap);

        // Add added cells to list of cells to refine in next iteration
        addCells(splitMap, cellLabels_);

        // Update refinement direction for added cells.
        if (dirField.size() != 1)
        {
            forAll(cellDirections, i)
            {
                update(splitMap, cellDirections[i]);
            }
        }

        if (debug)
        {
            Info<< "multiDirRefinement : Done refining direction " << dirI
                << " resulting in " << cellLabels_.size() << " cells" << nl
                << endl;
        }
    }
}


void Foam::multiDirRefinement::refineFromDict
(
    morphMesh& mesh,
    List<vectorField>& cellDirections,
    const dictionary& dict,
    const bool writeMesh
)
{
    // How to walk cell circumference.
    Switch pureGeomCut(dict.lookup("geometricCut"));

    autoPtr<cellLooper> cellWalker(NULL);
    if (pureGeomCut)
    {
        cellWalker.reset(new geomCellLooper(mesh));
    }
    else
    {
        cellWalker.reset(new hexCellLooper(mesh));
    }

    // Construct undoable refinement topology modifier. 
    //Note: undoability switched off since hexRef8 not undoable anyway.
    // Might want to reconsider if needs to be possible. But then can always
    // use other constructor.
    undoableMeshCutter cutter(mesh, false);

    refineAllDirs(mesh, cellDirections, cellWalker(), cutter, writeMesh);
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiDirRefinement::multiDirRefinement
(
    morphMesh& mesh,
    const labelList& cellLabels,        // list of cells to refine
    const dictionary& dict
)
:
    cellLabels_(cellLabels),
    addedCells_(mesh.nCells())
{
    Switch useHex(dict.lookup("useHexTopology"));

    Switch writeMesh(dict.lookup("writeMesh"));

    wordList dirNames(dict.lookup("directions"));

    if (useHex && dirNames.size() == 3)
    {
        // Filter out hexes from cellLabels_
        labelList hexCells(splitOffHex(mesh));

        refineHex8(mesh, hexCells, writeMesh);
    }

    label nRemainingCells = cellLabels_.size();

    reduce(nRemainingCells, sumOp<label>());

    if (nRemainingCells > 0)
    {
        // Any cells to refine using meshCutter

        // Determine directions for every cell. Can either be uniform
        // (size = 1) or non-uniform (one vector per cell)
        directions cellDirections(mesh, dict);

        refineFromDict(mesh, cellDirections, dict, writeMesh);
    }
}


// Construct from directionary and directions to refine.
Foam::multiDirRefinement::multiDirRefinement
(
    morphMesh& mesh,
    const labelList& cellLabels,        // list of cells to refine
    const List<vectorField>& cellDirs,  // Explicitly provided directions
    const dictionary& dict
)
:
    cellLabels_(cellLabels),
    addedCells_(mesh.nCells())
{
    Switch useHex(dict.lookup("useHexTopology"));

    Switch writeMesh(dict.lookup("writeMesh"));

    wordList dirNames(dict.lookup("directions"));

    if (useHex && dirNames.size() == 3)
    {
        // Filter out hexes from cellLabels_
        labelList hexCells(splitOffHex(mesh));

        refineHex8(mesh, hexCells, writeMesh);
    }

    label nRemainingCells = cellLabels_.size();

    reduce(nRemainingCells, sumOp<label>());

    if (nRemainingCells > 0)
    {
        // Any cells to refine using meshCutter

        // Working copy of cells to refine
        List<vectorField> cellDirections(cellDirs);

        refineFromDict(mesh, cellDirections, dict, writeMesh);
    }
}


// Construct from components. Implies meshCutter since directions provided.
Foam::multiDirRefinement::multiDirRefinement
(
    morphMesh& mesh,
    undoableMeshCutter& cutter,     // actual mesh modifier
    const cellLooper& cellWalker,   // how to cut a single cell with
                                    // a plane
    const labelList& cellLabels,    // list of cells to refine
    const List<vectorField>& cellDirs,
    const bool writeMesh            // write intermediate meshes
)
:
    cellLabels_(cellLabels),
    addedCells_(mesh.nCells())
{
    // Working copy of cells to refine
    List<vectorField> cellDirections(cellDirs);

    refineAllDirs(mesh, cellDirections, cellWalker, cutter, writeMesh);
}


// ************************************************************************* //
