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
    Calculate agglomeration clusters given a fine matrix.
    Agglomeration clusters give a coarse cluster label for each fine cell.
    Hook and return true if the number of cell on coarse level is larger than
    required minimum (nCellsInTopLevel_).

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "amgSymSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool amgSymSolver::calcAgglomeration(const label fineLevelIndex)
{
    // Algorithm:
    // 1) Create temporary cell-face addressing using a double-pass algorithm.
    //    Other options here include a guessed number of neighbours ("do-if"
    //    problem or singly-linked lists. If you are sure there is another
    //    option which is faster and more memory-efficient than mine, plase
    //    replace the cell-face search (currently it is two-pass)
    // 2) loop through all cells and for each cell find the best fit neighbour.
    //    If all neighbours are grouped, "group" the cell on its own.
    // 4) If the number of coarse cells is grater than minimum, hook onto the
    //    list and return true; otherwise return false

    // Get addressing
    const lduMatrix& fineMatrix = matrixLevel(fineLevelIndex);

    const label nEqns = fineMatrix.lduAddr().size();

    const unallocLabelList& upperAddr = fineMatrix.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = fineMatrix.lduAddr().lowerAddr();

    // Get magnitudes of off-diagonal matrix coefficients
    const scalarField magUpper = mag(fineMatrix.upper());

    // For each cell calculate faces
    labelListList cellFaces(nEqns);

    // memory management
    {
        labelList nNbrs(nEqns, 0);

        forAll (upperAddr, faceI)
        {
            nNbrs[upperAddr[faceI]]++;
        }

        forAll (lowerAddr, faceI)
        {
            nNbrs[lowerAddr[faceI]]++;
        }

        forAll (cellFaces, cellI)
        {
            cellFaces[cellI].setSize(nNbrs[cellI]);
        }

        // reset the whole list to use as counter
        nNbrs = 0;

        forAll (upperAddr, faceI)
        {
            cellFaces[upperAddr[faceI]][nNbrs[upperAddr[faceI]]] = faceI;

            nNbrs[upperAddr[faceI]]++;
        }

        forAll (lowerAddr, faceI)
        {
            cellFaces[lowerAddr[faceI]][nNbrs[lowerAddr[faceI]]] = faceI;

            nNbrs[lowerAddr[faceI]]++;
        }
    }


    // go through the faces and create clusters

    labelList coarseCellMap(nEqns, -1);

    label nCoarseCells = 0;

    forAll (cellFaces, cellI)
    {
        if (coarseCellMap[cellI] < 0)
        {
            // cell not merged. Find the best neighbour
            const labelList& curFaces = cellFaces[cellI];

            label matchFaceNo = -1;

            scalar maxFaceCoeff = -GREAT;

            // check all faces to find ungrouped neighbour with largest face
            // coefficient
            forAll (curFaces, curFaceI)
            {
                // I don't know whether the current cell is owner or neighbour.
                // Therefore I'll check both sides
                if
                (
                    coarseCellMap[upperAddr[curFaces[curFaceI]]] < 0
                 && coarseCellMap[lowerAddr[curFaces[curFaceI]]] < 0
                 && magUpper[curFaces[curFaceI]] > maxFaceCoeff
                )
                {
                    // Match found. Pick up all the necessary data
                    matchFaceNo = curFaceI;

                    maxFaceCoeff = magUpper[curFaces[curFaceI]];
                }
            }

            if (matchFaceNo >= 0)
            {
                // Make a new group
                coarseCellMap[upperAddr[curFaces[matchFaceNo]]] = nCoarseCells;
                coarseCellMap[lowerAddr[curFaces[matchFaceNo]]] = nCoarseCells;

                nCoarseCells++;
            }
            else
            {
                // No match. Find the best neighbouring cluster and
                // put the cell there
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = -GREAT;

                forAll (curFaces, curFaceI)
                {
                    if (magUpper[curFaces[curFaceI]] > maxFaceCoeff)
                    {
                        clusterMatchFaceNo = curFaceI;
                        clusterMaxFaceCoeff = magUpper[curFaces[curFaceI]];
                    }
                }

                if (clusterMatchFaceNo >= 0)
                {
                    // Add the cell to the best cluster
                    coarseCellMap[cellI] =
                        max
                        (
                            coarseCellMap
                                [upperAddr[curFaces[clusterMatchFaceNo]]],
                            coarseCellMap
                                [lowerAddr[curFaces[clusterMatchFaceNo]]]
                        );
                }
            }
        }
    }

    // reverse the map to facilitate later agglomeration

    // for easier substitutions, decrement nCoarseCells by one
    nCoarseCells--;

    forAll (coarseCellMap, cellI)
    {
        coarseCellMap[cellI] = nCoarseCells - coarseCellMap[cellI];
    }

    // The decision on parallel agglomeration needs to be made for the
    // whole gang of processes; otherwise I may end up with a different
    // number of agglomeration levels on different processors.
    // 
    bool moreAgglom = false;

    if (nCoarseCells >= nCellsInTopLevel_)
    {
        moreAgglom = true;
    }

    reduce(moreAgglom, andOp<bool>());

    if (moreAgglom)
    {
        restrictAddressing_.hook(new labelField(coarseCellMap));
    }

    return moreAgglom;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
