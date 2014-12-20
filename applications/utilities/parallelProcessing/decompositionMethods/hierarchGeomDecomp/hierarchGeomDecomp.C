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

#include "hierarchGeomDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hierarchGeomDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        hierarchGeomDecomp,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        hierarchGeomDecomp,
        dictionaryMesh
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::hierarchGeomDecomp::setDecompOrder()
{
    word order(geomDecomDict_.lookup("order"));

    if (order.size() != 3)
    {
        FatalIOErrorIn
        (
            "hierarchGeomDecomp::hierarchGeomDecomp"
            "(const dictionary& decompositionDict)",
            decompositionDict_
        )   << "number of characters in order (" << order << ") != 3"
            << exit(FatalIOError);
    }

    for (label i = 0; i < 3; i++)
    {
        if (order[i] == 'x')
        {
            decompOrder_[i] = 0;
        }
        else if (order[i] == 'y')
        {
            decompOrder_[i] = 1;
        }
        else if (order[i] == 'z')
        {
            decompOrder_[i] = 2;
        }
        else
        {
            FatalIOErrorIn
            (
                "hierarchGeomDecomp::hierarchGeomDecomp"
                "(const dictionary& decompositionDict)",
                decompositionDict_
            )   << "Illegal decomposition order " << order << endl
                << "It should only contain x, y or z" << exit(FatalError);
        }
    }
}


// Sort points into bins according to one component. Recurses to next component.
void Foam::hierarchGeomDecomp::sortComponent
(
    const pointField& points,
    const labelList& current,       // slice of points to decompose
    const direction componentIndex, // index in decompOrder_
    const label mult,               // multiplication factor for finalDecomp
    labelList& finalDecomp
)
{
    // Current component
    label compI = decompOrder_[componentIndex];

    //Info<< "Sorting slice of size " << current.size() << " in component "
    //    << compI << endl;

    // Storage for sorted component compI
    SortableList<scalar> sortedCoord(current.size());

    forAll(current, i)
    {
        label pointI = current[i];

        sortedCoord[i] = points[pointI][compI];
    }
    sortedCoord.sort();
    

    // Index in current.
    label index = 0;

    // Sort bins of size n
    for (label bin = 0; bin < n_[compI]; bin++)
    {
        // Size of bin
        label dx = label(current.size()/n_[compI]);

        if (bin == n_[compI]-1)
        {
            // Make final slice a bit bigger
            dx = current.size() - index;
        }
        labelList slice(dx);

        forAll(slice, i)
        {
            label pointI = current[sortedCoord.indices()[index++]];

            // Mark point into correct bin
            finalDecomp[pointI] += bin*mult;

            // And collect for next sorting action
            slice[i] = pointI;
        }

        // Sort slice in next component
        if (componentIndex < 2)
        {
            sortComponent
            (
                points,
                slice,
                componentIndex+1,
                mult*n_[compI],     // Multiplier to apply to decomposition.
                finalDecomp
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hierarchGeomDecomp::hierarchGeomDecomp
(
    const dictionary& decompositionDict
)
:
    geomDecomp(decompositionDict, typeName),
    decompOrder_()
{
    setDecompOrder();
}


Foam::hierarchGeomDecomp::hierarchGeomDecomp
(
    const dictionary& decompositionDict,
    const primitiveMesh&
)
:
    geomDecomp(decompositionDict, typeName),
    decompOrder_()
{
    setDecompOrder();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hierarchGeomDecomp::decompose(const pointField& points)
{
    // construct a list for the final result
    labelList finalDecomp(points.size(), 0);

    // Start off with every point sorted onto itself.
    labelList slice(points.size());
    forAll(slice, i)
    {
        slice[i] = i;
    }

    pointField rotatedPoints = rotDelta_ & points;

    // Sort recursive
    sortComponent
    (
        rotatedPoints,
        slice,
        0,              // Sort first component in decompOrder.
        1,              // Offset for different x bins.
        finalDecomp
    );

    return finalDecomp;
}


// ************************************************************************* //
