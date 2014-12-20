/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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

#include "checkEdges.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Merge points and create mapping array. Return mapping.
bool checkEdges
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar tol,
    labelHashSet* setPtr
)
{
    const pointField& pts = mesh.points();
    const edgeList& edges = mesh.edges();

    scalar minLen = GREAT;
    label nSmall = 0;
    scalar maxLen = -GREAT;

    forAll (edges, edgeI)
    {
        const edge& e = edges[edgeI];
        scalar magE = e.mag(pts);

        if (magE < tol)
        {
            if (report)
            {
                Pout<< "Zero size or very small edge size detected for "
                    << "edge " << edgeI << " vertices " << e
                    << ".  Length = " << magE << endl;
            }

            if (setPtr)
            {
                setPtr->insert(e[0]);
                setPtr->insert(e[1]);
            }
            nSmall++;
        }

        minLen = min(minLen, magE);
        maxLen = max(maxLen, magE);
    }

    reduce(minLen, minOp<scalar>());
    reduce(maxLen, maxOp<scalar>());
    reduce(nSmall, sumOp<label>());

    if (report)
    {
        Info<< "Minumum edge length = " << minLen
            << ". Maximum edge length = " << maxLen
            << '.' << nl << endl;
    }

    if (nSmall > 0)
    {
        if (Pstream::master())
        {
            WarningIn
            (
                "checkEdges"
                "(const primitiveMesh& mesh, const bool report,"
                "const scalar tol, labelHashSet* setPtr"
            )  << nSmall  << " small edges found" << nl << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
