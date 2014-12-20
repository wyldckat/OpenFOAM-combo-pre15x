/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Routines to check for point-nearness.

\*---------------------------------------------------------------------------*/

#include "checkCoords.H"
#include "pointSet.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Merge points and create mapping array. Return mapping.
bool checkCoords
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar reportDistSqr,
    labelHashSet* setPtr
)
{
    const pointField& pts = mesh.points();

    // Sort points
    SortableList<scalar> sortedMag(magSqr(pts));

    label nClose = 0;

    for (label i = 1; i < sortedMag.size(); i++)
    {
        label ptI = sortedMag.indices()[i];

        // Compare ptI to any previous points with similar sortedMag
        for
        (
            label j = i-1;
            j >= 0 && (sortedMag[j] > sortedMag[i]-reportDistSqr);
            --j
        )
        {
            label prevPtI = sortedMag.indices()[j];

            if (magSqr(pts[ptI] - pts[prevPtI]) < reportDistSqr)
            {
                // Check if unconnected.
                const labelList& pEdges = mesh.pointEdges()[ptI];

                bool connected = false;

                forAll(pEdges, pEdgeI)
                {
                    if (mesh.edges()[pEdges[pEdgeI]].otherVertex(prevPtI) != -1)
                    {
                        connected = true;
                        break;
                    }
                }

                if (!connected)
                {
                    nClose++;

                    if (report)
                    {
                        Pout<< "checkCoords : unconnected points "
                            << ptI << " and " << prevPtI
                            << " with coordinates:" << pts[ptI]
                            << " and " << pts[prevPtI] << endl
                            << " are relatively close to each other."
                            << " This might be correct in case of e.g. baffles."
                            << endl;
                    }

                    if (setPtr)
                    {
                        setPtr->insert(ptI);
                        setPtr->insert(prevPtI);
                    }
                }
            }
        }
    }

    reduce(nClose, sumOp<label>());

    if (nClose > 0)
    {
        if (Pstream::master())
        {
            WarningIn
            (
                "checkCoords(const primitiveMesh&, "
                "const bool, const scalar, labelHashSet*)"
            )   << nClose  << " points that are close together found." << nl
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
