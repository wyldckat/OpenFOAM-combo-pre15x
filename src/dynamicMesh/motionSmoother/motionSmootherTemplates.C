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

#include "motionSmoother.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Average of connected points. Unweighted.
template <class Type>
Type Foam::motionSmoother::avg
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const edgeList& edges,
    const pointField& points,
    const labelList& edgeLabels,
    const label pointI
)
{
    // Get avg scale of neighbours
    Type avg(pTraits<Type>::zero);

    forAll(edgeLabels, i)
    {
        const edge& e = edges[edgeLabels[i]];

        avg += fld[e.otherVertex(pointI)];
    }

    return avg/edgeLabels.size();
}


// Distance weighted average with varying diffusivity.
template <class Type>
Type Foam::motionSmoother::avg
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const scalarField& edgeGamma,
    const edgeList& edges,
    const pointField& points,
    const labelList& edgeLabels,
    const label pointI
)
{
    // Get avg scale of neighbours
    Type avg(pTraits<Type>::zero);

    scalar sumWeight = 0.0;

    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];

        const edge& e = edges[edgeI];

        scalar weight = min(GREAT, edgeGamma[edgeI] / Foam::sqr(e.mag(points)));

        avg += weight*fld[e.otherVertex(pointI)];
        sumWeight += weight;
    }

    return avg/sumWeight;
}


// smooth field (Gauss Seidel i.e. inline update)
template <class Type>
void Foam::motionSmoother::smooth
(
    GeometricField<Type, pointPatchField, pointMesh>& fld
) const
{
    const polyMesh& mesh = fld.mesh()();

    const labelListList& pointEdges = mesh.pointEdges();
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    forAll(fld, pointI)
    {
        if (isInternalPoint(pointI))
        {
            fld[pointI] =
                0.5*fld[pointI]
              + 0.5*avg(fld, edges, points, pointEdges[pointI], pointI);
        }
    }
    fld.correctBoundaryConditions();    
}


// smooth field (point-jacobi)
template <class Type>
void Foam::motionSmoother::smooth
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const scalarField& edgeGamma,
    GeometricField<Type, pointPatchField, pointMesh>& newFld
) const
{
    const polyMesh& mesh = fld.mesh()();

    const labelListList& pointEdges = mesh.pointEdges();
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    forAll(fld, pointI)
    {
        if (isInternalPoint(pointI))
        {
            newFld[pointI] =
                0.5*fld[pointI]
              + 0.5*avg
                    (
                        fld,
                        edgeGamma,
                        edges,
                        points,
                        pointEdges[pointI],
                        pointI
                    );
        }
    }
    newFld.correctBoundaryConditions();    
}



// ************************************************************************* //
