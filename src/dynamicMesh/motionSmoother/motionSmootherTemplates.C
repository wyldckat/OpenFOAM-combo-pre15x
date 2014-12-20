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

#include "motionSmoother.H"
#include "meshTools.H"
#include "processorPointPatchFields.H"
#include "globalProcessorPointPatchFields.H"

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


//- Sychronizes patch points on pointField
template<class Type, class CombineOp>
void Foam::motionSmoother::syncField
(
    GeometricField<Type, pointPatchField, pointMesh>& fld,
    const Type& zero,
    const CombineOp& cop
)
{
    const polyMesh& mesh = fld.mesh()();


    typedef ProcessorPointPatchField
        <pointPatchField, pointPatch, processorPointPatch, Type>
        ProcessorField;

    if (Pstream::parRun() && mesh.parallelData().parallel())
    {
        forAll(fld.boundaryField(), patchI)
        {
            if
            (
                isA<ProcessorField>
                (
                    fld.boundaryField()[patchI]
                )
            )
            {
                const ProcessorField& pFld =
                    refCast<const ProcessorField>
                    (
                        fld.boundaryField()[patchI]
                    );

                const processorPointPatch& procPatch =
                    refCast<const processorPointPatch>
                    (
                        pFld.patch()
                    );

                OPstream toNeighbProc
                (
                    procPatch.neighbProcNo(),
                    pFld.size()*sizeof(Type)
                );

                toNeighbProc << pFld.patchInternalField();
            }
        }

        forAll(fld.boundaryField(), patchI)
        {
            if
            (
                isA<ProcessorField>
                (
                    fld.boundaryField()[patchI]
                )
            )
            {
                ProcessorField& pFld =
                    refCast<ProcessorField>
                    (
                        fld.boundaryField()[patchI]
                    );

                const processorPointPatch& procPatch =
                    refCast<const processorPointPatch>
                    (
                        pFld.patch()
                    );

                IPstream fromNeighbProc
                (
                    procPatch.neighbProcNo(),
                    pFld.size()*sizeof(Type)
                );

                Field<Type> nbrFld(fromNeighbProc);

                // Combing neighbouring values with my values.
                Field<Type> myFld(pFld.patchInternalField());

                cop(myFld, nbrFld);

                // Get the addressing
                const labelList& mp = procPatch.meshPoints();

                forAll(pFld, i)
                {
                    fld[mp[i]] = myFld[i];
                }
            }
        }
    }


    typedef GlobalProcessorPointPatchField
        <pointPatchField, pointPatch, globalProcessorPointPatch, Type>
        GlobalProcessorField;

    // Shared points.
    if (Pstream::parRun() && mesh.parallelData().parallel())
    {
        forAll(fld.boundaryField(), patchI)
        {
            if
            (
                isA<GlobalProcessorField>
                (
                    fld.boundaryField()[patchI]
                )
            )
            {
                const GlobalProcessorField& pFld =
                    refCast<const GlobalProcessorField>
                    (
                        fld.boundaryField()[patchI]
                    );

                const globalProcessorPointPatch& procPatch =
                    refCast<const globalProcessorPointPatch>
                    (
                        pFld.patch()
                    );



                Field<Type> gpf(procPatch.globalPointSize(), zero);

                // Get the addressing
                const labelList& mp = procPatch.meshPoints();
                const labelList& addr = procPatch.sharedPointAddr();

                forAll (addr, i)
                {
                    gpf[addr[i]] = fld[mp[i]];
                }

                combineReduce(gpf, cop);

                forAll (addr, i)
                {
                    fld[mp[i]] = gpf[addr[i]];
                }
            }
        }
    }    
}


// ************************************************************************* //
