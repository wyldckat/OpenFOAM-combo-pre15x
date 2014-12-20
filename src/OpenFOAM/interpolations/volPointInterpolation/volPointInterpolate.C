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

Class
    volPointInterpolation

Description

\*---------------------------------------------------------------------------*/

#include "volPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "primitiveMesh.H"
#include "emptyFvPatch.H"
#include "parallelInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "volPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const FieldField<Field, scalar>& weights = pointWeights();
    const labelListList& pointCells = vf.mesh().pointCells();

    // Multiply volField by weighting factor matrix to create pointField
    forAll(pointCells, pointi)
    {
        pf[pointi] = pTraits<Type>::zero;

        forAll(pointCells[pointi], pointCelli)
        {
            pf[pointi] +=
                weights[pointi][pointCelli]*vf[pointCells[pointi][pointCelli]];
        }
    }


    // Interpolate patch values: over-ride the internal values for the points
    // on the patch with the interpolated point values from the faces
    const fvBoundaryMesh& bm = vMesh().boundary();

    const PtrList<primitivePatchInterpolation>& pi = patchInterpolators();
    forAll (bm, patchI)
    {
        // If the patch is empty, skip it
        // If the patch is coupled, and there are no cyclic parallels, skip it
        if
        (
            bm[patchI].type() != emptyFvPatch::typeName
//         && !bm[patchI].coupled()
        && !(
                bm[patchI].coupled()
             && Pstream::parRun()
             && !vMesh().parallelData().cyclicParallel()
            )
        )
        {
            pf.boundaryField()[patchI].setInInternalField
            (
                pf.internalField(),
                pi[patchI].faceToPointInterpolate
                (
                    vf.boundaryField()[patchI]
                )()
            );
        }
    }

    // Do edge correction

    const labelList& ptc = boundaryPoints();

    const FieldField<Field, scalar>& pbw = pointBoundaryWeights();

    const labelListList& PointFaces = vMesh().pointFaces();

    forAll (ptc, pointI)
    {
        const label curPoint = ptc[pointI];

        const labelList& curFaces = PointFaces[curPoint];

        label fI = 0;

        // Reset the boundary value before accumulation
        pf[curPoint] = pTraits<Type>::zero;

        // Go through all the faces
        forAll (curFaces, faceI)
        {
            if (!vMesh().isInternalFace(curFaces[faceI]))
            {
                // This is a boundary face.  If not in the empty patch
                // or coupled calculate the extrapolation vector
                label patchID =
                    vMesh().boundaryMesh().whichPatch(curFaces[faceI]);

                if
                (
                    vMesh().boundary()[patchID].type()
                 != emptyFvPatch::typeName
//                 && !vMesh().boundary()[patchID].coupled()
                && !(
                        vMesh().boundary()[patchID].coupled()
                     && Pstream::parRun()
                     && !vMesh().parallelData().cyclicParallel()
                    )
                )
                {
                    label faceInPatchID =
                        bm[patchID].patch().whichFace(curFaces[faceI]);

                    pf[curPoint] +=
                        pbw[pointI][fI]*
                        vf.boundaryField()[patchID][faceInPatchID];

                    fI++;
                }
            }
        }
    }

    // Update coupled boundaries
    // Work-around for cyclic parallels.  
    if (Pstream::parRun() && !vMesh().parallelData().cyclicParallel())
    {
        forAll (pf.boundaryField(), patchI)
        {
            if (pf.boundaryField()[patchI].coupled())
            {
                pf.boundaryField()[patchI].initAddField();
            }
        }

        forAll (pf.boundaryField(), patchI)
        {
            if (pf.boundaryField()[patchI].coupled())
            {
                pf.boundaryField()[patchI].addField(pf.internalField());
            }
        }
    }

    pf.correctBoundaryConditions();

    if (debug)
    {
        Info<< "volPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "finished interpolating field from cells to points"
            << endl;
    }
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "volPointInterpolate(" + vf.name() + ')',
                vf.instance(),
                vf.db()
            ),
            pointMesh_,
            vf.dimensions()
        )
    );

    // Perform interpolation
    interpolate(vf, tpf());

    return tpf;
}


template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
volPointInterpolation::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpf =
        interpolate(tvf());
    tvf.clear();
    return tpf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
