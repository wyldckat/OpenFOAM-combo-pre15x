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

\*---------------------------------------------------------------------------*/

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
    const GeometricField<Type, fvPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                "0",
                ssf.dimensions()/dimVol,
                pTraits<Type>::zero
            ),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();


    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        vf[owner[faceI]] += ssf[faceI];
        vf[neighbour[faceI]] -= ssf[faceI];
    }

    forAll(mesh.boundary(), patchI)
    {
        const labelList::subList pFaceCells =
            mesh.boundary()[patchI].faceCells();

        const fvPatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            vf[pFaceCells[faceI]] += pssf[faceI];
        }
    }

    vf.internalField() /= mesh.V();
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, surfaceMesh> >& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        fvc::surfaceIntegrate(tssf())
    );
    tssf.clear();
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
surfaceSum
(
    const GeometricField<Type, fvPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", ssf.dimensions(), pTraits<Type>::zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        vf[owner[faceI]] += ssf[faceI];
        vf[neighbour[faceI]] += ssf[faceI];
    }

    forAll(mesh.boundary(), patchI)
    {
        const labelList::subList pFaceCells =
            mesh.boundary()[patchI].faceCells();

        const fvPatchField<Type>& pssf = ssf.boundaryField()[patchI];

        forAll(mesh.boundary()[patchI], faceI)
        {
            vf[pFaceCells[faceI]] += pssf[faceI];
        }
    }

    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > surfaceSum
(
    const tmp<GeometricField<Type, fvPatchField, surfaceMesh> >& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
