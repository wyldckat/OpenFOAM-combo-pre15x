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

\*---------------------------------------------------------------------------*/

#include "meshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > meshSubset::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Get reference to the subset mesh
    const fvMesh& sMesh = subMesh();

    // Create and map the internal-field values
    Field<Type> internalField(vf.internalField(), cellMap());

    // Create and map the patch field values
    const labelList& pm = patchMap();

    ptrList<fvPatchField<Type> > patchFields(pm.size());

    forAll (patchFields, patchI)
    {
        // Hook the first one by hand as it corresponds to the
        // internal field.  Additional interpolation can be put here
        // as necessary.  
        if (pm[patchI] == -1)
        {
            patchFields.hook
            (
                new emptyFvPatchField<Type>
                (
                    sMesh.boundary()[patchI],
                    internalField
                )
            );
        }
        else
        {
            patchFields.hook
            (
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[pm[patchI]],
                    sMesh.boundary()[patchI],
                    internalField,
                    patchFieldSubset(*this, patchI)
                )
            );
        }
    }


    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > tresF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "subset"+vf.name(),
                sMesh.time().timeName(),
                sMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sMesh,
            vf.dimensions(),
            internalField,
            patchFields
        )
    );

    return tresF;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
