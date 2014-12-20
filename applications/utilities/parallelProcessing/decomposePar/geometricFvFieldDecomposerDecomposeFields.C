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

#include "geometricFvFieldDecomposer.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "processorFvPatchField.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoField>
void geometricFvFieldDecomposer::readFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    ptrList<GeoField>& fields
)
{
    // Search list of objects for volScalarFields
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Remove the cellDist field
    IOobjectList::iterator celDistIter = fieldObjects.find("cellDist");
    if (celDistIter != fieldObjects.end())
    {
        fieldObjects.erase(celDistIter);
    }

    // Construct the vol scalar fields
    fields.setSize(fieldObjects.size());

    for
    (
        IOobjectList::iterator iter = fieldObjects.begin();
        iter != fieldObjects.end();
        ++iter
    )
    {
        fields.hook
        (
            new GeoField
            (
                *iter(),
                mesh
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
geometricFvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.internalField(), cellAddressing_);

    // Create and map the patch field values
    ptrList<fvPatchField<Type> > patchFields(boundaryAddressing_.size());

    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFields.hook
            (
                fvPatchField<Type>::New
                (
                    field.boundaryField()
                        [boundaryAddressing_[patchI]],
                    processorMesh_.boundary()[patchI],
                    internalField,
                    fvPatchFieldDecomposer
                    (
                        processorMesh_.boundary()[patchI].patchSlice
                        (
                            faceAddressing_
                        ),
                        field.mesh().boundaryMesh()
                        [boundaryAddressing_[patchI]].start()
                    )
                )
            );
        }
        else
        {
            patchFields.hook
            (
                new processorFvPatchField<Type>
                (
                    processorMesh_.boundary()[patchI],
                    internalField,
                    Field<Type>
                    (
                        field.internalField(),
                        processorFvPatchVolFieldDecomposer
                        (
                            field.mesh(),
                            processorMesh_.boundary()[patchI].patchSlice
                            (
                                faceAddressing_
                            )
                        )
                    )
                )
            );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                field.name(),
                processorMesh_.time().timeName(),
                processorMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            processorMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
geometricFvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvPatchField, surfaceMesh>& field
) const
{
    labelList mapAddr
    (
        labelList::subList
        (
            faceAddressing_,
            processorMesh_.nInternalFaces()
        )
    );
    forAll (mapAddr, i)
    {
        mapAddr[i] -= 1;
    }

    // Create and map the internal field values
    Field<Type> internalField
    (
        field.internalField(),
        mapAddr
    );

    // Problem with addressing when a processor patch picks up both internal
    // faces and faces from cyclic boundaries. This is a bit of a hack, but
    // I cannot find a better solution without making the internal storage
    // mechanism for surfaceFields correspond to the one of faces in polyMesh
    // (i.e. using slices)
    Field<Type> allFaceField(field.mesh().nFaces());

    forAll (field.internalField(), i)
    {
        allFaceField[i] = field.internalField()[i];
    }

    forAll (field.boundaryField(), patchI)
    {
        const Field<Type> & p = field.boundaryField()[patchI];

        const label patchStart = field.mesh().boundaryMesh()[patchI].start();

        forAll (p, i)
        {
            allFaceField[patchStart + i] = p[i];
        }
    }

    // Create and map the patch field values
    ptrList<fvPatchField<Type> > patchFields(boundaryAddressing_.size());

    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFields.hook
            (
                fvPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchI]],
                    processorMesh_.boundary()[patchI],
                    internalField,
                    fvPatchFieldDecomposer
                    (
                        processorMesh_.boundary()[patchI].patchSlice
                        (
                            faceAddressing_
                        ),
                        field.mesh().boundaryMesh()
                            [boundaryAddressing_[patchI]].start()
                    )
                )
            );
        }
        else
        {
            patchFields.hook
            (
                new processorFvPatchField<Type>
                (
                    processorMesh_.boundary()[patchI],
                    internalField,
                    Field<Type>
                    (
                        allFaceField,
                        processorFvPatchSurfaceFieldDecomposer
                        (
                            (const unallocLabelList&)
                            processorMesh_.boundary()[patchI].patchSlice
                            (
                                faceAddressing_
                            )
                        )
                    )
                )
            );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
    (
        new GeometricField<Type, fvPatchField, surfaceMesh>
        (
            IOobject
            (
                field.name(),
                processorMesh_.time().timeName(),
                processorMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            processorMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class GeoField>
void geometricFvFieldDecomposer::decomposeFields
(
    const ptrList<GeoField>& fields
) const
{
    forAll (fields, fieldI)
    {
        decomposeField(fields[fieldI])().write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
