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

#include "geometricTetPointFieldReconstructor.H"
#include "ptrList.H"
#include "tetPolyPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
geometricTetPointFieldReconstructor::reconstructTetPointField
(
    const IOobject& fieldIoObject
)
{
    // Read the field for all the processors
    ptrList<GeometricField<Type, tetPolyPatchField, tetPointMesh> > procFields
    (
        procMeshes_.size()
    );

    forAll (procMeshes_, procI)
    {
        procFields.hook
        (
            new GeometricField<Type, tetPolyPatchField, tetPointMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[procI]().time().timeName(),
                    procMeshes_[procI](),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[procI]
            )
        );
    }


    // Create the internalField
    Field<Type> internalField(mesh_.nPoints());

    // Create the patch fields
    ptrList<tetPolyPatchField<Type> > patchFields(mesh_.boundary().size());


    forAll (procMeshes_, procI)
    {
        const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
            procField = procFields[procI];

        // Get processor-to-global addressing for use in rmap
        labelList procToGlobalAddr = procAddressing(procI);

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.internalField(),
            procToGlobalAddr
        );

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[procI], patchI)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[procI][patchI];

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                if (!patchFields(curBPatch))
                {
                    patchFields.set(curBPatch) =
                        tetPolyPatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[curBPatch],
                            internalField,
                            tetPolyPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        ).ptr();
                }

                // If the field stores values, do the rmap
                if (patchFields[curBPatch].storesFieldData())
                {
                    patchFields[curBPatch].rmap
                    (
                        procField.boundaryField()[patchI],
                        procPatchAddressing
                        (
                            procToGlobalAddr,
                            procI,
                            patchI
                        )
                    );
                }
            }
        }
    }

    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
    (
        new GeometricField<Type, tetPolyPatchField, tetPointMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_().time().timeName(),
                mesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


// Reconstruct and write all volume fields
template<class Type>
void geometricTetPointFieldReconstructor::reconstructTetPointFields
(
    const IOobjectList& objects
)
{
    word fieldClassName
    (
        GeometricField<Type, tetPolyPatchField, tetPointMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        for
        (
            IOobjectList::iterator fieldIter = fields.begin();
            fieldIter != fields.end();
            ++fieldIter
        )
        {
            Info<< "        " << fieldIter()->name() << endl;

            reconstructTetPointField<Type>(*fieldIter())().write();
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
