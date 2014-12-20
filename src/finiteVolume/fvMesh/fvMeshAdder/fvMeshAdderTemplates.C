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


#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMeshAdder::map
(
    const Field<Type>& oldFld,
    const labelList& oldToNew,
    Field<Type>& fld
)
{
    forAll(oldToNew, cellI)
    {
        label newCellI = oldToNew[cellI];

        if (newCellI != -1)
        {
            fld[newCellI] = oldFld[cellI];
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::mapPatchField
(
    const Field<Type>& oldFld,
    const label oldStart,
    const labelList& oldToNew,
    fvPatchField<Type>& fld
)
{
    // Start of new patch.
    const label newStart = fld.patch().patch().start();
    const label newSize = fld.patch().patch().size();

    forAll(oldFld, i)
    {
        // From old patch face to mesh face to new mesh face.
        label newFaceI = oldToNew[i + oldStart];

        if (newFaceI >= newStart+newSize)
        {
            FatalErrorIn("mapPatchField")
                << "oldPatchface:" << i << " newMeshFace:" << newFaceI
                << " oldSize:" << oldFld.size()
                << " oldStart:" << oldStart
                << " newStart:" << newStart
                << " newSize:" << newSize
                << abort(FatalError);
        }

        if (newFaceI >= newStart)
        {
            // So  -face is still there  -has not become internal
            fld[newFaceI - newStart] = oldFld[i];
        }
    }
}


template<class Type>
void Foam::fvMeshAdder::MapVolField
(
    const labelList& oldPatchStarts,
    GeometricField<Type, fvPatchField, volMesh>& fld,
    const GeometricField<Type, fvPatchField, volMesh>& fldToAdd
) const
{
    const fvMesh& mesh = fld.mesh();

    const mapAddedPolyMesh& meshMap = polyMeshAdder::map();

    // Internal field
    // ~~~~~~~~~~~~~~

    // Store old internal field
    Field<Type> oldInternalField(fld.internalField());

    // Modify internal field
    Field<Type>& intFld = fld.internalField();

    intFld.setSize(meshMap.cellMap().size());

    map(oldInternalField, meshMap.oldCellMap(), intFld);
    map(fldToAdd.internalField(), meshMap.addedCellMap(), intFld);


    // Patch fields
    // ~~~~~~~~~~~~

    // Problem with the patch fields is that we need to reset both
    // the fvPatch reference and the internal field reference and
    // the .clone function only allows resetting the internal
    // field so we have to recreate completely new patch fields.

    const labelList& oldPatchMap = meshMap.oldPatchMap();
    const labelList& addedPatchMap = meshMap.addedPatchMap();

    // Create new patchFields.
    PtrList<fvPatchField<Type> > patchFields
    (
        mesh.boundary().size()
    );

    // Add old mesh patches
    forAll(oldPatchMap, patchI)
    {
        label newPatchI = oldPatchMap[patchI];

        if (newPatchI != -1)
        {
            // Create new patchField with same type as existing one.
            patchFields.set(newPatchI) =
                fvPatchField<Type>::New
                (
                    fld.boundaryField()[patchI].type(),
                    mesh.boundary()[newPatchI],
                    fld.internalField()
                ).ptr();

            // Map values.
            mapPatchField
            (
                fld.boundaryField()[patchI],
                oldPatchStarts[patchI],
                meshMap.oldFaceMap(),
                patchFields[newPatchI]
            );
        }
    }

    // Add addedMesh patches
    forAll(addedPatchMap, patchI)
    {
        label newPatchI = addedPatchMap[patchI];

        if (newPatchI != -1)
        {
            if (!patchFields(newPatchI))
            {
                patchFields.set(newPatchI) =
                    fvPatchField<Type>::New
                    (
                        fldToAdd.boundaryField()[patchI].type(),
                        mesh.boundary()[newPatchI],
                        fld.internalField()
                    ).ptr();
            }

            // Map values.
            mapPatchField
            (
                fldToAdd.boundaryField()[patchI],
                fldToAdd.boundaryField()[patchI].patch().patch().start(),
                meshMap.addedFaceMap(),
                patchFields[newPatchI]
            );
        }
    }

    // Patch fields should now have correct
    // -patch -internalField reference.

    // Modify the boundaryField. Keep the reference to the
    // boundary mesh intact.
    fld.boundaryField().clear();
    fld.boundaryField().setSize(patchFields.size());
    forAll(patchFields, i)
    {
        fld.boundaryField().hook(patchFields[i].clone());
    }
}


template<class Type>
void Foam::fvMeshAdder::MapVolFields
(
    const labelList& oldPatchStarts,
    const fvMesh& mesh,
    const fvMesh& meshToAdd
) const
{
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mesh.objectRegistry::lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >
        ()
    );

    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fieldsToAdd
    (
        meshToAdd.objectRegistry::lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >
        ()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.  

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        const_cast<GeometricField<Type, fvPatchField, volMesh>*>(fieldIter())
            ->storeOldTimes();
    }


    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, fvPatchField, volMesh>& fld =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (
                *fieldIter()
            );

        if (fieldsToAdd.found(fld.name()))
        {
            Pout<< "Mapping field " << fld.name() << endl;

            const GeometricField<Type, fvPatchField, volMesh>& fldToAdd =
                *fieldsToAdd[fld.name()];

            MapVolField<Type>(oldPatchStarts, fld, fldToAdd);
        }
        else
        {
            Pout<< "Not mapping field " << fld.name()
                << " since not present on mesh to add"
                << endl;
        }
    }
}


//template<class Type>
//Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
//Foam::fvMeshAdder::add
//(
//    const fvMesh& mergedMesh,
//    const GeometricField<Type, fvPatchField, volMesh>& fld0,
//    const GeometricField<Type, fvPatchField, volMesh>& fld1
//) const
//{
//    // Create and map the internal-field values. Note no need to initialize
//    // internal field since all cells guaranteed to originate from one of either
//    // meshes.
//    Field<Type> internalField(mergedMesh.nCells());
//
//    map(fld0, meshMap.oldCellMap(), internalField);
//    map(fld1, meshMap.addedCellMap(), internalField);
//
//    // Create and map the patch field values
//    PtrList<fvPatchField<Type> > patchFields(mergedMesh.boundary().size());
//
//    // Add mesh0 patches
//    forAll(meshMap.oldPatchMap(), patchI)
//    {
//        label newPatchI = meshMap.oldPatchMap()[patchI];
//
//        if (newPatchI != -1)
//        {
//            // Create new patchField with same type as fld0.
//
//            Pout<< "    Copying old mesh field at "
//                << fld0.boundaryField()[patchI].patch().name()
//                << "  type " << fld0.boundaryField()[patchI].type()
//                << " to position " << newPatchI << endl;
//
//            patchFields.set(newPatchI) = fvPatchField<Type>::New
//            (
//                fld0.boundaryField()[patchI].type(),
//                mergedMesh.boundary()[newPatchI],
//                internalField
//            ).ptr();
//
//            // Copy values.
//            mapPatch
//            (
//                fld0.boundaryField()[patchI],
//                fld0.boundaryField()[patchI].patch().patch().start(),
//                meshMap.oldFaceMap(),
//                patchFields[newPatchI]
//            );
//        }
//    }
//
//    // Add mesh1 patches
//    forAll(meshMap.addedPatchMap(), patchI)
//    {
//        label newPatchI = meshMap.addedPatchMap()[patchI];
//
//        if (newPatchI != -1)
//        {
//            if (!patchFields(newPatchI))
//            {
//                Pout<< "    Copying added mesh field at "
//                    << fld0.boundaryField()[patchI].patch().name()
//                    << "  type " << fld0.boundaryField()[patchI].type()
//                    << " to position " << newPatchI << endl;
//
//                patchFields.set(newPatchI) = fvPatchField<Type>::New
//                (
//                    fld1.boundaryField()[patchI].type(),
//                    mergedMesh.boundary()[newPatchI],
//                    internalField
//                ).ptr();
//            }
//
//            // Copy values.
//            mapPatch
//            (
//                fld1.boundaryField()[patchI],
//                fld1.boundaryField()[patchI].patch().patch().start(),
//                meshMap.addedFaceMap(),
//                patchFields[newPatchI]
//            );
//        }
//    }
//
//    // Create the complete field from the pieces
//    tmp<GeometricField<Type, fvPatchField, volMesh> > tresF
//    (
//        new GeometricField<Type, fvPatchField, volMesh>
//        (
//            IOobject
//            (
//                "merged"+fld0.name(),
//                mergedMesh.time().timeName(),
//                mergedMesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mergedMesh,
//            fld0.dimensions(),
//            internalField,
//            patchFields
//        )
//    );
//
//    return tresF;
//}


//template<class Type>
//void Foam::fvMeshAdder::readAndAddTypedVolFields
//(
//    const fvMesh& mesh,
//    const IOobjectList& objects,
//    const fvMesh& meshToAdd,
//    const IOobjectList& objectsToAdd,
//    const fvMesh& mergedMesh
//) const
//{
//    // Get objects of Type.
//    wordList fieldNames
//    (
//        objects.names
//        (
//            GeometricField<Type, fvPatchField, volMesh>::typeName
//        )
//    );
//
//
//    // Merge _0 fields in iteration 2, all others in iteration 1.
//    // Bit of a hack. GeometricField will try to load old time step field
//    // ('readOldTimeIfPresent') so make sure it isn't mapped before the
//    // this-timestep field since the _0 field will have the new mesh size.
//
//    // We enforce this by first doing all fields ending not in _0, then those
//    // ending in _0, then those ending in _0_0 etc. until all have been done.
//
//    word postFix = "";
//
//    label nFiles = 0;
//
//    do
//    {
//        // Merge all fields ending in postFix.
//
//        forAll(fieldNames, i)
//        {
//            const word& fieldName = fieldNames[i];
//
//            string::size_type pos = fieldName.rfind(postFix);
//
//            if (pos != string::npos && pos == fieldName.size() - postFix.size())
//            {
//                // Check whether another level of _0
//                word baseName = fieldName.substr(0, pos);
//
//                string::size_type pos2 = baseName.rfind("_0");
//
//                if (pos2 == string::npos || pos2 != baseName.size() - 2)
//                {
//                    nFiles++;
//
//                    Info<< "Merging field " << fieldName << " ..." << endl;
//
//                    if (objectsToAdd.found(fieldName))
//                    {
//                        // Load both fields
//
//                        GeometricField<Type, fvPatchField, volMesh> fld
//                        (
//                            *objects[fieldName],
//                            mesh
//                        );
//
//                        GeometricField<Type, fvPatchField, volMesh> fldToAdd
//                        (
//                            *objectsToAdd[fieldName],
//                            meshToAdd
//                        );
//
//                        // Merge the two fields onto the new mergedMesh
//                        tmp<GeometricField<Type, fvPatchField, volMesh> > tfld
//                        (
//                            add(mergedMesh, fld, fldToAdd)
//                        );
//
//                        tfld().rename(fieldName);
//                        tfld().write();
//                    }
//                    else
//                    {
//                        Warning<< "Cannot find field " << fieldName
//                            << " on mesh to add" << endl;
//                    }
//                }
//            }
//        }
//
//        // New iteration going through older timelevel fields
//        postFix += "_0";
//
//    } while (nFiles < fieldNames.size());
//
//}


// ************************************************************************* //
