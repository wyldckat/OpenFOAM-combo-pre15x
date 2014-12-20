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

Class
    fvMeshDistribute

\*----------------------------------------------------------------------------*/

#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class T, class Mesh>
void Foam::fvMeshDistribute::printFieldInfo(const fvMesh& mesh)
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        mesh.objectRegistry::lookupClass<fldType>()
    );

    for
    (
        typename HashTable<const fldType*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const fldType& fld = *iter();

        Pout<< "Field:" << iter.key() << " internalsize:" << fld.size()
            //<< " value:" << fld
            << endl;

        forAll(fld.boundaryField(), patchI)
        {
            Pout<< "    " << patchI
                << ' ' << fld.boundaryField()[patchI].patch().name()
                << ' ' << fld.boundaryField()[patchI].type()
                << ' ' << fld.boundaryField()[patchI].size()
                //<< ' ' << fld.boundaryField()[patchI]
                << endl;
        }
    }
}


// Add patch field
template <class T, class Mesh>
void Foam::fvMeshDistribute::addPatchFields(const word& patchFieldType)
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    for
    (
        typename HashTable<const fldType*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const fldType& fld = *iter();

        typename fldType::GeometricBoundaryField& bfld =
            const_cast<typename fldType::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        label sz = bfld.size();
        bfld.setSize(sz+1);
        bfld.set
        (
            sz,
            fvPatchField<T>::New
            (
                patchFieldType,
                mesh_.boundary()[sz],
                fld.internalField()
            )
        );
    }
}


// Delete trailing patch fields
template<class T, class Mesh>
void Foam::fvMeshDistribute::deleteTrailingPatchFields()
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    for
    (
        typename HashTable<const fldType*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const fldType& fld = *iter();

        typename fldType::GeometricBoundaryField& bfld =
            const_cast<typename fldType::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        label sz = bfld.size();

        // Shrink patchFields
        bfld.setSize(sz-1);
    }
}


// Save whole boundary field
template <class T, class Mesh>
void Foam::fvMeshDistribute::saveBoundaryFields
(
    PtrList<FieldField<fvPatchField, T> >& bflds
) const
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    bflds.setSize(flds.size());

    label i = 0;

    for
    (
        typename HashTable<const fldType*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const fldType& fld = *iter();

        bflds.set(i, fld.boundaryField().clone().ptr());

        i++;
    }
}


// Map boundary field
template <class T, class Mesh>
void Foam::fvMeshDistribute::mapBoundaryFields
(
    const mapPolyMesh& map,
    const PtrList<FieldField<fvPatchField, T> >& oldBflds
)
{
    const labelList& oldPatchStarts = map.oldPatchStarts();
    const labelList& faceMap = map.faceMap();

    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    if (flds.size() != oldBflds.size())
    {
        FatalErrorIn("fvMeshDistribute::mapBoundaryFields") << "problem"
            << abort(FatalError);
    }

    label fieldI = 0;

    for
    (
        typename HashTable<const fldType*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const fldType& fld = *iter();
        typename fldType::GeometricBoundaryField& bfld =
            const_cast<typename fldType::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );


        const FieldField<fvPatchField, T>& oldBfld = oldBflds[fieldI++];

        // Pull from old boundary field into bfld.

        forAll(bfld, patchI)
        {
            fvPatchField<T>& patchFld = bfld[patchI];
            label faceI = patchFld.patch().patch().start();

            forAll(patchFld, i)
            {
                label oldFaceI = faceMap[faceI++];

                // Find patch and local patch face oldFaceI was in.
                forAll(oldPatchStarts, oldPatchI)
                {
                    label oldLocalI = oldFaceI - oldPatchStarts[oldPatchI];

                    if (oldLocalI >= 0 && oldLocalI < oldBfld[oldPatchI].size())
                    {
                        patchFld[i] = oldBfld[oldPatchI][oldLocalI];
                    }
                }
            }
        }
    }
}


// Init patch fields of certain type
template <class T, class Mesh>
void Foam::fvMeshDistribute::initPatchFields
(
    const word& patchFieldType,
    const T& initVal
)
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    HashTable<const fldType*> flds
    (
        mesh_.objectRegistry::lookupClass<fldType>()
    );

    for
    (
        typename HashTable<const fldType*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const fldType& fld = *iter();

        typename fldType::GeometricBoundaryField& bfld =
            const_cast<typename fldType::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        forAll(bfld, patchI)
        {
            if (bfld[patchI].type() == patchFieldType)
            {
                bfld[patchI] == initVal;
            }
        }
    }
}


// Send fields. Note order supplied so we can receive in exactly the same order.
template <class T, class Mesh>
void Foam::fvMeshDistribute::sendFields
(
    const label domain,
    const wordList& fieldNames,
    const fvMeshSubset& subsetter
)
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    forAll(fieldNames, i)
    {
        //Pout<< "Subsetting field " << fieldNames[i]
        //    << " for domain:" << domain
        //    << endl;

        // Send all fieldNames. This has to be exactly the same set as is
        // being received!
        const fldType& fld =
            subsetter.baseMesh().lookupObject<fldType>(fieldNames[i]);

        tmp<fldType> tsubfld = subsetter.interpolate(fld);

        // Send
        OPstream toNbr(domain);
        toNbr << tsubfld();
    }
}


// Opposite of sendFields
template<class T, class Mesh>
void Foam::fvMeshDistribute::receiveFields
(
    const label domain,
    const wordList& fieldNames,
    fvMesh& mesh,
    PtrList<GeometricField<T, fvPatchField, Mesh> >& fields
)
{
    typedef GeometricField<T, fvPatchField, Mesh> fldType;

    fields.setSize(fieldNames.size());

    forAll(fieldNames, i)
    {
        //Pout<< "Receiving field " << fieldNames[i]
        //    << " from domain:" << domain
        //    << endl;

        IPstream fromNbr(domain);

        fields.set
        (
            i,
            new fldType
            (
                IOobject
                (
                    fieldNames[i],
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                fromNbr
            )
        );
    }
}


// ************************************************************************* //
