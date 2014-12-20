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

#include "ProcessorPointPatchField.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Raw field sending and receiving
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::sendPointField
(
    const tmp<Field<Type2> >& tf
) const
{
    OPstream::write
    (
        procPatch_.neighbProcNo(),
        reinterpret_cast<const char*>(tf().begin()),
        tf().size()*sizeof(Type2)
    );
    tf.clear();
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::sendEdgeField
(
    const tmp<Field<Type2> >& tf
) const
{
    OPstream::write
    (
        procPatch_.neighbProcNo(),
        reinterpret_cast<const char*>(tf().begin()),
        procPatch_.localEdgeIndices().size()*sizeof(Type2)
    );
    tf.clear();
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
tmp<Field<Type2> >
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
receivePointField() const
{
    tmp<Field<Type2> > tf(new Field<Type2>(this->size()));
    Field<Type2>& f = tf();

    IPstream::read
    (
        procPatch_.neighbProcNo(),
        reinterpret_cast<char*>(f.begin()),
        this->size()*sizeof(Type2)
    );

    return tf;
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
tmp<Field<Type2> >
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
receiveEdgeField() const
{
    tmp<Field<Type2> > tf
    (
        new Field<Type2>(procPatch_.localEdgeIndices().size())
    );
    Field<Type2>& f = tf();

    IPstream::read
    (
        procPatch_.neighbProcNo(),
        reinterpret_cast<char*>(f.begin()),
        procPatch_.localEdgeIndices().size()*sizeof(Type2)
    );

    return tf;
}


// Initialise diagonal/source update.
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddFieldTempl
(
    const Field<Type2>& pField
) const
{
    sendPointField(patchInternalField(pField));
}


// Add the diagonal/source to the internal field.
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::
addFieldTempl
(
    Field<Type2>& pField
) const
{
    // Get the neighbour side values
    tmp<Field<Type2> > tpNeighbour = receivePointField<Type2>();
    addToInternalField(pField, tpNeighbour());
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p))
{
    if (!isType<ProcessorPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "ProcessorPointPatchField"
            "<PatchField, PointPatch, ProcessorPointPatch, Type>::"
            "ProcessorPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not processor type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
        <PatchField, PointPatch, ProcessorPointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch()))
{
    if (!isType<ProcessorPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "ProcessorPointPatchField"
            "<PatchField, PointPatch, ProcessorPointPatch, Type>::"
            "ProcessorPointPatchField\n"
            "(\n"
            "    const ProcessorPointPatchField"
                 "<PatchField, PointPatch, ProcessorPointPatch, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const PointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
        <PatchField, PointPatch, ProcessorPointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(ptf, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
~ProcessorPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Initialise field transfer
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddField() const
{
    initAddFieldTempl(this->internalField());
}


// Add field
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
addField(Field<Type>& f) const
{
    addFieldTempl(f);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
