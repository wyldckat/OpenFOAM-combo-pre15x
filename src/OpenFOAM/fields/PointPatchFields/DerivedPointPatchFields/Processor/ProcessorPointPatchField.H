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
    ProcessorPointPatchField

Description

SourceFiles
    ProcessorPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessorPointPatchField_H
#define ProcessorPointPatchField_H

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ProcessorPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
class ProcessorPointPatchField
:
    public CoupledPointPatchField<PatchField, PointPatch, Type>
{
    // Private data

        //- Local reference to processor patch
        const ProcessorPointPatch& procPatch_;


    // Private member functions

        //- Raw field send function
        template<class Type2>
        void sendPointField(const tmp<Field<Type2> >&) const;

        template<class Type2>
        void sendEdgeField(const tmp<Field<Type2> >&) const;


        //- Raw field receive function
        template<class Type2>
        tmp<Field<Type2> > receivePointField() const;

        template<class Type2>
        tmp<Field<Type2> > receiveEdgeField() const;

        //- Templated add diagonal/source functions
        template<class Type2>
        void initAddFieldTempl(const Field<Type2>&) const;

        template<class Type2>
        void addFieldTempl(Field<Type2>&) const;



public:

    //- Runtime type information
    TypeName(ProcessorPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        ProcessorPointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        ProcessorPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        ProcessorPointPatchField
        (
            const ProcessorPointPatchField
                <PatchField, PointPatch, ProcessorPointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new ProcessorPointPatchField
                <PatchField, PointPatch, ProcessorPointPatch, Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        ProcessorPointPatchField
        (
            const ProcessorPointPatchField
                <PatchField, PointPatch, ProcessorPointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone(const Field<Type>& iF) const
        {
            return autoPtr<PatchField<Type> >
            (
                new ProcessorPointPatchField
                <PatchField, PointPatch, ProcessorPointPatch, Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Destructor

        ~ProcessorPointPatchField();


    // Member functions

        // Access

            //- Evaluate the patch field
            virtual void evaluate()
            {}

            //- Init add field to internal field
            virtual void initAddField() const;

            //- Add field to internal field
            virtual void addField(Field<Type>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "ProcessorPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
