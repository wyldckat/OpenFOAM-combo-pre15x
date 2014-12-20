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
    CyclicPointPatchField

Description
    Cyclic front and back plane patch field

SourceFiles
    CyclicPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef CyclicPointPatchField_H
#define CyclicPointPatchField_H

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CyclicPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class PointPatch,
    class CyclicPointPatch,
    class Type
>
class CyclicPointPatchField
:
    public CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        Type
    >
{
    // Private data

        //- Local reference cast into the cyclic patch
        const CyclicPointPatch& cyclicPatch_;


public:

    //- Runtime type information
    TypeName(CyclicPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        CyclicPointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        CyclicPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        CyclicPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        CyclicPointPatchField
        (
            const CyclicPointPatchField
                <PatchField, PointPatch, CyclicPointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new CyclicPointPatchField
                <PatchField, PointPatch, CyclicPointPatch, Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        CyclicPointPatchField
        (
            const CyclicPointPatchField
                <PatchField, PointPatch, CyclicPointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const Field<Type>& iF
        ) const
        {
            return autoPtr<PatchField<Type> >
            (
                new CyclicPointPatchField
                <PatchField,PointPatch, CyclicPointPatch, Type>
                (
                    *this, iF
                )
            );
        }


    // Member functions

        //- Cyclic coupled interface functions

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                return !(cyclicPatch_.parallel() || pTraits<Type>::rank == 0);
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPatch_.reverseT();
            }


        //- Evaluate the patch field
        virtual void evaluate()
        {}

        //- Complete swap of patch point values and add to local values
        virtual void swapAdd(Field<Type>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "CyclicPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
