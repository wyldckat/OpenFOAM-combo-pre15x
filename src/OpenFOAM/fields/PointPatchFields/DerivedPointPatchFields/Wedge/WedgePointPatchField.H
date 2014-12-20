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
    WedgePointPatchField

Description
    Wedge front and back plane patch field

SourceFiles
    WedgePointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef WedgePointPatchField_H
#define WedgePointPatchField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class WedgePointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class PointPatch,
    class WedgePointPatch,
    class Type
>
class WedgePointPatchField
:
    public PatchField<Type>
{

public:

    //- Runtime type information
    TypeName(WedgePointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        WedgePointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        WedgePointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        WedgePointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        WedgePointPatchField
        (
            const WedgePointPatchField
                <PatchField, PointPatch, WedgePointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new WedgePointPatchField
                <PatchField, PointPatch, WedgePointPatch, Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        WedgePointPatchField
        (
            const WedgePointPatchField
                <PatchField, PointPatch, WedgePointPatch, Type>&,
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
                new WedgePointPatchField
                <PatchField,PointPatch, WedgePointPatch, Type>
                (
                    *this, iF
                )
            );
        }


    // Member functions

        //- Update the patch field
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "WedgePointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
