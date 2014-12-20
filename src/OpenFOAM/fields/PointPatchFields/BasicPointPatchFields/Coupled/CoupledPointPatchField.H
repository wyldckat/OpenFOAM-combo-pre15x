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
    CoupledPointPatchField

Description

SourceFiles
    CoupledPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef CoupledPointPatchField_H
#define CoupledPointPatchField_H

#include "PointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class CoupledPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class PointPatch,
    class CoupledPointPatch,
    class Type
>
class CoupledPointPatchField
:
    public PatchField<Type>
{

public:

    //- Runtime type information
    TypeName(CoupledPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        CoupledPointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        CoupledPointPatchField
        (
            const CoupledPointPatchField
                <PatchField, PointPatch, CoupledPointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        CoupledPointPatchField
        (
            const CoupledPointPatchField
                <PatchField, PointPatch, CoupledPointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const Field<Type>& iF
        ) const = 0;


    // Member functions

        // Access

            //- Return true if this patch field is derived from
            //  coupledFvPatchField<Type>.
            virtual bool coupled() const
            {
                return true;
            }

            //- Evaluate the patch field
            virtual void evaluate() = 0;

            //- Initialise swap of patch point values
            virtual void initSwapAdd(Field<Type>&) const
            {}

            //- Complete swap of patch point values and add to local values
            virtual void swapAdd(Field<Type>&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "CoupledPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
