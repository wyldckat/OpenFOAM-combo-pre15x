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
    BasicSymmetryTetPointPatchField

Description

SourceFiles
    BasicSymmetryTetPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef BasicSymmetryTetPointPatchField_H
#define BasicSymmetryTetPointPatchField_H

#include "ValueStoredTetPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class BasicSymmetryTetPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class PointPatch, class Type>
class BasicSymmetryTetPointPatchField
:
    public ValueStoredTetPointPatchField<PatchField, PointPatch, Type>
{
    // Private memeber functions

        //- Calculate implicitness matrix
        Type fixFlagMatrix(const vector&) const;


public:

    // Constructors

        //- Construct from patch and internal field
        BasicSymmetryTetPointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        BasicSymmetryTetPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        BasicSymmetryTetPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        BasicSymmetryTetPointPatchField
        (
            const BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const TetPointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        BasicSymmetryTetPointPatchField
        (
            const BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone(const Field<Type>& iF) const
        {
            return autoPtr<PatchField<Type> >
            (
                new BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the patch field
        virtual void updateBoundaryField();

        //- Set boundary condition to matrix
        virtual void setBoundaryCondition
        (
            Map<constraint<Type> > &
        ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "BasicSymmetryTetPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
