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

Class
    FieldField<T>

Description
    Generic field type.

SourceFiles
    FieldField.C

\*---------------------------------------------------------------------------*/

#ifndef FieldField_H
#define FieldField_H

#include "tmp.H"
#include "ptrList.H"
#include "scalar.H"
#include "direction.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<template<class> class Field, class Type>
class FieldField;

template<template<class> class Field, class Type>
Ostream& operator<<
(
    Ostream&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
Ostream& operator<<
(
    Ostream&,
    const tmp<FieldField<Field, Type> >&
);


/*---------------------------------------------------------------------------*\
                           Class FieldField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class Field, class Type>
class FieldField
:
    public refCount,
    public ptrList<Field<Type> >
{

public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        FieldField();

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit FieldField(const label);

        //- Construct using the Field sizes from the given FieldField
        //  and the given Field type.
        //  Used for temporary fields which are initialised after construction
        FieldField(const word&, const FieldField<Field, Type>&);

        //- Construct as copy
        FieldField(const FieldField<Field, Type>&);

        //- Construct as copy of a ptrList<Field, Type>
        FieldField(const ptrList<Field<Type> >&);

        //- Construct as copy of tmp<FieldField>
#       ifdef ConstructFromTmp
        FieldField(const tmp<FieldField<Field, Type> >&);
#       endif

        //- Construct from Istream
        FieldField(Istream&);

        //- Clone
        tmp<FieldField<Field, Type> > clone() const;

        //- Return a pointer to a new calculatedFvPatchFieldField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<FieldField<Field, Type> > NewCalculatedType
        (
            const FieldField<Field, Type2>& ff
        )
        {
            FieldField<Field, Type>* nffPtr
            (
                new FieldField<Field, Type>(ff.size())
            );

            forAll(*nffPtr, i)
            { 
                nffPtr->hook(Field<Type>::NewCalculatedType(ff[i]).ptr());
            }

            return tmp<FieldField<Field, Type> >(nffPtr);
        }


    // Member functions

        //- Negate this field
        void negate();

        //- Return a component field of the field
        tmp<FieldField<Field, cmptType> > component(const direction) const;

        //- Replace a component field of the field
        void replace(const direction, const FieldField<Field, cmptType>&);

        //- Return the field transpose (only defined for second rank tensors)
        tmp<FieldField<Field, Type> > T() const;


    // Member operators

        void operator=(const FieldField<Field, Type>&);
        void operator=(const tmp<FieldField<Field, Type> >&);
        void operator=(const Type&);

        void operator+=(const FieldField<Field, Type>&);
        void operator+=(const tmp<FieldField<Field, Type> >&);

        void operator-=(const FieldField<Field, Type>&);
        void operator-=(const tmp<FieldField<Field, Type> >&);

        void operator*=(const FieldField<Field, scalar>&);
        void operator*=(const tmp<FieldField<Field, scalar> >&);

        void operator/=(const FieldField<Field, scalar>&);
        void operator/=(const tmp<FieldField<Field, scalar> >&);

        void operator+=(const Type&);
        void operator-=(const Type&);

        void operator*=(const scalar&);
        void operator/=(const scalar&);


    // IOstream operators

        friend Ostream& operator<<
#       ifndef KAI
            <Field, Type>
#       endif
        (
            Ostream&,
            const FieldField<Field, Type>&
        );

        friend Ostream& operator<<
#       ifndef KAI
            <Field, Type>
#       endif
        (
            Ostream&,
            const tmp<FieldField<Field, Type> >&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldFieldFunctions.H"

#ifdef NoRepository
#   include "FieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
