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
    PointPatchField

Description
    Abstract base class for tetrahedral mesh patch fields. Note the special
    mechanism at work here: the field itself holds no values, as the point
    values belong to the internal field. However, the field will create a
    list of values if required.

SourceFiles
    PointPatchField.C
    newPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef PointPatchField_H
#define PointPatchField_H

#include "primitiveFields.H"
#include "lduCoupledInterface.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

template<template<class> class PatchField, class PointPatch, class Type>
class PointPatchField;

class objectRegistry;
class dictionary;
class PointPatchFieldMapper;
class lduMatrix;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
Ostream& operator<<
(
    Ostream&,
    const PointPatchField<PatchField, PointPatch, Type>&
);


/*---------------------------------------------------------------------------*\
                       Class PointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class PointPatch, class Type>
class PointPatchField
:
    public lduCoupledInterface
{
    // Private data

        //- Reference to patch
        const PointPatch& patch_;

        //- Reference to internal field
        const Field<Type>& internalField_;


    // Private member functions

        void checkInternalField() const;


public:

    typedef PointPatch Patch;
    typedef PatchField<Type> pointPatchTypeField;

    //- Runtime type information
    TypeName("PointPatchField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchTypeField,
            PointPatch,
            (const PointPatch& p, const Field<Type>& iF),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchTypeField,
            patchMapper,
            (
                const PointPatchField<PatchField, PointPatch, Type>& ptf,
                const PointPatch& p,
                const Field<Type>& iF,
                const PointPatchFieldMapper& m
            ),
            (dynamic_cast<const pointPatchTypeFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchTypeField,
            dictionary,
            (
                const PointPatch& p,
                const Field<Type>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        PointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct as copy
        PointPatchField(const PointPatchField<PatchField, PointPatch, Type>&);

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const = 0;

        //- Construct as copy setting internal field reference
        PointPatchField
        (
            const PointPatchField<PatchField, PointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone
        (
            const Field<Type>& iF
        ) const = 0;


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static autoPtr<PatchField<Type> > New
        (
            const word&,
            const PointPatch&,
            const Field<Type>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given PointPatchField mapped onto a new patch
        static autoPtr<PatchField<Type> > New
        (
            const PointPatchField<PatchField, PointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static autoPtr<PatchField<Type> > New
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );


        //- Return a pointer to a new CalculatedPointPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static autoPtr<PointPatchField<PatchField, PointPatch, Type> >
        NewCalculatedType
        (
            const PointPatchField<PatchField, PointPatch, Type2>&
        );


    // Destructor

        virtual ~PointPatchField<PatchField, PointPatch, Type>()
        {}


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return size
            label size() const
            {
                return patch().size();
            }

            //- Return patch
            const PointPatch& patch() const
            {
                return patch_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return field created from appropriate internal field values
            tmp<Field<Type> > patchInternalField() const;

            //- Return field created from appropriate internal field values
            // given internal field reference
            template<class Type1>
            tmp<Field<Type1> > patchInternalField
            (
                const Field<Type1>& iF
            ) const;

            //- Given the internal field and a patch field,
            //  add the patch field to the internal field
            template<class Type1> void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1> void setInInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Return the type of the calculated for of PointPatchField
            static const word& calculatedType();


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
            )
            {}

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const PointPatchField<PatchField, PointPatch, Type>&,
                const labelList&
            )
            {}


        // Evaluation functions

            //- Does it store field data
            virtual bool storesFieldData() const
            {
                return false;
            }

            //- Does this patchField correspond to a pointTypeField
            bool isPointField() const;

            //- Check that this patchField corresponds to a pointTypeField,
            // if not abort!
            void checkPointField() const;

            //- Initialise evaluation of the patch field (do nothing)
            virtual void initEvaluate(const bool=false)
            {}

            //- Evaluate the patch field
            virtual void evaluate()
            {}

            //- Update boundary value
            virtual void updateBoundaryField()
            {}

            //- Init add field to internal field
            virtual void initAddField() const
            {}

            //- Add field to internal field
            virtual void addField(Field<Type>&) const
            {
                notImplemented(type() + "::addField(Field<Type>&) const");
            }

            //- Return neighbour colouring
            virtual tmp<labelField> nbrColour(const labelField& cField) const
            {
                notImplemented(type() + "::nbrColour(const labelField&) const");

                return tmp<labelField>(labelField::null());
            }


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=
        (
            const PointPatchField<PatchField, PointPatch, Type>&
        )
        {}

        virtual void operator+=
        (
            const PointPatchField<PatchField, PointPatch, Type>&
        )
        {}

        virtual void operator-=
        (
            const PointPatchField<PatchField, PointPatch, Type>&
        )
        {}

        virtual void operator*=
        (
            const PointPatchField<PatchField, PointPatch, scalar>&
        )
        {}

        virtual void operator/=
        (
            const PointPatchField<PatchField, PointPatch, scalar>&
        )
        {}

        virtual void operator=(const Field<Type>&){}
        virtual void operator+=(const Field<Type>&){}
        virtual void operator-=(const Field<Type>&){}

        virtual void operator*=(const Field<scalar>&){}
        virtual void operator/=(const Field<scalar>&){}

        virtual void operator=(const Type&){}
        virtual void operator+=(const Type&){}
        virtual void operator-=(const Type&){}
        virtual void operator*=(const scalar){}
        virtual void operator/=(const scalar){}


        // Force an assignment irrespective of form of patch

        virtual void operator==
        (
            const PointPatchField<PatchField, PointPatch, Type>&
        )
        {}

        virtual void operator==(const Field<Type>&){}
        virtual void operator==(const Type&){}


    // Ostream operator

        friend Ostream& operator<< <PatchField, PointPatch, Type>
        (
            Ostream&,
            const PointPatchField<PatchField, PointPatch, Type>&
        );
};


// This function is added as a hack to enable simple backward compatability
// with verions using referenceLevel in GeometicField
template<template<class> class PatchField, class PointPatch, class Type>
const PointPatchField<PatchField, PointPatch, Type>& operator+
(
    const PointPatchField<PatchField, PointPatch, Type>& ppf,
    const Type&
)
{
    return ppf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PointPatchFieldFunctions.H"

#ifdef NoRepository
#   include "PointPatchField.C"
#   include "CalculatedPointPatchField.H"
#endif


#define makePointPatchTypeFieldTypeName(typePatchTypeField)                \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);

#define makePointPatchFieldsTypeName(typePatchField, pointPatchType)       \
                                                                           \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##ScalarField); \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##VectorField); \
makePointPatchTypeFieldTypeName(typePatchField##pointPatchType##TensorField);


#define makePointPatchTypeField(PatchTypeField, typePatchTypeField)        \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, PointPatch                         \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField,                                                        \
    typePatchTypeField,                                                    \
    patchMapper                                                            \
);                                                                         \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, dictionary                         \
);


#define makePointPatchFields(typePatchField, pointPatch, pointPatchType)   \
                                                                           \
makePointPatchTypeField                                                    \
(                                                                          \
    pointPatch##ScalarField,                                               \
    typePatchField##pointPatchType##ScalarField                            \
);                                                                         \
                                                                           \
makePointPatchTypeField                                                    \
(                                                                          \
    pointPatch##VectorField,                                               \
    typePatchField##pointPatchType##VectorField                            \
);                                                                         \
                                                                           \
makePointPatchTypeField                                                    \
(                                                                          \
    pointPatch##TensorField,                                               \
    typePatchField##pointPatchType##TensorField                            \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
