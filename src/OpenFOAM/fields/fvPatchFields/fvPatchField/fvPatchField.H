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
    fvPatchField

Description
    fvPatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation proceedures, particularly with respect to specific fields.

SourceFiles
    fvPatchField.C
    newPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchField_H
#define fvPatchField_H

#include "fvPatch.H"
#include "primitiveFields.H"
#include "lduCoupledInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;
class dictionary;
class fvPatchFieldMapper;
class lduMatrix;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class fvPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const fvPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                           Class patch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvPatchField
:
    virtual public lduCoupledInterface,
    public Field<Type>
{
    // Private data

        const fvPatch& patch_;
        const Field<Type>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;


    // Private member functions

        void checkPatch() const;
        void checkInternalField() const;


public:

    typedef fvPatch Patch;


    //- Runtime type information
    TypeName("fvPatchField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patch,
            (const fvPatch& p, const Field<Type>& iF),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            patchMapper,
            (
                const fvPatchField<Type>& ptf,
                const fvPatch& p,
                const Field<Type>& iF,
                const fvPatchFieldMapper& m
            ),
            (dynamic_cast<const fvPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            fvPatchField,
            dictionary,
            (const fvPatch& p, const Field<Type>& iF, const dictionary& dict),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        fvPatchField
        (
            const fvPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        fvPatchField
        (
            const fvPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        fvPatchField
        (
            const fvPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping the given fvPatchField onto a new patch
        fvPatchField
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const Field<Type>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvPatchField(const fvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >(new fvPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        fvPatchField(const fvPatchField<Type>&, const Field<Type>&);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone(const Field<Type>& iF) const
        {
            return tmp<fvPatchField<Type> >(new fvPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<fvPatchField<Type> > New
        (
            const word&,
            const fvPatch&,
            const Field<Type>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given fvPatchField mapped onto a new patch
        static tmp<fvPatchField<Type> > New
        (
            const fvPatchField<Type>&,
            const fvPatch&,
            const Field<Type>&,
            const fvPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<fvPatchField<Type> > New
        (
            const fvPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFvPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<fvPatchField<Type> > NewCalculatedType
        (
            const fvPatchField<Type2>&
        );


    // Destructor

        virtual ~fvPatchField<Type>()
        {}


    // Member functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const fvPatch& patch() const
            {
                return patch_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return the type of the calculated for of fvPatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }

            //- Does this patchField correspond to a volTypeField
            bool isVolField() const;

            //- Check that this patchField corresponds to a volTypeField,
            // if not abort!
            void checkVolField() const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type2>
            const fvPatchField<Type2>& patchField(const GeometricField&) const;

            //- Lookup and return the patchField of the named field from the
            //  local objectRegistry.
            //  N.B.  The dummy pointer arguments are used if this function is 
            //  instantiated within a templated function to avoid a bug in gcc.
            //  See inletOutletFvPatchField.C and outletInletFvPatchField.C
            template<class GeometricField, class Type2>
            const fvPatchField<Type2>& lookupPatchField
            (
                const word& name,
                const GeometricField* = NULL,
                const Type2* = NULL
            ) const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad() const;

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Return given internal field next to patch as patch field
            template<class Type2>
            tmp<Field<Type2> > patchInternalField
            (
                const Field<Type2>&
            ) const;

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type> > patchInternalField() const;

            //- Return neighbour coupled given internal cell data
            virtual tmp<Field<Type> > patchNeighbourField
            (
                const Field<Type>& iField
            ) const;

            //- Return neighbour colouring.
            //  Needed for AMG solver agglomeration
            virtual tmp<labelField> nbrColour
            (
                const labelField& iColour
            ) const;

            //- Return patchField of the values on the patch or on the
            //  opposite patch
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const bool=false)
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate();


            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<Field<scalar> >&
            ) const
            {
                notImplemented
                (
                    type()
                  + "::valueInternalCoeffs(const tmp<Field<scalar> >&)"
                );
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<Field<scalar> >&
            ) const
            {
                notImplemented
                (
                    type()
                  + "::valueBoundaryCoeffs(const tmp<Field<scalar> >&)"
                );
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const
            {
                notImplemented(type() + "::gradientInternalCoeffs()");
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const
            {
                notImplemented(type() + "::gradientBoundaryCoeffs()");
                return *this;
            }


        //- Write
        virtual void write(Ostream&) const;


        // Check

            //- Check fvPatchField<Type> against given fvPatchField<Type>
            void check(const fvPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const fvPatchField<Type>&);
        virtual void operator+=(const fvPatchField<Type>&);
        virtual void operator-=(const fvPatchField<Type>&);
        virtual void operator*=(const fvPatchField<scalar>&);
        virtual void operator/=(const fvPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const fvPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const fvPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvPatchField.C"
#   include "calculatedFvPatchField.H"
#endif


#define makePatchTypeFieldTypeName(typePatchTypeField)                     \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);

#define makePatchFieldsTypeName(typePatchField)                            \
                                                                           \
makePatchTypeFieldTypeName(typePatchField##FvPatchScalarField);            \
makePatchTypeFieldTypeName(typePatchField##FvPatchVectorField);            \
makePatchTypeFieldTypeName(typePatchField##FvPatchTensorField);


#define makePatchTypeField(PatchTypeField, typePatchTypeField)             \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
                                                                           \
addToRunTimeSelectionTable                                                 \
(                                                                          \
    PatchTypeField, typePatchTypeField, patch                              \
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


#define makePatchFields(typePatchField)                                    \
                                                                           \
makePatchTypeField(fvPatchScalarField, typePatchField##FvPatchScalarField);\
makePatchTypeField(fvPatchVectorField, typePatchField##FvPatchVectorField);\
makePatchTypeField(fvPatchTensorField, typePatchField##FvPatchTensorField);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
