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
    faPatchField

Description
    faPatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation proceedures, particularly with respect to specific fields.

SourceFiles
    faPatchField.C
    newPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef faPatchField_H
#define faPatchField_H

#include "faPatch.H"
#include "primitiveFields.H"
#include "lduCoupledInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;
class dictionary;
class faPatchFieldMapper;
class lduMatrix;

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class faPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const faPatchField<Type>&);

template<class Type>
class calculatedFaPatchField;


/*---------------------------------------------------------------------------*\
                           Class patch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faPatchField
:
    public Field<Type>,
    virtual public lduCoupledInterface
{
    // Private data

        const faPatch& patch_;
        const Field<Type>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;

    // Private member functions

        void checkPatch() const;
        void checkInternalField() const;


public:

    typedef faPatch Patch;


    //- Runtime type information
    TypeName("faPatchField");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patch,
            (const faPatch& p, const Field<Type>& iF),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patchMapper,
            (
                const faPatchField<Type>& ptf,
                const faPatch& p,
                const Field<Type>& iF,
                const faPatchFieldMapper& m
            ),
            ((const faPatchFieldType&)ptf, p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            dictionary,
            (const faPatch& p, const Field<Type>& iF, const dictionary& dict),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        faPatchField
        (
            const faPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        faPatchField
        (
            const faPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        faPatchField
        (
            const faPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping the given faPatchField onto a new patch
        faPatchField
        (
            const faPatchField<Type>&,
            const faPatch&,
            const Field<Type>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        faPatchField(const faPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >(new faPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        faPatchField(const faPatchField<Type>&, const Field<Type>&);

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone(const Field<Type>& iF) const
        {
            return tmp<faPatchField<Type> >(new faPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static tmp<faPatchField<Type> > New
        (
            const word&,
            const faPatch&,
            const Field<Type>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given faPatchField mapped onto a new patch
        static tmp<faPatchField<Type> > New
        (
            const faPatchField<Type>&,
            const faPatch&,
            const Field<Type>&,
            const faPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static tmp<faPatchField<Type> > New
        (
            const faPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFaPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static tmp<calculatedFaPatchField<Type> > NewCalculatedType
        (
            const faPatchField<Type2>&
        );


    // Destructor

        virtual ~faPatchField<Type>()
        {}


    // Member functions

        // Access

            //- Return database
            const objectRegistry& db() const;

            //- Return patch
            const faPatch& patchMesh() const
            {
                return patch_;
            }

            //- Return internal field reference
            const Field<Type>& internalField() const
            {
                return internalField_;
            }

            //- Return the type of the calculated for of faPatchField
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

            //- Does this patchField correspond to a areaTypeField
            bool areaPatchField() const;

            //- Check that this patchField corresponds to a areaTypeField,
            // if not abort!
            void checkAreaField() const;

            //- Return the corresponding patchField of the named field
            template<class GeometricField, class Type2>
            const faPatchField<Type2>& patchField(const GeometricField&) const;

            //- Lookup and return the patchField of the named field from the DB
            template<class GeometricField, class Type2>
            const faPatchField<Type2>& lookupPatchField
            (
                const word& name,
                const GeometricField* = NULL,
                const Type2* = NULL
            ) const;


        // Mapping functions

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
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
            virtual void initEvaluate(const bool)
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

            //- Check faPatchField<Type> against given faPatchField<Type>
            void check(const faPatchField<Type>&) const;


    // Member operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const faPatchField<Type>&);
        virtual void operator+=(const faPatchField<Type>&);
        virtual void operator-=(const faPatchField<Type>&);
        virtual void operator*=(const faPatchField<scalar>&);
        virtual void operator/=(const faPatchField<scalar>&);

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

        virtual void operator==(const faPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const faPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faPatchField.C"
#   include "calculatedFaPatchField.H"
#endif


#define makeFaPatchTypeFieldTypeName(typePatchTypeField)                   \
                                                                           \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);

#define makeFaPatchFieldsTypeName(typePatchField)                          \
                                                                           \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchScalarField);          \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchVectorField);          \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchTensorField);


#define makeFaPatchTypeField(PatchTypeField, typePatchTypeField)           \
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


#define makeFaPatchFields(typePatchField)                                     \
                                                                              \
makeFaPatchTypeField(faPatchScalarField, typePatchField##FaPatchScalarField); \
makeFaPatchTypeField(faPatchVectorField, typePatchField##FaPatchVectorField); \
makeFaPatchTypeField(faPatchTensorField, typePatchField##FaPatchTensorField);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
