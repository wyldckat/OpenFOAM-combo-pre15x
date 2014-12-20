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
    ComponentMixedTetPointPatchVectorField. The boundary condition is a mix
    of a fixedValue and a zeroGradient boundary condition, where a
    fixedValue/zeroGradient mix may be different for each direction.

    I am still not sure how to do the fixedValue-fixedGradient
    combination. 

Description

SourceFiles
    ComponentMixedTetPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef ComponentMixedTetPointPatchVectorField_H
#define ComponentMixedTetPointPatchVectorField_H

#include "tetPolyPatchFields.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class ComponentMixedTetPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class PatchField,
    class PointPatch
>
class ComponentMixedTetPointPatchVectorField
:
    public PatchField<vector>
{
    // Private data

        //- Fraction vector of value used for boundary condition
        vectorField refValue_;

        //- Fraction vector of value used for boundary condition
        vectorField valueFraction_;


    // Private member functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("componentMixed");


    // Constructors

        //- Construct from patch and internal field
        ComponentMixedTetPointPatchVectorField
        (
            const PointPatch&,
            const vectorField&
        );

        //- Construct from patch and internal field, value
        // and value fraction
        ComponentMixedTetPointPatchVectorField
        (
            const PointPatch&,
            const vectorField&,
            const vectorField&,
            const vectorField&
        );

        //- Construct from patch, internal field and dictionary
        ComponentMixedTetPointPatchVectorField
        (
            const PointPatch&,
            const vectorField&,
            const dictionary&
        );

        //- Construct by mapping given patchVectorField onto a new patch
        ComponentMixedTetPointPatchVectorField
        (
            const ComponentMixedTetPointPatchVectorField<PatchField, PointPatch>&,
            const PointPatch&,
            const vectorField&,
            const TetPointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<vector> > clone() const
        {
            return autoPtr<PatchField<vector> >
            (
                new ComponentMixedTetPointPatchVectorField
                <PatchField, PointPatch>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        ComponentMixedTetPointPatchVectorField
        (
            const ComponentMixedTetPointPatchVectorField<PatchField, PointPatch>&,
            const vectorField&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<vector> > clone(const vectorField& iF) const
        {
            return autoPtr<PatchField<vector> >
            (
                new ComponentMixedTetPointPatchVectorField
                <PatchField, PointPatch>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Return defining fields

            virtual vectorField& refValue()
            {
                return refValue_;
            }

            virtual const vectorField& refValue() const
            {
                return refValue_;
            }


            virtual vectorField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const vectorField& valueFraction() const
            {
                return valueFraction_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const TetPointPatchFieldMapper&
            );

            //- Reverse map the given TetPointPatchVectorField onto
            //  this TetPointPatchVectorField
            virtual void rmap
            (
                const TetPointPatchField<PatchField, PointPatch, vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Insert boundary value into the internal field
            virtual void evaluate();

            //- Set boundary condition to matrix
            virtual void setBoundaryCondition
            (
                Map<constraint<vector> > &
            ) const;


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "ComponentMixedTetPointPatchVectorField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
