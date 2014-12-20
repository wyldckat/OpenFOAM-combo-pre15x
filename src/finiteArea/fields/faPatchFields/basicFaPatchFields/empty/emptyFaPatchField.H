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
    emptyFaPatchField

Description

SourceFiles
    emptyFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef emptyFaPatchField_H
#define emptyFaPatchField_H

#include "faPatchField.H"
#include "emptyFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class emptyFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class emptyFaPatchField
:
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("empty");


    // Constructors

        //- Construct from patch and internal field
        emptyFaPatchField
        (
            const faPatch&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        emptyFaPatchField
        (
            const faPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given emptyFaPatchField onto a new patch
        emptyFaPatchField
        (
            const emptyFaPatchField<Type>&,
            const faPatch&,
            const Field<Type>&,
            const faPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type> > clone() const
        {
            return tmp<faPatchField<Type> >
            (
                new emptyFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        emptyFaPatchField
        (
            const emptyFaPatchField<Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type> > clone(const Field<Type>& iF) const
        {
            return tmp<faPatchField<Type> >
            (
                new emptyFaPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type> > gradientInternalCoeffs() const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            tmp<Field<Type> > gradientBoundaryCoeffs() const
            {
                return tmp<Field<Type> >(new Field<Type>(0));
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "emptyFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
