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

Description
    Specialisation of Field<T> for sphericalTensor.

\*---------------------------------------------------------------------------*/

#include "sphericalTensorField.H"
#include "transformField.H"
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION_N(scalar, sphericalTensor, tr)
UNARY_FUNCTION_N(scalar, sphericalTensor, det)
UNARY_FUNCTION_R(sphericalTensor, sphericalTensor, inv)

BINARY_OPERATOR_NR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR_NR(sphericalTensor, scalar, sphericalTensor, /, divide)


template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tensorField& tf
)
{
    return tmp<Field<sphericalTensor> >
    (
        new Field<sphericalTensor>
        (
            tf.size(), sphericalTensor::one
        )
    );
}

template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<sphericalTensor> > ret =
        transformFieldMask<sphericalTensor>(ttf());
    ttf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
