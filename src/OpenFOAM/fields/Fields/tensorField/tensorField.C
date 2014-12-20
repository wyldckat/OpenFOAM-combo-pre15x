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
    Specialisation of Field<T> for tensor.

\*---------------------------------------------------------------------------*/

#include "tensorField.H"
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION_N(scalar, tensor, tr)
UNARY_FUNCTION_R(tensor, tensor, dev)
UNARY_FUNCTION_R(tensor, tensor, dev2)
UNARY_FUNCTION_N(scalar, tensor, det)

void inv(Field<tensor>& tf, const UList<tensor>& tf1)
{
    if (!tf.size())
    {
        return;
    }

    tensorField tf1Plus(tf1);

    if (mag(tf1[0].xx()) < SMALL)
    {
        tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
    }

    if (mag(tf1[0].yy()) < SMALL)
    {
        tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
    }

    if (mag(tf1[0].zz()) < SMALL)
    {
        tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
    }

    TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1Plus)

    if (mag(tf1[0].xx()) < SMALL)
    {
        tf -= tensor(1,0,0,0,0,0,0,0,0);
    }

    if (mag(tf1[0].yy()) < SMALL)
    {
        tf -= tensor(0,0,0,0,1,0,0,0,0);
    }

    if (mag(tf1[0].zz()) < SMALL)
    {
        tf -= tensor(0,0,0,0,0,0,0,0,1);
    }
}

tmp<tensorField> inv(const UList<tensor>& tf)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result(), tf);
    return result;
}

tmp<tensorField> inv(const tmp<tensorField>& tf)
{
    tmp<tensorField> result(tf.ptr());
    inv(result(), result());
    return result;
}


UNARY_FUNCTION_R(tensor, tensor, hinv)
UNARY_FUNCTION_R(tensor, tensor, symm)
UNARY_FUNCTION_R(tensor, tensor, skew)
UNARY_FUNCTION_N(vector, tensor, eigenValues)
UNARY_FUNCTION_R(tensor, tensor, eigenVectors)


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR_N(vector, tensor, *, hdual)
UNARY_OPERATOR_N(tensor, vector, *, hdual)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
