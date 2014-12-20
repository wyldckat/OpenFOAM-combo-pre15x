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

Description
    Specialisation of Field<T> for tensor.

\*---------------------------------------------------------------------------*/

#include "tensorField.H"
#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

void hdual(vectorField& vf, const UList<tensor>& tf)
{
    TFOR_ALL_F_OP_OP_F(vector, vf, =, *, tensor, tf)
}

tmp<vectorField> operator*(const tmp<tensorField>& tf)
{
    tmp<vectorField> hDual(new vectorField(tf().size()));
    hdual(hDual(), tf);
    tf.clear();
    return hDual;
}


void hdual(tensorField& vf, const UList<vector>& tf)
{
    TFOR_ALL_F_OP_OP_F(tensor, vf, =, *, vector, tf)
}

tmp<tensorField> operator*(const tmp<vectorField>& tf)
{
    tmp<tensorField> hDual(new tensorField(tf().size()));
    hdual(hDual(), tf);
    tf.clear();
    return hDual;
}


// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

#define tensorFunc(func)                                                      \
void func(Field<tensor>& tf, const UList<tensor>& tf1)                        \
{                                                                             \
    TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, func, tensor, tf1)                    \
}                                                                             \
                                                                              \
tmp<tensorField> func(const UList<tensor>& tf)                                \
{                                                                             \
    tmp<tensorField> result(new tensorField(tf.size()));                      \
    func(result(), tf);                                                       \
    return result;                                                            \
}                                                                             \
                                                                              \
tmp<tensorField> func(const tmp<tensorField>& tf)                             \
{                                                                             \
    tmp<tensorField> result(tf.ptr());                                        \
    func(result(), result());                                                 \
    return result;                                                            \
}


void diag(vectorField& vf, const UList<tensor>& tf)
{
    TFOR_ALL_F_OP_FUNC_F(vector, vf, =, diag, tensor, tf)
}

tmp<vectorField> diag(const tmp<tensorField>& tf)
{
    tmp<vectorField> result(new vectorField(tf().size()));
    diag(result(), tf);
    tf.clear();
    return result;
}


void tr(scalarField& sf, const UList<tensor>& tf)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, sf, =, tr, tensor, tf)
}

tmp<scalarField> tr(const tmp<tensorField>& tf)
{
    tmp<scalarField> result(new scalarField(tf().size()));
    tr(result(), tf);
    tf.clear();
    return result;
}


tensorFunc(dev)
tensorFunc(dev2)

void det(scalarField& sf, const UList<tensor>& tf)
{
    TFOR_ALL_F_OP_FUNC_F(scalar, sf, =, det, tensor, tf)
}

tmp<scalarField> det(const tmp<tensorField>& tf)
{
    tmp<scalarField> result(new scalarField(tf().size()));
    det(result(), tf);
    tf.clear();
    return result;
}


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


tensorFunc(hinv)
tensorFunc(symm)
tensorFunc(skew)

void eigenValues(vectorField& vf, const UList<tensor>& tf)
{
    TFOR_ALL_F_OP_FUNC_F(vector, vf, =, eigenValues, tensor, tf)
}

tmp<vectorField> eigenValues(const tmp<tensorField>& tf)
{
    tmp<vectorField> result(new vectorField(tf().size()));
    eigenValues(result(), tf);
    tf.clear();
    return result;
}


tensorFunc(eigenVectors)

#undef tensorFunc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
