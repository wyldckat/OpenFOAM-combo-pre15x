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

\*---------------------------------------------------------------------------*/

#include "tetPolyPatchInterpolationFaceDecomp.H"
#include "faceTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolationFaceDecomp::faceToPointInterpolate
(
    const Field<Type>& ff
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>(patch_.size())
    );
    Field<Type>& result = tresult();

    // Insert the point values first
    label i = 0;

    const Field<Type> pointRes = interpolator_.faceToPointInterpolate(ff);

    forAll (pointRes, pointI)
    {
        result[i] = pointRes[pointI];
        i++;
    }

    // Insert the face centre values; no interpolation necessary
    forAll (ff, faceI)
    {
        result[i] = ff[faceI];
        i++;
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::tetPolyPatchInterpolationFaceDecomp::faceToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = faceToEdgeInterpolate(tff());
    tff.clear();
    return tint;
}


// ************************************************************************* //
