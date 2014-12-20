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

\*---------------------------------------------------------------------------*/

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const dimensionSet& dims,
    const Field<Type>& field
)
:
    regIOobject(io),
    Field<Type>(field),
    dimensions_(dims)
{}


// Construct as copy
template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const DimensionedField<Type, GeoMesh>& df
)
:
    regIOobject(df),
    Field<Type>(df),
    dimensions_(df.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::~DimensionedField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate and return sum value average
template<class Type, class GeoMesh>
dimensioned<Type> DimensionedField<Type, GeoMesh>::average() const
{
    return dimensioned<Type>
    (
        name() + ".average()",
        dimensions_,
        gAverage(*this)
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::operator=
(
    const DimensionedField<Type, GeoMesh>& df
)
{
    // Check for assignment to self
    if (this == &df)
    {
        FatalErrorIn
        (
            "DimensionedField<Type, GeoMesh>::operator="
            "(const DimensionedField<Type, GeoMesh>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    dimensions() = df.dimensions();
    Field<Type>::operator=(df);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldIO.C"

// ************************************************************************* //
