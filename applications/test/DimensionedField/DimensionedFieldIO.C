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

#include "DimensionedField.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from IOobject
template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const word& fieldDictEntry
)
:
    regIOobject(io),
    Field<Type>(0),
    dimensions_(dimless)
{
    Istream& is = readStream(typeName);

    if (is.version() < 2.0)
    {
        dimensions_.reset(dimensionSet(is));
        operator>>(is, static_cast<Field<Type>&>(*this));
    }
    else
    {
        dictionary fieldDict(is);

        dimensions_.reset(dimensionSet(fieldDict.lookup("dimensions")));

        Field<Type>::operator=
        (
            tmp<Field<Type> >
            (
                new Field<Type>
                (
                    GeoMesh().size(),
                    pTraits<Type>(fieldDict.lookup(fieldDictEntry))
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
bool DimensionedField<Type, GeoMesh>::writeData(Ostream& os) const
{
    os  << dimensions() << nl
        << static_cast<const Field<Type>&>(*this)
        << endl;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const DimensionedField<Type, GeoMesh>&)"
    );

    return (os.good());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Ostream& operator<<(Ostream& os, const DimensionedField<Type, GeoMesh>& df)
{
    df.writeData(os);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
