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

Type
    typeInfo

Description
    Basic run-time type information using word as the type's name. 
    Used to enhance the standard RTTI to cover I/O.

    The user can get the type's type name using the type info access function

        type()

    The reference type cast template function:

        refCast<T>(r)

    wraps dynamic_cast to handle the bad_cast exception and generate a
    FatalError.

\*---------------------------------------------------------------------------*/

#ifndef typeInfo_H
#define typeInfo_H

#include "error.H"
#include "className.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TypeName(TypeNameString)                              \
    ClassName(TypeNameString);                                \
    virtual const word& type() const { return typeName; }


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Reference type cast template function wraps dynamic_cast to handle the
// bad_cast exception and generate a FatalError.

template<class To, class From>
inline To& dynamicCast(From& r)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (std::bad_cast)
    {
        FatalErrorIn("dynamicCast<To>(From&)")
            << "Attempt to cast type " << typeid(r).name()
            << " to type " << typeid(To).name()
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


template<class To, class From>
inline To& refCast(From& r)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (std::bad_cast)
    {
        FatalErrorIn("refCast<To>(From&)")
            << "Attempt to cast type " << r.type()
            << " to type " << To::typeName
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
