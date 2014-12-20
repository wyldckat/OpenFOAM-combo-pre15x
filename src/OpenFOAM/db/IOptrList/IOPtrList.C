/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    An IOPtrList of a given type is a PtrList of that type which supports
    automated input and output.

\*---------------------------------------------------------------------------*/

#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
template<class INew>
IOPtrList<T>::IOPtrList(const IOobject& io, const INew& inewt)
:
    regIOobject(io),
    PtrList<T>(readStream(typeName), inewt)
{
    close();
}


template<class T>
IOPtrList<T>::IOPtrList(const IOobject& io)
:
    regIOobject(io),
    PtrList<T>(readStream(typeName))
{
    close();
}


template<class T>
IOPtrList<T>::IOPtrList(const IOobject& io, const PtrList<T>& list)
:
    regIOobject(io),
    PtrList<T>(list)
{}


template<class T>
void IOPtrList<T>::operator=(const IOPtrList<T>& rhs)
{
    PtrList<T>::operator=(rhs);
}


template<class T>
bool IOPtrList<T>::writeData(Ostream& os) const
{
    return (os << *this).good();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
