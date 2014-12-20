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

#include "DictionaryBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
template<class IDLListType, class T>
DictionaryBase<IDLListType, T>::DictionaryBase()
{}


// Copy constructor
template<class IDLListType, class T>
DictionaryBase<IDLListType, T>::DictionaryBase(const DictionaryBase& dict)
:
    IDLListType(dict)
{
    for
    (
        typename IDLListType::iterator iter = this->begin();
        iter != this->end();
        ++iter
    )
    {
        this->hashedEntries_.insert((*iter).keyword(), &(*iter));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Find and return T
template<class IDLListType, class T>
bool DictionaryBase<IDLListType, T>::found(const word& keyword) const
{
    return hashedTs_.found(keyword);
}


// Find and return T*
template<class IDLListType, class T>
const T* DictionaryBase<IDLListType, T>::lookup(const word& keyword) const
{
    typename HashTable<T*>::const_iterator iter = hashedTs_.find(keyword);

    if (iter == hashedTs_.end())
    {
        // If keyword not found print error message ...
        FatalErrorIn
        (
            "DictionaryBase<IDLListType, T>::"
            "lookup(const word& keyword) const"
        )   << keyword << " is undefined"
            << exit(FatalError);
    }

    return *iter;
}


// Find and return T*
template<class IDLListType, class T>
T* DictionaryBase<IDLListType, T>::lookup(const word& keyword)
{
    typename HashTable<T*>::iterator iter = hashedTs_.find(keyword);

    if (iter == hashedTs_.end())
    {
        // If keyword not found print error message ...
        FatalErrorIn
        (
            "DictionaryBase<IDLListType, T>::lookup(const word& keyword)"
        )   << keyword << " is undefined"
            << exit(FatalError);
    }

    return *iter;
}


// Return the table of contents
template<class IDLListType, class T>
wordList DictionaryBase<IDLListType, T>::toc() const
{
    return hashedTs_.toc();
}


// Add at head of dictionary
template<class IDLListType, class T>
void DictionaryBase<IDLListType, T>::insert(const word& keyword, T* tPtr)
{
    IDLListType::insert(tPtr);
    hashedTs_.insert(keyword, tPtr);
}


// Add at tail of dictionary
template<class IDLListType, class T>
void DictionaryBase<IDLListType, T>::append(const word& keyword, T* tPtr)
{
    IDLListType::append(tPtr);
    hashedTs_.insert(keyword, tPtr);
}


template<class IDLListType, class T>
T* DictionaryBase<IDLListType, T>::remove(const word& Keyword)
{
    typename HashTable<T*>::iterator iter = hashedTs_.find(Keyword);

    if (iter != hashedTs_.end())
    {
        T* tPtr = IDLListType::remove(iter());
        hashedTs_.erase(iter);
        return tPtr;
    }
    else
    {
        return NULL;
    }
}


//- Clear the dictionary
template<class IDLListType, class T>
void DictionaryBase<IDLListType, T>::clear()
{
    IDLListType::clear();
    hashedTs_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class IDLListType, class T>
void DictionaryBase<IDLListType, T>::operator=
(
    const DictionaryBase<IDLListType, T>& dict
)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("DictionaryBase::operator=(const DictionaryBase&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    IDLListType::operator=(dict);

    this->hashedEntries_.clear();

    for
    (
        typename IDLListType::iterator iter = this->begin();
        iter != this->end();
        ++iter
    )
    {
        this->hashedEntries_.insert((*iter).keyword(), &(*iter));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DictionaryBaseIO.C"

// ************************************************************************* //
