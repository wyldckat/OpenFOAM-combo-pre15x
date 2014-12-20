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
    HashTable

Description
    STL conforming hash table.

\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline HashTable<T, Key, Hash>::hashedEntry::hashedEntry
(
    const Key& key,
    hashedEntry* next,
    T newEntry
)
:
    key_(key),
    next_(next),
    obj_(newEntry)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline label HashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& HashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("HashTable<T, Key, Hash>::operator[](const Key&)")
            << key << " not found in table.  Valid entries are "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}

template<class T, class Key, class Hash>
inline const T& HashTable<T, Key, Hash>::operator[](const Key& key) const
{
    const_iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("HashTable<T, Key, Hash>::operator[](const Key&) const")
            << key << " not found in table.  Valid entries are "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
Iterator
(
    TableRef curHashTable,
    HashedEntryPtr elmt,
    HashedEntryPtr prev,
    label hashIndex
)
:
    curHashTable_(curHashTable),
    elmtPtr_(elmt),
    prevElmtPtr_(prev),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
Iterator(const iterator& iter)
:
    curHashTable_(iter.curHashTable_),
    elmtPtr_(iter.elmtPtr_),
    prevElmtPtr_(iter.prevElmtPtr_),
    hashIndex_(iter.hashIndex_)
{}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline void HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator=(const iterator& iter)
{
    this->elmtPtr_ = iter.elmtPtr_;
    this->prevElmtPtr_ = iter.prevElmtPtr_;
    this->hashIndex_ = iter.hashIndex_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline bool HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator==(const iterator& iter) const
{
    return elmtPtr_ == iter.elmtPtr_;
}

template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline bool HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator==(const const_iterator& iter) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline bool HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator!=(const iterator& iter) const
{
    return elmtPtr_ != iter.elmtPtr_;
}

template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline bool HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator!=(const const_iterator& iter) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline TRef HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator*()
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline TRef HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
operator()()
{
    return operator*();
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline
typename HashTable<T, Key, Hash>::template Iterator
<
    TRef,
    TableRef,
    HashedEntryPtr
>&
HashTable<T, Key, Hash>::Iterator
<
    TRef,
    TableRef,
    HashedEntryPtr
>::operator++()
{
    prevElmtPtr_ = elmtPtr_;

    if
    (
        !(elmtPtr_ = elmtPtr_->next_)
     && ++hashIndex_ < curHashTable_.tableSize_
     && !(elmtPtr_ = curHashTable_.table_[hashIndex_])
    )
    {
        prevElmtPtr_ = 0;

        while
        (
            ++hashIndex_ < curHashTable_.tableSize_
         && !(elmtPtr_ = curHashTable_.table_[hashIndex_])
        );
    }

    return *this;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline
typename HashTable<T, Key, Hash>::template Iterator
<
    TRef,
    TableRef,
    HashedEntryPtr
>
HashTable<T, Key, Hash>::Iterator
<
    TRef,
    TableRef,
    HashedEntryPtr
>::operator++
(
    int
)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef, class HashedEntryPtr>
inline
const Key& HashTable<T, Key, Hash>::Iterator<TRef, TableRef, HashedEntryPtr>::
key()
{
    return elmtPtr_->key_;
}


template<class T, class Key, class Hash>
inline typename HashTable<T, Key, Hash>::iterator
HashTable<T, Key, Hash>::begin()
{
    label i = 0;

    while (table_ && !table_[i] && ++i < tableSize_);

    if (i == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable is empty\n";
        }
#       endif

        return HashTable<T, Key, Hash>::endIter_;
    }
    else
    {
        return iterator(*this, table_[i], 0, i);
    }
}


template<class T, class Key, class Hash>
inline const typename HashTable<T, Key, Hash>::iterator&
HashTable<T, Key, Hash>::end()
{
    return HashTable<T, Key, Hash>::endIter_;
}


template<class T, class Key, class Hash>
inline typename HashTable<T, Key, Hash>::const_iterator
HashTable<T, Key, Hash>::begin() const
{
    label i = 0;

    while (table_ && !table_[i] && ++i < tableSize_);

    if (i == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable is empty\n";
        }
#       endif

        return HashTable<T, Key, Hash>::endConstIter_;
    }
    else
    {
        return const_iterator(*this, table_[i], 0, i);
    }
}


template<class T, class Key, class Hash>
inline const typename HashTable<T, Key, Hash>::const_iterator&
HashTable<T, Key, Hash>::end() const
{
    return HashTable<T, Key, Hash>::endConstIter_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
