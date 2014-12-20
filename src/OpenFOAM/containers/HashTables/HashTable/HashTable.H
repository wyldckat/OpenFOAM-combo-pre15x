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

SourceFiles
    HashTableI.H
    HashTable.C
    HashTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashTable_H
#define HashTable_H

#include "word.H"
#include "hashString.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class T>
class List;

template<class T, class Key, class Hash> class HashTable;

template<class T, class Key, class Hash> Istream& operator>>
(
    Istream&,
    HashTable<T, Key, Hash>&
);

template<class T, class Key, class Hash> Ostream& operator<<
(
    Ostream&,
    const HashTable<T, Key, Hash>&
);


/*---------------------------------------------------------------------------*\
                        Class HashTableName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(HashTable);


/*---------------------------------------------------------------------------*\
                          Class HashTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=hashString>
class HashTable
:
    public HashTableName
{
    // Private data type for table entries

        struct hashedEntry
        {
            //- The lookup key
            Key key_;

            //- Pointer to next hashedEntry in sub-list
            hashedEntry* next_;

            //- The data object
            T obj_;

            //- Constructors

                //- Construct given key, next pointer and object
                inline hashedEntry
                (
                    const Key& key,
                    hashedEntry* next,
                    T newEntry
                );

                //- Dissallow construction as copy
                hashedEntry(const hashedEntry&);
        };


    // Private data: size of table, the table and current number of elements

        //- Number of primary entries allocated in table (not necessarily used)
        label tableSize_;

        //- The table of primary entries
        hashedEntry** table_;

        //- The current number of elements in table
        label nElmts_;

        //- Hashing function class
        Hash hash_;


public:

        //- Declare friendship with the HashPtrTable class
        template<class T2, class Key2, class Hash2>
        friend class HashPtrTable;


    // Forward declaration of STL iterators

        template<class TRef, class TableRef, class HashedEntryPtr>
        class Iterator;

        typedef Iterator
        <
            T&,
            HashTable<T, Key, Hash>&,
            hashedEntry*
        > iterator;

        typedef Iterator
        <
            const T&,
            const HashTable<T, Key, Hash>&,
            const hashedEntry*
        > const_iterator;


    // Declare friendship with the iterators

        friend class Iterator
        <
            T&,
            HashTable<T, Key, Hash>&,
            hashedEntry*
        >;

        friend class Iterator
        <
            const T&,
            const HashTable<T, Key, Hash>&,
            const hashedEntry*
        >;


    // Constructors

        //- Construct given initial table size
        HashTable(const label size = 100);

        //- Construct from Istream
        HashTable(Istream&, const label size = 100);

        //- Construct as copy
        HashTable(const HashTable<T, Key, Hash>&);


    // Destructor

        ~HashTable();


    // Member Functions

        // Access

            //- Return number of elements in table.
            inline label size() const;

            //- Return true if hashedEntry is found in table
            bool found(const Key& key) const;

            //- Find and return an iterator set at the hashedEntry
            //  If not found iterator = end()
            iterator find(const Key& key);

            //- Find and return an const_iterator set at the hashedEntry
            //  If not found iterator = end()
            const_iterator find(const Key& key) const;

            //- Return the table of contents
            List<Key> toc() const;


        // Edit

            //- Insert a new hashedEntry
            bool insert(const Key& key, T newElmt);

            //- Erase an hashedEntry specified by given iterator
            bool erase(const iterator& it);

            //- Erase an hashedEntry specified by given key if in table
            bool erase(const Key& key);

            //- Resize the hash table for efficiency
            void resize(const label newSize);

            //- Clear all entries from table
            void clear();

            //- Transfer the contents of the argument table into this table
            //  and annull the argument table.
            void transfer(HashTable<T, Key, Hash>&);


    // Member Operators

        //- Find and return an hashedEntry
        inline T& operator[](const Key& key);

        //- Find and return an hashedEntry
        inline const T& operator[](const Key& key) const;

        //- Assignment
        void operator=(const HashTable<T, Key, Hash>&);


    // STL type definitions

        //- Type of values the HashTable contains.
        typedef T value_type;

        //- Type that can be used for storing into HashTable::value_type
        //  objects.  This type is usually List::value_type&.
        typedef T& reference;

        //- Type that can be used for storing into constant
        //  HashTable::value_type objects.  This type is usually const
        //  HashTable::value_type&.
        typedef const T& const_reference;

        //- The type that can represent the size of a HashTable.
        typedef label size_type;


    // STL iterator

        template<class TRef, class TableRef, class HashedEntryPtr>
        class Iterator
        {
            friend class HashTable;

            template<class TRef2, class TableRef2, class HashedEntryPtr2>
            friend class Iterator;

            // Private data

                //- Reference to the HashTable this is an iterator for
                TableRef curHashTable_;

                //- Current element
                HashedEntryPtr elmtPtr_;

                //- Previous element
                HashedEntryPtr prevElmtPtr_;

                //- Current hash index
                label hashIndex_;


        public:

            // Constructors

                //- Construct from hash table, element and hash index
                inline Iterator
                (
                    TableRef curHashTable,
                    HashedEntryPtr elmt,
                    HashedEntryPtr prev,
                    label hashIndex
                );

                //- Construct from the non-const iterator
                inline Iterator(const iterator&);


            // Member operators

                inline void operator=(const iterator& iter);

                inline bool operator==(const iterator& iter) const;
                inline bool operator==(const const_iterator& iter) const;

                inline bool operator!=(const iterator& iter) const;
                inline bool operator!=(const const_iterator& iter) const;

                inline TRef operator*();
                inline TRef operator()();

                inline Iterator& operator++();
                inline Iterator operator++(int);

                inline const Key& key();
        };


        //- iterator set to the begining of the HashTable
        inline iterator begin();

        //- iterator set to beyond the end of the HashTable
        inline const iterator& end();

        //- const_iterator set to the begining of the HashTable
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the HashTable
        inline const const_iterator& end() const;


    // IOstream Operator

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream&,
            HashTable<T, Key, Hash>&
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream&,
            const HashTable<T, Key, Hash>&
        );


private:

        //- iterator returned by end()
        iterator endIter_;

        //- const_iterator returned by end()
        const_iterator endConstIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "HashTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef NoHashTableC
#ifdef NoRepository
#   include "HashTable.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
