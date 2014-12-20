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

#include "dictionary.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

defineTypeNameAndDebug(dictionary, 0);

//- Null dictionary
const dictionary dictionary::null;


// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

// Add a new entry
void dictionary::add(entry* ePtr)
{
    if (!hashedEntries_.found(ePtr->keyword()))
    {
        append(ePtr);
        hashedEntries_.insert(ePtr->keyword(), ePtr);
    }
    else
    {
        Warning
            << "dictionary::add(entry* ePtr) : "
               "attempt to add an entry already in dictionary " << name()
            << endl;

        delete ePtr;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
dictionary::dictionary()
{}


// Construct as copy
dictionary::dictionary(const dictionary& dict)
:
    IDLList<entry>(dict),
    name_(dict.name())
{
    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert((*iter).keyword(), &(*iter));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dictionary::~dictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    
// Return line number of first token in dictionary
label dictionary::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}


// Return line number of last token in dictionary
label dictionary::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


// Find and return entry
bool dictionary::found(const word& keyword) const
{
    return hashedEntries_.found(keyword);
}


// Find and return an entry
const entry& dictionary::lookupEntry(const word& keyword) const
{
    HashTable<entry*>::const_iterator iter = hashedEntries_.find(keyword);

    if (iter == hashedEntries_.end())
    {
        // If keyword not found print error message ...
        FatalIOErrorIn
        (
            "dictionary::lookupEntry(const word& keyword) const",
            *this
        )   << " keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return *(*iter);
}


// Find and return an entry
ITstream& dictionary::lookup(const word& keyword) const
{
    return lookupEntry(keyword).stream();
}


// Find and return a sub-dictionary
const dictionary& dictionary::subDict(const word& keyword) const
{
    return lookupEntry(keyword).dict();
}


// Return the table of contents
wordList dictionary::toc() const
{
    wordList keywords(size());

    label i = 0;
    for
    (
        IDLList<entry>::const_iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        keywords[i++] = iter().keyword();
    }

    return keywords;
}


// Add an entry
void dictionary::add(const entry& e)
{
    add(e.clone().ptr());
}


// Add a token entry
void dictionary::add(const word& keyword, const token& t)
{
    add(new primitiveEntry(keyword, t));
}

// Add a word entry
void dictionary::add(const word& keyword, const word& w)
{
    add(new primitiveEntry(keyword, token(w)));
}

// Add a string entry
void dictionary::add(const word& keyword, const Foam::string& s)
{
    add(new primitiveEntry(keyword, token(s)));
}

// Add a label entry
void dictionary::add(const word& keyword, const label l)
{
    add(new primitiveEntry(keyword, token(l)));
}

// Add a word entry
void dictionary::add(const word& keyword, const scalar s)
{
    add(new primitiveEntry(keyword, token(s)));
}

// Add an entry constructed from ITstream
void dictionary::add(const word& keyword, const ITstream& tokens)
{
    add(new primitiveEntry(keyword, tokens));
}

// Add an entry constructed from tokenList
void dictionary::add(const word& keyword, const tokenList& tokens)
{
    add(new primitiveEntry(keyword, tokens));
}

// Add a dictionary entry
void dictionary::add(const word& keyword, const dictionary& dict)
{
    add(new dictionaryEntry(keyword, dict));
}



bool dictionary::remove(const word& Keyword)
{
    HashTable<entry*>::iterator iter = hashedEntries_.find(Keyword);

    if (iter != hashedEntries_.end())
    {
        IDLList<entry>::remove(iter());
        delete iter();
        hashedEntries_.erase(iter);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Find and return entry
ITstream& dictionary::operator[](const word& keyword) const
{
    return lookup(keyword);
}


void dictionary::operator=(const dictionary& dict)
{
    // Check for assignment to self
    if (this == &dict)
    {
        FatalErrorIn("dictionary::operator=(const dictionary&)")
            << "attempted assignment to self for dictionary " << name()
            << abort(FatalError);
    }

    IDLList<entry>::operator=(dict);
    name_ = dict.name();

    hashedEntries_.clear();

    for
    (
        IDLList<entry>::iterator iter = begin();
        iter != end();
        ++iter
    )
    {
        hashedEntries_.insert((*iter).keyword(), &(*iter));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
