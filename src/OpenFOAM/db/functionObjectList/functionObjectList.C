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

\*----------------------------------------------------------------------------*/

#include "functionObjectList.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList(const Time& t)
:
    HashPtrTable<functionObject>(),
    time_(t)
{
    if (t.controlDict().found("functions"))
    {
        HashPtrTable<functionObject> functions
        (
            t.controlDict().lookup("functions"),
            functionObject::iNew(time_)
        );
        
        transfer(functions);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectList::~functionObjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjectList::start()
{
    bool ok = false;

    forAllIter(HashPtrTable<functionObject>, *this, iter)
    {
        ok = iter()->start() && ok;
    }

    return ok;
}


bool Foam::functionObjectList::execute()
{
    bool ok = false;

    forAllIter(HashPtrTable<functionObject>, *this, iter)
    {
        ok = iter()->execute() && ok;
    }

    return ok;
}


bool Foam::functionObjectList::read()
{
    bool read = false;

    if (time_.controlDict().found("functions"))
    {
        HashPtrTable<dictionary> functionDicts
        (
            time_.controlDict().lookup("functions")
        );

        // Update existing and add new functionObjects
        forAllConstIter(HashPtrTable<dictionary>, functionDicts, iter)
        {
            if (found(iter.key()))
            {
                read = find(iter.key())()->read(*iter()) && read;
            }
            else
            {
                functionObject* functionObjectPtr = 
                    functionObject::New(time_, *iter()).ptr();

                functionObjectPtr->start();

                insert(iter.key(), functionObjectPtr);
            }
        }

        // Remove deleted functionObjects
        forAllIter(HashPtrTable<functionObject>, *this, iter)
        {
            if (!functionDicts.found(iter.key()))
            {
                erase(iter);
            }
        }
    }
    else
    {
        clear();
        read = true;
    }

    return read;
}


// ************************************************************************* //
