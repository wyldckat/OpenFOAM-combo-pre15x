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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectRegistry, 0);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// Add an regIOobject to list
bool objectRegistry::checkIn(regIOobject& io) const
{
    if
    (
        this == dynamic_cast<const objectRegistry*>(&time_)
     && io.name() != "time"
     && io.name() != "region0"
    )
    {
        Warning
            << "objectRegistry::checkIn(regIOobject&): \n"
            << "    Registering object '" << io.name() << '\''
            << " with objectRegistry 'time'" << nl
            << "    This is only appropriate for registering the regions"
               " of a multi-region computation\n"
               "    or in other special circumstances.\n"
               "    Otherwise please register this object with it's region"
               " (mesh)"
            << endl;
    }

    if (objectRegistry::debug)
    {
        Info<< "objectRegistry::checkIn(regIOobject&) : "
            << "checking in " << io.name()
            << endl;
    }

    return ((regIOobjectTable&)(*this)).insert(io.name(), &io);
}


// Remove an regIOobject from list
bool objectRegistry::checkOut(regIOobject& io) const
{
    iterator iter = ((regIOobjectTable&)(*this)).find(io.name());

    if (iter != end())
    {
        if (objectRegistry::debug)
        {
            Info<< "objectRegistry::checkOut(regIOobject&) : "
                << "checking out " << io.name()
                << endl;
        }

        if (iter() != &io)
        {
            if (objectRegistry::debug)
            {
                Warning
                    << "objectRegistry::checkOut(regIOobject&) : "
                    << "attempt to checkOut copy of " << io.name()
                    << endl;
            }

            return false;
        }
        else
        {
            if (io.registries())
            {
                delete iter();
            }

            return ((regIOobjectTable&)(*this)).erase(iter);
        }
    }
    else
    {
        if (objectRegistry::debug)
        {
            Info<< "objectRegistry::checkOut(regIOobject&) : "
                << "could not find " << io.name()
                << endl;
        }

        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors *  * * * * * * * * * * * * * //

objectRegistry::objectRegistry
(
    const Time& t,
    const label nIoObjects
)
:
    regIOobject
    (
        IOobject
        (
            "",
            "",
            t,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        )
    ),
    regIOobjectTable(nIoObjects),
    time_(t),
    parent_(t),
    dbDir_(name())
{}


objectRegistry::objectRegistry
(
    const IOobject& io,
    const label nIoObjects
)
:
    regIOobject(io),
    regIOobjectTable(nIoObjects),
    time_(io.time()),
    parent_(io.db()),
    dbDir_(parent_.dbDir()/name())
{
    writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

objectRegistry::~objectRegistry()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->registries())
        {
            delete iter();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return the list of names of the IOobjects
wordList objectRegistry::names() const
{
    wordList objectNames(size());

    label count=0;
    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        objectNames[count++] = iter()->name();
    }

    return objectNames;
}


// Return the list of names of the IOobjects of given class
wordList objectRegistry::names(const word& ClassName) const
{
    wordList objectNames(size());

    label count=0;
    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->type() == ClassName)
        {
            objectNames[count++] = iter()->name();
        }
    }

    objectNames.setSize(count);

    return objectNames;
}


const objectRegistry& objectRegistry::subRegistry(const word& name) const
{
    return lookupObject<objectRegistry>(name);
}


bool objectRegistry::update()
{
    if (!updated())
    {
        bool ok = false;

        for (iterator iter = begin(); iter != end(); ++iter)
        {
            if (objectRegistry::debug)
            {
                Info<< "objectRegistry::update() : "
                    << "Considering updating object "
                    << iter()->name()
                    << endl;
            }

            ok = iter()->update() && ok;
        }

        regIOobject::update();

        return ok;
    }

    return true;
}


void objectRegistry::resetUpdate()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        iter()->resetUpdate();
    }

    regIOobject::resetUpdate();
}


void objectRegistry::readModifiedObjects()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Info<< "objectRegistry::readModifiedObjects() : "
                << "Considering reading object "
                << iter()->name()
                << endl;
        }

        iter()->readIfModified();
    }
}


bool objectRegistry::readIfModified()
{
    readModifiedObjects();
    return true;
}


bool objectRegistry::write
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
) const
{
    bool ok = true;

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Info<< "objectRegistry::write() : "
                << "Considering writing object "
                << iter()->name()
                << endl;
        }

        if (iter()->writeOpt() != NO_WRITE)
        {
            ok = iter()->write() && ok;
        }
    }

    return ok;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
