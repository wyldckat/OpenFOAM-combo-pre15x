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
    Constructors & destructor for regIOobject.

\*---------------------------------------------------------------------------*/

#include "regIOobject.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(regIOobject, 0);

int regIOobject::fileModificationSkew
(
    debug::optimisationSwitch("fileModificationSkew", 30)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOobject
regIOobject::regIOobject(const IOobject& io)
:
    IOobject(io),
    registered_(false),
    registries_(false),
    lastModified_(0),
    updated_(false),
    isPtr_(NULL)
{
    // Register with objectRegistry if requested
    if (registerObject())
    {
        checkIn();
    }
}


// Construct as copy
regIOobject::regIOobject(const regIOobject& rio)
:
    IOobject(rio),
    registered_(false),
    registries_(false),
    lastModified_(rio.lastModified_),
    updated_(false),
    isPtr_(NULL)
{
    // Do not register copy with objectRegistry
}


// Construct as copy, and transfering objectRegistry registration to copy
// if registerCopy is true
regIOobject::regIOobject(const regIOobject& rio, bool registerCopy)
:
    IOobject(rio),
    registered_(false),
    registries_(false),
    lastModified_(rio.lastModified_),
    updated_(false),
    isPtr_(NULL)
{
    if (registerCopy && rio.registered_)
    {
        ((regIOobject&)rio).checkOut();
        checkIn();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Delete read stream, checkout from objectRegistry and destroy
regIOobject::~regIOobject()
{
    if (objectRegistry::debug)
    {
        Info<< "Destroying regIOobject called " << name()
            << " of type " << type()
            << " in directory " << path()
            << endl;
    }

    if (isPtr_)
    {
        delete isPtr_;
    }

    // Check out of objectRegistry
    checkOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void regIOobject::checkIn()
{
    if (!registered_)
    {
        // Attempt to register object with objectRegistry
        if (!db().checkIn(*this))
        {
            if (objectRegistry::debug)
            {
                Warning
                    << "regIOobject::checkIn() : "
                    << "failed to register object " << name()
                    << " the name already exists in the objectRegistry"
                    << endl;
            }
        }
        else
        {
            registered_ = true;
        }
    }
}


void regIOobject::checkOut()
{
    if (registered_)
    {
        db().checkOut(*this);
        registered_ = false;
    }
}


// Rename object and re-register with objectRegistry under new name
void regIOobject::rename(const word& newName)
{
    // Check out of objectRegistry
    checkOut();

    IOobject::rename(newName);

    // Re-register object with objectRegistry
    checkIn();
}


// Reset the updated flag
void regIOobject::resetUpdate()
{
    updated_ = false;
}


// Update this object if it hasn't already been updated
bool regIOobject::update()
{
    updated_ = true;

    return true;
}


// Assign to IOobject
void regIOobject::operator=(const IOobject& io)
{
    if (isPtr_)
    {
        delete isPtr_;
        isPtr_ = NULL;
    }

    // Check out of objectRegistry
    checkOut();

    IOobject::operator=(io);

    // Re-register object with objectRegistry
    checkIn();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
