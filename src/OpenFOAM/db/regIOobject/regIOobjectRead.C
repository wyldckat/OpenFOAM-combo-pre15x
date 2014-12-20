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

#include "regIOobject.H"
#include "IFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read header and return stream
Istream& regIOobject::readStream()
{
    if (IFstream::debug)
    {
        Info<< "regIOobject::readStream() : "
            << "reading object " << name()
            << " from file " << objectPath()
            << endl;
    }

    if (readOpt() == NO_READ)
    {
        FatalErrorIn("regIOobject::readStream(const word&)")
            << "NO_READ specified for read-constructor of object " << name()
            << " of class " << headerClassName()
            << abort(FatalError);
    }

    // Construct object stream and read header if not already constructed
    if (!isPtr_)
    {
        if (!(isPtr_ = objectStream()))
        {
            FatalIOError
            (
                "regIOobject::readStream(const word&)",
                __FILE__,
                __LINE__,
                objectPath(),
                0
            )   << "cannot open file"
                << exit(FatalIOError);
        }
        else if (!readHeader(*isPtr_))
        {
            FatalIOErrorIn("regIOobject::readStream(const word&)", *isPtr_)
                << "problem while reading header for object " << name()
                << exit(FatalIOError);
        }
    }

    if (!lastModified_)
    {
        lastModified_ = lastModified(objectPath());
    }

    return *isPtr_;
}


// Read header for given object type
Istream& regIOobject::readStream(const word& aType)
{
    if (IFstream::debug)
    {
        Info<< "regIOobject::readStream(const word&) : "
            << "reading object " << name()
            << " from file " << objectPath()
            << endl;
    }

    // Construct IFstream if not already constructed
    if (!isPtr_)
    {
        readStream();

        // Check the className of the regIOobject
        // dictionary is an allowable name in case the actual class
        // instantiated is a dictionary
        if
        (
            headerClassName() != aType
         && headerClassName() != "dictionary"
        )
        {
            FatalIOErrorIn("regIOobject::readStream(const word&)", *isPtr_)
                << "unexpected class name " << headerClassName()
                << " expected " << aType << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }
    }

    return *isPtr_;
}


// Close Istream
void regIOobject::close()
{
    if (IFstream::debug)
    {
        Info<< "regIOobject::close() : "
            << "finished reading " << objectPath()
            << endl;
    }

    if (isPtr_)
    {
        delete isPtr_;
        isPtr_ = NULL;
    }
}


bool regIOobject::readData(Istream&)
{
    return false;
}


bool regIOobject::read()
{
    bool ok = readData(readStream(type()));
    close();
    return ok;
}


bool regIOobject::readIfModified()
{
    bool modified = false;

    if (lastModified_)
    {
        time_t newTimeStamp = lastModified(objectPath());

        if (newTimeStamp > lastModified_ + fileModificationSkew)
        {
            lastModified_ = newTimeStamp;
            modified = read();
        }
        else
        {
            modified = false;
        }
    }

    return modified;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
