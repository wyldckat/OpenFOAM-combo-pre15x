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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Reads the header information of a File up to and including the class name.

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "dictionary.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read header for given object type
bool IOobject::readHeader(Istream& is)
{
    if (IOobject::debug)
    {
        Info<< "IOobject::readHeader(Istream&) : reading header for file "
            << is.name() << endl;
    }

    // Check Istream not already bad
    if (!is.good())
    {
        if (rOpt_ == MUST_READ)
        {
            FatalIOErrorIn("IOobject::readHeader(Istream&)", is)
                << " stream not open for reading essential object"
                << exit(FatalIOError);
        }

        if (IOobject::debug)
        {
            SeriousIOErrorIn("IOobject::readHeader(Istream&)", is)
                << "IOobject::readHeader(Istream&) :"
                << " stream not open for reading from file "
                << is.name() << endl;
        }

        return false;
    }

    dictionary headerDict;

    token firstToken(is);

    if (is.good() && firstToken.isWord())
    {
        if (firstToken.wordToken() == "FoamFile")
        {
            headerDict = dictionary(is);
        }
        else if
        (
            firstToken.wordToken() == "version" 
         || firstToken.wordToken() == "format"
        )
        {
            is.putBack(firstToken);
            headerDict = dictionary(is, "object");
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }


    is.setVersion(headerDict.lookup("version"));
    is.setFormat(headerDict.lookup("format"));

    fileName headerRoot(headerDict.lookup("root"));
    if (IOobject::debug && headerRoot != rootPath())
    {
        IOWarningIn("IOobject::readHeader(Istream&)", is)
            << " root moved from " << headerRoot
            << " to " << rootPath()
            << " for file " << is.name() << endl;
    }

    fileName headerCase(headerDict.lookup("case"));
    if (IOobject::debug && headerCase != caseName())
    {
        IOWarningIn("IOobject::readHeader(Istream&)", is)
            << " case moved from " << headerCase
            << " to " << caseName()
            << " for file " << is.name() << endl;
    }

    // For backward compatibility, the local and instance paths and form
    // are optional.
    if (headerDict.found("instance"))
    {
        fileName headerInstance(headerDict.lookup("instance"));

        if (IOobject::debug && headerInstance != instance())
        {
            IOWarningIn("IOobject::readHeader(Istream&)", is)
                << " instance moved from " << headerInstance
                << " to " << instance()
                << " for file " << is.name() << endl;
        }
    }

    if (headerDict.found("local"))
    {
        fileName headerLocal(headerDict.lookup("local"));

        if (IOobject::debug && headerLocal != local())
        {
            IOWarningIn("IOobject::readHeader(Istream&)", is)
                << " local moved from " << headerLocal
                << " to " << local()
                << " for file " << is.name() << endl;
        }
    }

    if (headerDict.found("form"))
    {
        if (word(headerDict.lookup("form")) == "dictionary")
        {
            is.setVersion(2.0);
        }
    }

    // "object" should be the last entry in the header, if it isn't
    // then "class" might be as this was used in older versions of FOAM files.
    // To handle this case a special read opperation is used here to pickup
    // the class entry.  If this fails we have to give up.
    if (headerDict.found("class"))
    {
        headerClassName_ = word(headerDict.lookup("class"));
    }
    else
    {
        word keyClass(is);
        if (keyClass == "class")
        {
            char dummySemicolon;
            is >> headerClassName_ >> dummySemicolon;

            if (dummySemicolon != token::END_STATEMENT)
            {
                FatalIOErrorIn("IOobject::readHeader(Istream&)", is)
                    << " illformed class entry in header"
                    << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorIn("IOobject::readHeader(Istream&)", is)
                << " class is undefined"
                << exit(FatalIOError);
        }
    }

    // The note entry is optional
    if (headerDict.found("note"))
    {
        note_ = string(headerDict.lookup("note"));
    }

    word headerObject(headerDict.lookup("object"));
    if (IOobject::debug && headerObject != name())
    {
        IOWarningIn("IOobject::readHeader(Istream&)", is)
            << " object renamed from "
            << name() << " to " << headerObject
            << " for file " << is.name() << endl;
    }


    // Check stream is still OK
    if (!is.good())
    {
        if (rOpt_ == MUST_READ)
        {
            FatalIOErrorIn("IOobject::readHeader(Istream&)", is)
                << " stream failure while reading header for essential object"
                << name()
                << exit(FatalIOError);
        }

        if (IOobject::debug)
        {
            Info<< "IOobject::readHeader(Istream&) :"
                << " stream failure while reading header"
                << " on line " << is.lineNumber()
                << " of file " << is.name() << endl;
        }

        objState_ = BAD;

        return false;
    }

    objState_ = GOOD;

    if (IOobject::debug)
    {
        Info<< " .... read" << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
