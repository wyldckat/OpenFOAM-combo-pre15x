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
// Foam header files.
#include "word.H"
#include "string.H"
#include "IFstream.H"
#include "wordList.H"
#include "stringList.H"
#include "fileNameList.H"
#include "OFstream.H"
#include "dimensionSet.H"
#include "Time.H"
#include "banner.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "DictionaryWriter.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::DictionaryWriter::DictionaryWriter(const fileName& fileName)
:
    pFileStream_(NULL),
    fileName_(fileName),
    caseRoot_(""),
    caseName_(""),
    instance_(fileName_.path()),
    object_(fileName_.name()),
    entryNameSize_(20)
{
    static const char* functionName =
        "FoamX::DictionaryWriter::DictionaryWriter(const fileName& fileName)";

    try
    {
        // Make directory if required.
        if (!dir(fileName_.path()) && !mkDir(fileName_.path()))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to create directory '" + fileName_.path() + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Backup existing file. TODO : Make this behaviour configurable.
        if (Foam::file(fileName_))
        {
            mv(fileName_, fileName_ + ".bak");
        }

        // Open the file stream.
        pFileStream_ = new OFstream(fileName_);

        if (!pFileStream_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to open file '" + fileName_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::DictionaryWriter::DictionaryWriter
(
    const fileName& caseRoot,
    const fileName& caseName,
    const fileName& instance,
    const fileName& object
)
:
    pFileStream_(NULL),
    caseRoot_(caseRoot),
    caseName_(caseName),
    instance_(instance),
    object_(object),
    entryNameSize_(20)
{
    static const char* functionName =
        "FoamX::DictionaryWriter::DictionaryWriter"
        "(const fileName& caseRoot, const fileName& caseName, "
        "const fileName& instance, const fileName& object)";

    try
    {
        // Make directory if required.
        fileName_ = caseRoot_ / caseName_ / instance_ / object_;
        if (!dir(fileName_.path()) && !mkDir(fileName_.path()))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to create directory '" + fileName_.path() + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Backup existing file. TODO : Make this behaviour configurable.
        if (Foam::file(fileName_))
        {
            mv(fileName_, fileName_ + ".bak");
        }

        // Open the file stream.
        pFileStream_ = new OFstream(fileName_);

        if (!pFileStream_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to open file '" + fileName_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::DictionaryWriter::~DictionaryWriter()
{
    if (pFileStream_ != NULL)
    {
        delete pFileStream_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeHeader
(
    const Foam::string& title,
    const word& className
)
{
    writeBanner(file());
    file() << endl;
    writeComment(title);
    writeEndl();
    startSubDict("FoamFile");
    writeEntry("version", IOstream::currentVersion.str().c_str());
    writeEntry("format", word("ascii"));
    writeEndl();
    writeEntry("root", caseRoot_);
    writeEntry("case", caseName_);
    writeEntry("instance", instance_);
    writeEntry("local", Foam::string(""));
    writeEndl();
    writeEntry("class", className);
    writeEntry("object", fileName_.name());
    endSubDict();
    writeEndl();
    writeBar();
    writeEndl();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeSectionHeader(const Foam::string& name)
{
    file() << endl << "// " << name.c_str() << endl << "// ";
    writeChars('~', name.size());
    writeEndl();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeKeyword(const word& name)
{
    file().writeKeyword(name);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeKeywordOnly(const word& name)
{
    file().indent();
    file() << name;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const char* value
)
{
    writeKeyword(name);
    file() << value;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const FoamXString& value
)
{
    writeKeyword(name);
    value.write(file());
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const Foam::string& value
)
{
    writeKeyword(name);
    file() << value;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry(const word& name, const word& value)
{
    writeKeyword(name);
    file() << value;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry(const word& name, const label& value)
{
    writeKeyword(name);
    file() << value;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const Foam::word& name,
    const Foam::scalar& value
)
{
    writeKeyword(name);
    file() << value;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const Foam::word& name,
    const bool& value
)
{
    writeKeyword(name);
    file() << value;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const Foam::word& name,
    const FoamX::FoamXAny& value
)
{
    writeKeyword(name);
    writeValue(value);
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const FoamXServer::DimensionSet& entry
)
{
    writeKeyword(name);
    file() << token::BEGIN_SQR
           <<        entry.mass
           << " " << entry.length
           << " " << entry.time
           << " " << entry.temperature
           << " " << entry.moles
           << " " << entry.current
           << " " << entry.luminousIntensity
           << token::END_SQR;
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::startVectorSpace()
{
    file() << token::BEGIN_LIST;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::endVectorSpace()
{
    file() << token::END_LIST;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::startList(const Foam::label& numElements)
{
    file() << nl;

    if (numElements > 100)
    {
        file() << Foam::indent << numElements << nl;
    }

    file() << Foam::indent << token::BEGIN_LIST << incrIndent;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::endList()
{
    file() << nl << decrIndent << Foam::indent << token::END_LIST;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const FoamXServer::StringList& list
)
{
    writeKeyword(name);

    startList(list.length());

    for (unsigned int i = 0; i <list.length(); i++)
    {
        file() << Foam::indent << Foam::string(list[i]) << nl;
    }

    endList();
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const FoamX::FoamXStringList& list
)
{
    writeKeyword(name);

    startList(list.length());

    for (unsigned int i = 0; i <list.length(); i++)
    {
        file() << Foam::indent << Foam::string(list[i]) << nl;
    }

    endList();
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const word& name,
    const FoamX::FoamXWordList& list
)
{
    writeKeyword(name);

    startList(list.length());

    for (unsigned int i = 0; i <list.length(); i++)
    {
        file() << Foam::indent << word(list[i]) << nl;
    }

    endList();
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEntry
(
    const Foam::word& name,
    const FoamX::FoamXAnyList& list
)
{
    writeKeyword(name);

    startList(list.size());

    forAll(list, i)
    {
        file().indent();
        list[i].write(file());
        file() << nl;
    }

    endList();
    endEntry();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::startSubDict()
{
    file() << endl << Foam::indent << token::BEGIN_BLOCK << endl << incrIndent;
}

void FoamX::DictionaryWriter::startSubDict(const word& name)
{
    file() << Foam::indent << name;
    startSubDict();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::endDict()
{
    file() << decrIndent << Foam::indent << token::END_BLOCK;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::endSubDict()
{
    endDict();
    file() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeString(const Foam::string& text)
{
    file() << text.c_str();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeLine(const Foam::string& text)
{
    file() << Foam::indent << text.c_str() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeComment(const Foam::string& text)
{
    file() << Foam::indent << "// " << text.c_str() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeBar()
{
    writeLine
    (
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEndBar()
{
    writeLine
    (
        "\n// ************************************************************************* //"
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeEndl()
{
    file() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::endEntry()
{
    file() << token::END_STATEMENT << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::indent()
{
    file().indent();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::DictionaryWriter::writeChars(char ch, int count)
{
    for (register int i = 0; i < count; i++)
    {
        file() << ch;
    }
}


// ************************************************************************* //
