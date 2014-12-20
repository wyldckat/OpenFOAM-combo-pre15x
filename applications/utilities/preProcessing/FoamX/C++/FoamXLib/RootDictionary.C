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
#include "OSspecific.H"
#include "word.H"
#include "string.H"
#include "dictionary.H"
#include "wordList.H"
#include "stringList.H"
#include "dimensionSet.H"

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "LogEntry.H"
#include "RootDictionary.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::RootDictionary::RootDictionary
(
    ITypeDescriptor_ptr typeDesc,
    const Foam::fileName& caseRoot,
    const Foam::fileName& caseName
)
:
    IDictionaryEntryImpl(typeDesc),
    caseRoot_(caseRoot),
    caseName_(caseName)
{
    static const char* functionName = 
        "FoamX::RootDictionary::RootDictionary"
        "(ITypeDescriptor_ptr, const Foam::fileName&, const Foam::fileName&)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::RootDictionary::~RootDictionary()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::RootDictionary::load(Foam::Istream& is)
{
    static const char* functionName =
        "FoamX::RootDictionary::load(Foam::Istream&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Loading file for "
            << CORBA::String_var(typeDescriptor_->path()) << endl;

        token firstToken(is);

        if
        (
            is.good()
         && firstToken.isWord()
         && firstToken.wordToken() == "FoamFile"
        )
        {
            // Read the FoamFile header
            dictionary headerDict(is);
        }
        else
        {
            is.putBack(firstToken);
        }

        IDictionaryEntryImpl::load(is);

        // Check state of Istream.
        is.check(functionName);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::RootDictionary::save()
{
    static const char* functionName = "FoamX::RootDictionary::save()";

    // Overridden for root dictionary entries.
    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure that this entry is a top-level compound type.
        CORBA::String_var dictPath = typeDescriptor_->dictionaryPath();
        CORBA::String_var dictName = typeDescriptor_->name();

        if (!typeDescriptor_->isCompoundType() || strlen(dictName) == 0)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid root dictionary object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the dictionaries path and file name.
        Foam::fileName dictFilePath = caseRoot_/caseName_/fileName(dictPath);
        Foam::fileName dictFileName = dictFilePath/fileName(dictName);

        log << "Saving root dictionary " << dictFileName << "." << endl;

        DictionaryWriter dictWriter
        (
            caseRoot_,
            caseName_,
            (const char*)dictPath,
            (const char*)dictName
        );

        dictWriter.writeHeader
        (
            "FoamX Case Dictionary.",
            "dictionary"
        );

        // Write the comment if required.
        const char* comment = typeDescriptor_->comment();
        if (strlen(comment) > 0)
        {
            dictWriter.writeComment(comment);
            dictWriter.writeEndl();
        }

        bool isDict = true;

        if (typeDescriptor_->type() != Type_Dictionary)
        {
            isDict = false;
        }

        // Save all sub-elements.
        for
        (
            Foam::DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();
            iter != subElements_.end();
            ++iter
        )
        {
            // Write each sub-element as a dictionary entry.
            iter()->save(dictWriter, isDict);
            dictWriter.writeEndl();
        }

        dictWriter.writeEndBar();
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
