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
#include "dimensionSet.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "IBoundaryTypeDescriptorImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IBoundaryTypeDescriptorImpl::IBoundaryTypeDescriptorImpl
(
    const word& boundaryTypeName
)
:
    name_(boundaryTypeName),
    displayName_(boundaryTypeName),
    description_(boundaryTypeName),
    patchType_("patch"),
    superType_("")
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::IBoundaryTypeDescriptorImpl"
        "(const char* boundaryTypeName)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IBoundaryTypeDescriptorImpl::IBoundaryTypeDescriptorImpl
(
    const word& boundaryTypeName,
    const wordList& fieldList
)
:
    name_(boundaryTypeName),
    displayName_(boundaryTypeName),
    description_(boundaryTypeName + " Boundary Condition"),
    patchType_(boundaryTypeName),
    superType_("")
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::IBoundaryTypeDescriptorImpl"
        "(const char* boundaryTypeName, const wordList& fieldList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Get list of fields.
        patchFieldTypes_.length(fieldList.size());

        // Loop over all defined fields and determine the patch field types.
        for (int nField = 0; nField < fieldList.size(); nField++)
        {
            patchFieldTypes_[nField].name = fieldList[nField].c_str();
            patchFieldTypes_[nField].value = boundaryTypeName.c_str();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IBoundaryTypeDescriptorImpl::IBoundaryTypeDescriptorImpl
(
    const word& boundaryTypeName,
    const Foam::dictionary& boundaryTypeDict,
    const wordList& fieldList
)
:
    name_(boundaryTypeName),
    displayName_(boundaryTypeName),
    description_(boundaryTypeName)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::IBoundaryTypeDescriptorImpl"
        "(const char* boundaryTypeName, Foam::dictionary& boundaryTypeDict, "
        "const wordList& fieldList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check that the dictionary has the minimum required information.
        if 
        (
            !boundaryTypeDict.found("displayName")
         || !boundaryTypeDict.found("description")
         || !boundaryTypeDict.found("patchType")
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Malformed boundary type definition dictionary for "
                "boundary type '" + name_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the boundary type properties.
        boundaryTypeDict.lookup("displayName") >> displayName_;
        boundaryTypeDict.lookup("description") >> description_;
        boundaryTypeDict.lookup("patchType"   ) >> patchType_;

        // Get (optional) super type.
        if (boundaryTypeDict.found("superType"))
        {
            boundaryTypeDict.lookup("superType") >> superType_;
        }

        // Get list of fields.
        patchFieldTypes_.length(fieldList.size());

        // Loop over all defined fields and determine the patch field types.
        for (int nField = 0; nField < fieldList.size(); nField++)
        {
            // Make sure there is an entry for this field.
            if (!boundaryTypeDict.found(fieldList[nField]))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Malformed boundary type definition dictionary for "
                   + name_ + ".",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            patchFieldTypes_[nField].name = fieldList[nField].c_str();
            patchFieldTypes_[nField].value =
                word(boundaryTypeDict.lookup(word(fieldList[nField]))).c_str();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IBoundaryTypeDescriptorImpl::~IBoundaryTypeDescriptorImpl()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::~IBoundaryTypeDescriptorImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IBoundaryTypeDescriptorImpl::name()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(name_.c_str());
}

void FoamX::IBoundaryTypeDescriptorImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IBoundaryTypeDescriptorImpl::displayName()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::displayName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(displayName_.c_str());
}

void FoamX::IBoundaryTypeDescriptorImpl::displayName(const char* newName)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::displayName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    displayName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IBoundaryTypeDescriptorImpl::description()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(description_.c_str());
}

void FoamX::IBoundaryTypeDescriptorImpl::description(const char* newDescription)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::description"
        "(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IBoundaryTypeDescriptorImpl::patchType()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::patchType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(patchType_.c_str());
}

void FoamX::IBoundaryTypeDescriptorImpl::patchType(const char* newName)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::patchType(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    patchType_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IBoundaryTypeDescriptorImpl::superType()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::superType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(superType_.c_str());
}

void FoamX::IBoundaryTypeDescriptorImpl::superType(const char* newName)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::superType(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    superType_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringPairList*
FoamX::IBoundaryTypeDescriptorImpl::patchFieldTypes()
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::patchFieldTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new StringPairList(patchFieldTypes_);
}

void FoamX::IBoundaryTypeDescriptorImpl::patchFieldTypes
(
    const FoamXServer::StringPairList& newList
)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::patchFieldTypes";

    LogEntry log(functionName, __FILE__, __LINE__);

    patchFieldTypes_ = newList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IBoundaryTypeDescriptorImpl::save
(
    DictionaryWriter& dictWriter
)
{
    static const char* functionName =
        "FoamX::IBoundaryTypeDescriptorImpl::save"
        "(DictionaryWriter& dictWriter)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Save the boundary type information to a sub-dictionary.
        dictWriter.startSubDict(name_);

        dictWriter.writeEntry("displayName", displayName_);
        dictWriter.writeEntry("description", description_);
        dictWriter.writeEntry("patchType", patchType_);

        // Only write the super type if it has been set.
        if (superType_.size()> 0)
        {
            dictWriter.writeEntry("superType", superType_);
        }

        // Write the patch field types.
        for (unsigned int i = 0; i <patchFieldTypes_.length(); i++)
        {
            dictWriter.writeEntry
            (
                word(patchFieldTypes_[i].name),
                word(patchFieldTypes_[i].value)
            );
        }

        dictWriter.endSubDict();
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
