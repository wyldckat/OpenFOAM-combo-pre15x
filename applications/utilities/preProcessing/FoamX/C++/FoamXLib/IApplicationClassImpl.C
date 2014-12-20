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
#include "fvPatch.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "ITypeDescriptorImpl.H"
#include "IApplicationClassImpl.H"
#include "IPropertiesImpl.H"
#include "IGeometricFieldDescriptorImpl.H"
#include "IBoundaryTypeDescriptorImpl.H"
#include "DictionaryWriter.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IApplicationClassImpl::IApplicationClassImpl
(
    const FoamXServer::ApplicationClassDescriptor& appClassDesc,
    const IPropertiesImpl& foamSystemProperties
)
:
    name_(appClassDesc.name),
    description_(appClassDesc.name),
    category_(appClassDesc.category),
    systemClass_(appClassDesc.systemClass),
    foamSystemProperties_(foamSystemProperties)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::IApplicationClassImpl"
        "(const char* className, bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Set the default application class definition dictionary name.
    classDictName_ =
        (systemClass_ ? Paths::system : Paths::user)
       /"applications"/category_/name_/name_ + ".cfg";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IApplicationClassImpl::~IApplicationClassImpl()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::~IApplicationClassImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::validate()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::load()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::load()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // See if this application class has been overriden by the user.
        if
        (
            exists
            (
                Paths::user/"applications"/category_/name_/(name_ + ".cfg")
            )
        )
        {
            classDictName_ =
                Paths::user/"applications"/category_/name_/name_ + ".cfg";

            systemClass_ = false;
        }
        else
        {
            classDictName_ =
                Paths::system/"applications"/category_/name_/name_ + ".cfg";

            systemClass_ = true;
        }

        // Make sure that the definition dictionary exists somewhere.
        if (!exists(classDictName_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Application Class definition file '"
               + classDictName_ + "' could not be found.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Open the application class definition dictionary.
        dictionary classDict((IFstream(classDictName_)()));

        // Make sure we have all of the required entries.
        if
        (
            !classDict.found("description")
         || !classDict.found("modules")
         || !classDict.found("fields")
         || !classDict.found("boundaryTypes")
         || !classDict.found("dictionaries")
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid application class configuration dictionary '"
               + classDictName_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get application class info.
        classDict.lookup("description") >> description_;

        moduleList_.read(classDict.lookup("modules"));

        dictionaryList_.read(classDict.lookup("dictionaries"));

        // Initialise the dictionary objects.
        forAll(dictionaryList_, i)
        {
            word dictionaryName(dictionaryList_[i]);

            // See if we already have a dictionary with this name.
            if (dictTypeDescriptorMap_.found(dictionaryName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate dictionary name '" + dictionaryName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Check that dictionary has a definition dictionary and open it.
            fileName configDictName =
                classDictName_.path()/dictionaryName + ".cfg";

            if (!exists(configDictName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Dictionary definition file '" + configDictName.name()
                  + "' could not be found.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            dictionary configDict((IFstream(configDictName)()));

            log << "Reading dictionary type descriptor '" << configDict.name()
                << "'." << endl;

            // Create and initialise a new TypeDescriptor object for this
            // dictionary.
            ITypeDescriptorImpl* pTypeDescriptor = new ITypeDescriptorImpl
            (
                dictionaryName,
                "",
                configDict,
                foamSystemProperties_.foamTypesDict()
            );

            if (pTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create Dictionary TypeDescriptor object for "
                    "dictionary '" + dictionaryName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            dictTypeDescriptorMap_.insert(dictionaryName, pTypeDescriptor);
        }


        const dictionary& fieldsDict(classDict.subDict("fields"));

        // Initialise the field objects.
        for
        (
            dictionary::const_iterator iter = fieldsDict.begin();
            iter != fieldsDict.end();
            ++iter
        )
        {
            word fieldName(iter().keyword());

            // See if we already have a field with this name.
            if (fieldDescriptorMap_.found(fieldName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate field name '" + fieldName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new FieldDescriptor object for field.
            IGeometricFieldDescriptorImpl* pFieldDescriptor = new
            IGeometricFieldDescriptorImpl
            (
                fieldName,
                fieldsDict.subDict(fieldName),
                foamSystemProperties_.foamTypes(),
                foamSystemProperties_.geometryDescriptors()
            );

            if (pFieldDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create FieldDescriptor object for field '"
                   + fieldName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            fieldDescriptorMap_.insert(fieldName, pFieldDescriptor);
        }

        wordList fieldList = fieldDescriptorMap_.toc();


        wordList constranedPatchTypes = polyPatch::constraintTypes();

        forAll(constranedPatchTypes, i)
        {
            // See if we already have a boundary type with this name.
            if (boundaryTypeDescriptorMap_.found(constranedPatchTypes[i]))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate boundary type name '"
                  + constranedPatchTypes[i] + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new BoundaryTypeDescriptor object for
            // this field.
            IBoundaryTypeDescriptorImpl* pBoundaryDescriptor = new
            IBoundaryTypeDescriptorImpl
            (
                constranedPatchTypes[i],
                fieldList
            );

            if (pBoundaryDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create BoundaryTypeDescriptor object for "
                    "boundary type '" + constranedPatchTypes[i] + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            boundaryTypeDescriptorMap_.insert
            (
                constranedPatchTypes[i],
                pBoundaryDescriptor
            );
        }


        const dictionary& boundaryTypeDict(classDict.subDict("boundaryTypes"));

        // Initialise the boundary type information.
        for
        (
            dictionary::const_iterator iter = boundaryTypeDict.begin();
            iter != boundaryTypeDict.end();
            ++iter
        )
        {
            word boundaryTypeName(iter().keyword());

            // See if we already have a boundary type with this name.
            if (boundaryTypeDescriptorMap_.found(boundaryTypeName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate boundary type name '" + boundaryTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new BoundaryTypeDescriptor object for
            // this field.
            IBoundaryTypeDescriptorImpl* pBoundaryDescriptor = new
            IBoundaryTypeDescriptorImpl
            (
                boundaryTypeName,
                iter().dict(),
                fieldList
            );

            if (pBoundaryDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create BoundaryTypeDescriptor object for "
                    "boundary type '" + boundaryTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            boundaryTypeDescriptorMap_.insert
            (
                boundaryTypeName,
                pBoundaryDescriptor
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::save()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::save()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Cannot save system application classes.
        if (systemClass_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Cannot save system application class definition file '"
               + classDictName_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Write the application class data.
        DictionaryWriter dict(classDictName_);

        dict.writeHeader
        (
            "FoamX Application Class Configuration File for " + name_,
            "dictionary"
        );

        // Write execution info.
        dict.writeSectionHeader("Execution");
        dict.writeEntry("description", description_);
        dict.writeEndl();
        dict.writeEntry("arguments", arguments_);

        // Write module names.
        dict.writeSectionHeader("FoamX Modules");
        dict.writeEntry("modules", moduleList_);

        // Write dictionary names.
        dict.writeSectionHeader("Dictionaries");
        wordList dictNames;
        dictNames.setSize(dictTypeDescriptorMap_.size());
        label i = 0;
        for
        (
            HashTable<ITypeDescriptorImpl*>::iterator iter =
                dictTypeDescriptorMap_.begin() ;
            iter != dictTypeDescriptorMap_.end() ;
            ++iter
        )
        {
            CORBA::String_var name = iter()->name();
            dictNames[i++] = word(name);
        }

        // Use the top-level dictionary names instead of the keys.
        dict.writeEntry("dictionaries", dictNames);

        // Write field definitions.
        dict.writeSectionHeader("Fields");
        dict.startSubDict("fields");
        i = 0;

        for
        (
            HashTable<IGeometricFieldDescriptorImpl*>::iterator iter
                = fieldDescriptorMap_.begin();
            iter != fieldDescriptorMap_.end();
            ++iter
        )
        {
            iter()->save(dict);

            if (i++ < fieldDescriptorMap_.size()-1)
            {
                dict.writeEndl();
            }
        }
        dict.endSubDict();

        // Write boundary types.
        dict.writeSectionHeader("Boundary Types");
        dict.startSubDict("boundaryTypes");
        i = 0;

        for
        (
            HashTable<IBoundaryTypeDescriptorImpl*>::iterator iter
                = boundaryTypeDescriptorMap_.begin();
            iter != boundaryTypeDescriptorMap_.end();
            ++iter
        )
        {
            iter()->save(dict);

            if (i++ < boundaryTypeDescriptorMap_.size()-1)
            {
                dict.writeEndl();
            }
        }
        dict.endSubDict();

        dict.writeEndl();
        dict.writeEndBar();


        // Save the dictionary type definitions.
        forAll(dictionaryList_, i)
        {
            word dictKey(dictionaryList_[i]);

            ITypeDescriptorImpl* pDictType = dictTypeDescriptorMap_[dictKey];

            // Create a dictionary writer object for this dictionary.

            // Use the dictionary name, not the key.
            CORBA::String_var dictName = pDictType->name();

            // Same directory as master config file.
            fileName dictFilePath = classDictName_.path();

            fileName dictFileName = dictFilePath/word(dictName) + ".cfg";

            // Save dictionary type description.
            DictionaryWriter dictWriter(dictFileName);

            dictWriter.writeHeader
            (
                "FoamX Application Class Configuration File.",
                "dictionary"
            );

            // Top-level type descriptor.
            pDictType->save(dictWriter, true);

            dictWriter.writeEndl();
            dictWriter.writeEndBar();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IApplicationClassImpl::name()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(name_.c_str());
}

void FoamX::IApplicationClassImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IApplicationClassImpl::description()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(description_.c_str());
}

void FoamX::IApplicationClassImpl::description(const char* newDescription)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::description(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IApplicationClassImpl::category()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::category()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(category_.c_str());
}

void FoamX::IApplicationClassImpl::category(const char* newCategory)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::category(const char* newCategory)";

    LogEntry log(functionName, __FILE__, __LINE__);

    category_ = newCategory;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationClassImpl::modules()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::modules()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(moduleList_);
}

void FoamX::IApplicationClassImpl::modules
(
    const FoamXServer::StringList& newList
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::modules"
        "(const FoamXServer::StringList& newList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    moduleList_ = newList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IApplicationClassImpl::systemClass()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::systemClass()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return systemClass_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationClassImpl::fields()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::fields()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList
    (
        FoamXWordList((const wordList&)fieldDescriptorMap_.toc())
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::getField
(
    const char* fieldName,
    IGeometricFieldDescriptor_out fieldDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::getField"
        "(const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have this field's type descriptor cached.
        if (!fieldDescriptorMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the FieldDescriptor object.
        fieldDescriptor = fieldDescriptorMap_[fieldName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IApplicationClassImpl::findField
(
    const char* fieldName,
    IGeometricFieldDescriptor_out fieldDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::findField"
        "(const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Loop over all known field descriptors and find the one with the
        // specified name.
        for
        (
            ObjRefHashTable<IGeometricFieldDescriptorImpl*>::iterator iter =
                fieldDescriptorMap_.begin() ;
            iter != fieldDescriptorMap_.end() ;
            ++iter
        )
        {
            // Check for matching name.
            CORBA::String_var fieldDescName = iter()->name();
            if (strcmp(fieldDescName, fieldName) == 0)
            {
                fieldDescriptor = iter()->_this();
                break;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::addField
(
    const char* fieldName,
    IGeometricFieldDescriptor_out fieldDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::addField(const char* fieldName, "
        "IGeometricFieldDescriptor_out fieldDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // See if we already have a field with this name.
        if (fieldDescriptorMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new FieldDescriptor.
        IGeometricFieldDescriptorImpl* pFieldDescriptor =
            new IGeometricFieldDescriptorImpl(fieldName);

        if (pFieldDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create FieldDescriptor object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to field map.
        fieldDescriptorMap_.insert(fieldName, pFieldDescriptor);

        // Return a reference to the FieldDescriptor object.
        fieldDescriptor = pFieldDescriptor->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::deleteField(const char* fieldName)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::deleteField(const char* fieldName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have a field with this name.
        if (!fieldDescriptorMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from field map.
        HashTable<IGeometricFieldDescriptorImpl*>::iterator iter =
           fieldDescriptorMap_.find(fieldName);

        fieldDescriptorMap_.erase(iter);     // Releases object reference.
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationClassImpl::boundaryTypes()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::boundaryTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList
    (
        FoamXWordList((const wordList&)boundaryTypeDescriptorMap_.toc())
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::getBoundaryType
(
    const char* boundaryTypeName,
    IBoundaryTypeDescriptor_out boundaryDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::getBoundaryType"
        "(const char* boundaryTypeName, "
        "IBoundaryTypeDescriptor_out boundaryDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have this object cached.
        if (!boundaryTypeDescriptorMap_.found(boundaryTypeName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid boundary type name " + word(boundaryTypeName) + ".",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to object.
        boundaryDescriptor =
            boundaryTypeDescriptorMap_[boundaryTypeName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IApplicationClassImpl::findBoundaryType
(
    const char* boundaryTypeName,
    IBoundaryTypeDescriptor_out boundaryDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::findBoundaryType"
        "(const char* boundaryTypeName, "
        "IBoundaryTypeDescriptor_out boundaryDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Loop over all known boundary type descriptors and find the one
        // with the specified name.
        for
        (
            ObjRefHashTable<IBoundaryTypeDescriptorImpl*>::iterator iter =
                boundaryTypeDescriptorMap_.begin() ;
            iter != boundaryTypeDescriptorMap_.end() ;
            ++iter
        )
        {
            // Check for matching name.
            CORBA::String_var bndTypeName = iter()->name();
            if (strcmp(bndTypeName, boundaryTypeName) == 0)
            {
                boundaryDescriptor = iter()->_this();
                break;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::addBoundaryType
(
    const char* boundaryTypeName,
    IBoundaryTypeDescriptor_out boundaryDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::addBoundaryType"
        "(const char* boundaryTypeName, "
        "IBoundaryTypeDescriptor_out boundaryDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we already have a boundary type with this name.
        if (boundaryTypeDescriptorMap_.found(boundaryTypeName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid boundary type name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new BoundaryTypeDescriptor object.
        IBoundaryTypeDescriptorImpl* pBoundaryDescriptor =
            new IBoundaryTypeDescriptorImpl(boundaryTypeName);

        if (pBoundaryDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create BoundaryTypeDescriptor object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to map.
        boundaryTypeDescriptorMap_.insert
        (
            boundaryTypeName,
            pBoundaryDescriptor
        );

        // Return a reference to the BoundaryTypeDescriptor object.
        boundaryDescriptor = pBoundaryDescriptor->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::deleteBoundaryType
(
    const char* boundaryTypeName
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::deleteBoundaryType"
        "(const char* boundaryTypeName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "delete " << boundaryTypeName << endl;

        // See if we have a boundary type with this name.
        if (!boundaryTypeDescriptorMap_.found(boundaryTypeName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid boundary type name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from map.
        HashTable<IBoundaryTypeDescriptorImpl*>::iterator iter =
            boundaryTypeDescriptorMap_.find(boundaryTypeName);

        // Releases object reference.
        boundaryTypeDescriptorMap_.erase(iter);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationClassImpl::dictionaries()
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::dictionaries()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(dictionaryList_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::getDictionary
(
    const char* dictName,
    ITypeDescriptor_out dictTypeDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::getDictionary"
        "(const char* dictName, ITypeDescriptor_out dictTypeDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have this dictionary's type descriptor cached.
        if (!dictTypeDescriptorMap_.found(dictName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid dictionary name '" + word(dictName) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the TypeDescriptor object.
        dictTypeDescriptor = dictTypeDescriptorMap_[dictName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::addDictionary
(
    const char* dictName,
    ITypeDescriptor_out dictTypeDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::addDictionary"
        "(const char* dictName, ITypeDescriptor_out dictTypeDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we already have a dictionary with this name.
        if (dictTypeDescriptorMap_.found(dictName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid dictionary name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new TypeDescriptor object for a new
        // dictionary.
        ITypeDescriptorImpl* pTypeDescriptor = new ITypeDescriptorImpl
        (
            dictName,
            Type_Dictionary,
            dictName
        );

        if (pTypeDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create Dictionary TypeDescriptor object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to dictionary map.
        dictTypeDescriptorMap_.insert(dictName, pTypeDescriptor);

        // Add name to dictionary list.
        dictionaryList_.append(dictName);

        // Return a reference to the TypeDescriptor object.
        dictTypeDescriptor = pTypeDescriptor->_this();
    }
    CATCH_ALL(functionName);

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationClassImpl::deleteDictionary(const char* dictName)
{
    static const char* functionName =
        "FoamX::IApplicationClassImpl::deleteDictionary(const char* dictName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have a dictionary with this name.
        if (!dictTypeDescriptorMap_.found(dictName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid dictionary name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from dictionary map.
        HashTable<ITypeDescriptorImpl*>::iterator iter =
            dictTypeDescriptorMap_.find(dictName);

        // Releases object reference.
        dictTypeDescriptorMap_.erase(iter);

        // Remove name from dictionary list.
        dictionaryList_.remove(dictName);
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
