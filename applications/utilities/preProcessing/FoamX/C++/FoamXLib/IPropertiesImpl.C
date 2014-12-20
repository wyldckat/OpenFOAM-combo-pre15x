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
#include "IFstream.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "IPropertiesImpl.H"
#include "ITypeDescriptorImpl.H"
#include "IApplicationClassImpl.H"
#include "IDictionaryEntryImpl.H"
#include "RootDictionary.H"
#include "IGeometryDescriptorImpl.H"
#include "IPatchDescriptorImpl.H"
#include "ITypeDescriptorImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::ApplicationClassDescriptor*
FoamX::IPropertiesImpl::readAppClassDescriptor
(
    const Foam::word& name,
    const Foam::fileName& category,
    const bool systemClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::readAppClassDescriptor"
        "(const Foam::word& name, const Foam::fileName& category, "
        "const bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Create and initialise a new ApplicationClassDescriptor object.
        ApplicationClassDescriptor* pDesc = new ApplicationClassDescriptor();

        pDesc->name = name.c_str();
        pDesc->category = category.c_str();
        pDesc->systemClass = systemClass;

        return pDesc;
    }
    CATCH_ALL(functionName);
}


void FoamX::IPropertiesImpl::readAppClassDescriptors
(
    const Foam::fileName& dir,
    const Foam::fileName& category,
    const bool systemClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::readAppClassDescriptor"
        "(const Foam::fileName& dir, "
        "const Foam::fileName& category, "
        "const bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileNameList appClassDirs = readDir(dir, fileName::DIRECTORY);

        forAll (appClassDirs, i)
        {
            if
            (
                Foam::file(dir/appClassDirs[i]/(appClassDirs[i] + ".cfg"))
            && !appClassDescriptorMap_.found(appClassDirs[i])
            )
            {
                appClassDescriptorMap_.insert
                (
                    appClassDirs[i],
                    readAppClassDescriptor
                    (
                        appClassDirs[i],
                        category,
                        systemClass
                    )
                );
            }

            readAppClassDescriptors
            (
                dir/appClassDirs[i],
                category/appClassDirs[i],
                systemClass
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::UtilityDescriptor* FoamX::IPropertiesImpl::readUtilityDescriptor
(
    const Foam::word& name,
    const Foam::fileName& category,
    const bool systemClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::readUtilityDescriptor"
        "(const Foam::word& name, const Foam::fileName& category, "
        "const bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName utilityPath;

        if (systemClass)
        {
            utilityPath = Paths::system;
        }
        else
        {
            utilityPath = Paths::user;
        }

        utilityPath += "/utilities"/category;

        fileName utilityCfgPath = utilityPath/(name + ".cfg");

        // Create and initialise a new UtilityDescriptor object.
        UtilityDescriptor* pDesc = new UtilityDescriptor();

        pDesc->name = name.c_str();
        pDesc->category = category.c_str();
        pDesc->systemClass = systemClass;

        dictionary utilityConfigDict((IFstream(utilityCfgPath)()));

//        pDesc->changesMesh =
//            readBool(utilityConfigDict.lookup("changesMesh"));
        pDesc->changesMesh = false;

//        pDesc->changesFields =
//            readBool(utilityConfigDict.lookup("changesFields"));
        pDesc->changesFields = false;

        pDesc->description =
            Foam::string(utilityConfigDict.lookup("description")).c_str();

        pDesc->clientBean = 
            Foam::string(utilityConfigDict.lookup("clientBean")).c_str();

        word controlDictName = name + "Dict";

        if (utilityConfigDict.found(controlDictName))
        {
            ITypeDescriptorImpl* itdPtr = new ITypeDescriptorImpl
            (
                controlDictName,
                utilityCfgPath,
                utilityConfigDict.subDict(controlDictName),
                foamTypesDict_
            );

            pDesc->controlDict = itdPtr->_this();
        }
        else
        {
            pDesc->controlDict = ITypeDescriptor::_nil();
        }

        return pDesc;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::readUtilityDescriptors
(
    const Foam::fileName& dir,
    const Foam::fileName& category,
    const bool systemClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::readUtilityDescriptor"
        "(const Foam::fileName& dir, "
        "const Foam::fileName& category, "
        "const bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileNameList utilityCfgs = readDir(dir, fileName::FILE);

        // Initialise the utility descriptor objects.
        forAll (utilityCfgs, i)
        {
            word utilName = utilityCfgs[i].lessExt();

            if (!utilityDescriptorMap_.found(utilName))
            {
                utilityDescriptorMap_.insert
                (
                    utilName,
                    readUtilityDescriptor(utilName, category, systemClass)
                );
            }
        }

        fileNameList utilityDirs = readDir(dir, fileName::DIRECTORY);

        forAll (utilityDirs, i)
        {
            readUtilityDescriptors
            (
                dir/utilityDirs[i],
                category/utilityDirs[i],
                systemClass
            );
        }
    }
    CATCH_ALL(functionName);
}


FoamX::IPropertiesImpl::IPropertiesImpl(bool readOnly)
:
    readOnly_(readOnly)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::IPropertiesImpl(bool readOnly)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Open the FoamX system config dictionary.
        fileName fxConfigFileName = Paths::system/"FoamX.cfg";
        if (!exists(fxConfigFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find FoamX root configuration dictionary '"
              + fxConfigFileName + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }
        dictionary configDict((IFstream(fxConfigFileName)()));

        // Read list data from config dictionary.
        availableModules_.read(configDict.lookup("availableModules"));


        // Open the Foam types dictionary.
        Foam::fileName foamTypesDictFileName =
            Paths::system/"types/types.cfg";

        if (!exists(foamTypesDictFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find Foam types dictionary '"
              + foamTypesDictFileName + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }
        (IFstream(foamTypesDictFileName)()) >> foamTypesDict_;

        foamTypes_ = (const Foam::wordList&)foamTypesDict_.toc();

        // Initialise the foam type descriptor objects.
        forAll(foamTypes_, i)
        {
            word foamTypeName(foamTypes_[i]);

            // Check for a definition dictionary.
            if (!foamTypesDict_.found(foamTypeName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Foam type descriptor dictionary not found for field type '"
                  + foamTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new TypeDescriptor object.
            ITypeDescriptorImpl* pFieldTypeDescriptor = new ITypeDescriptorImpl
            (
                foamTypeName,
                "FoamX:foamTypes",
                foamTypesDict_.subDict(foamTypeName),
                foamTypesDict_
            );

            if (pFieldTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create Foam TypeDescriptor object for "
                    "field type '" + foamTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            foamTypeMap_.insert(foamTypeName, pFieldTypeDescriptor);
        }


        dictionary geometryDict
        (
            IFstream(Paths::system/"types/geometries.cfg")()
        );
        geometryTypes_ = (const Foam::wordList&)geometryDict.toc();

        // Initialise the geometry descriptor objects.
        forAll(geometryTypes_, i)
        {
            word geometricTypeName(geometryTypes_[i]);

            // Check for a definition dictionary.
            if (!geometryDict.found(geometricTypeName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Geometry descriptor dictionary not found for "
                    "geometry type '" + geometricTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new GeometryDescriptor object.
            IGeometryDescriptorImpl* pGeometryTypeDescriptor = new
            IGeometryDescriptorImpl
            (
                geometricTypeName,
                geometryDict
            );

            if (pGeometryTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create field GeometryDescriptor object for "
                    "geometry type '" + geometricTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            geometryTypeMap_.insert(geometricTypeName, pGeometryTypeDescriptor);
        }

        dictionary patchDict
        (
            IFstream(Paths::system/"types/patches.cfg")()
        );
        patchTypes_ = (const Foam::wordList&)patchDict.toc();

        // Initialise the patch descriptor objects.
        forAll(patchTypes_, i)
        {
            word patchTypeName(patchTypes_[i]);

            // Check for a definition dictionary.
            if (!patchDict.found(patchTypeName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Patch descriptor dictionary not found for patch type '"
                   + patchTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new PatchDescriptor object.
            IPatchDescriptorImpl* pPatchDescriptor =
                new IPatchDescriptorImpl(patchTypeName, patchDict);

            if (pPatchDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create PatchDescriptor object for patch type '"
                   + patchTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            patchMap_.insert(patchTypeName, pPatchDescriptor);
        }


        dictionary patchFieldDict
        (
            IFstream(Paths::system/"types/patchFields.cfg")()
        );
        patchFieldTypes_ = (const Foam::wordList&)patchFieldDict.toc();

        // Initialise the patch field descriptor objects.
        forAll(patchFieldTypes_, i)
        {
            word patchFieldTypeName(patchFieldTypes_[i]);

            // Check for a definition dictionary.
            if (!patchFieldDict.found(patchFieldTypeName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Patch field descriptor dictionary not found for "
                    "patch field type '" + patchFieldTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new PatchFieldDescriptor object.
            ITypeDescriptorImpl* pPatchFieldDescriptor = 
            new ITypeDescriptorImpl
            (
                patchFieldTypeName, 
                patchFieldDict.name(),
                patchFieldDict.subDict(patchFieldTypeName),
                foamTypesDict_
            );

            if (pPatchFieldDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create PatchFieldDescriptor object for "
                    "patch field type '" + patchFieldTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            patchFieldMap_.insert(patchFieldTypeName, pPatchFieldDescriptor);
        }


        // ---------------------------------------------------------------------
        // Open the user's Foam control dictionary.
        fileName controlDictFileName = Foam::dotFoam("controlDict");
        if (exists(controlDictFileName))
        {
            dictionary controlDict((IFstream(controlDictFileName)()));

            stringList rawRootDirs(1, ".");

            if (controlDict.found("caseRoots"))
            {
                // Read root directory list.
                // Filter out directories which do not exist.
                stringList controlDictRootDirs(controlDict.lookup("caseRoots"));

                if (controlDictRootDirs.size())
                {
                    rawRootDirs = controlDictRootDirs;
                }
            }

            // Check that each root directory actually exists.
            forAll (rawRootDirs, i)
            {
                fileName rootDir = rawRootDirs[i];
                rootDir.expand();

                if (exists(rootDir))
                {
                    rootDirectories_.length(rootDirectories_.length() + 1);
                    rootDirectories_[rootDirectories_.length() - 1]
                        = rootDir.c_str();

                    rawRootDirectories_.length
                    (
                        rawRootDirectories_.length() + 1
                    );

                    rawRootDirectories_[rawRootDirectories_.length() - 1]
                        = rawRootDirs[i].c_str();
                }
                else
                {
                    Warning
                        << "User specified root directory " << rawRootDirs[i];

                    if (rootDir.size() != rawRootDirs[i].size())
                    {
                        Info<< " (" << rootDir << ")";
                    }

                    Info<< " not found. Removing from root directories list."
                        << endl;
                }
            }
        }


        // ---------------------------------------------------------------------
        // Initialise the application class descriptor map.

        // Scan the user's utilities directory for valid utilities
        readAppClassDescriptors(Paths::user/"applications", "", false);

        // Scan the system utilities directory for valid utilities
        readAppClassDescriptors(Paths::system/"applications", "", true);


        // ---------------------------------------------------------------------
        // Initialise the utility descriptor map.

        // Scan the user's utilities directory for valid utilities
        readUtilityDescriptors(Paths::user/"utilities", "", false);

        // Scan the system utilities directory for valid utilities
        readUtilityDescriptors(Paths::system/"utilities", "", true);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPropertiesImpl::~IPropertiesImpl()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::~IPropertiesImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::availableModules()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::availableModules()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(availableModules_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::rootDirectories()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::rootDirectories()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(rootDirectories_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::rawRootDirectories()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::rawRootDirectories()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(rawRootDirectories_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::addRootDirectory(const char* rawRootDir)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::addRootDirectory(const char* rawRootDir)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName rootDir = rawRootDir;
    rootDir.expand();    

    if (!rootDirectories_.found(rootDir.c_str()))
    {
        rootDirectories_.append(rootDir.c_str());
        rawRootDirectories_.append(rawRootDir);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::deleteRootDirectory(const char* rawRootDir)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::deleteRootDirectory(const char* rawRootDir)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName rootDir = rawRootDir;
    rootDir.expand();    

    if (rootDirectories_.found(rootDir.c_str()))
    {
        rootDirectories_.remove(rootDir.c_str());
        rawRootDirectories_.remove(rawRootDir);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::foamTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::foamTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(foamTypes_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getFoamType
(
    const char* foamTypeKey,
    ITypeDescriptor_out typeDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getFoamType"
        "(const char* foamTypeKey, ITypeDescriptor_out typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Check for valid foam type name.
        if (!foamTypeMap_.found(foamTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid foam type name '" + word(foamTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the TypeDescriptor object.
        typeDesc = foamTypeMap_[foamTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::geometryTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::geometryTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(geometryTypes_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getGeometryType
(
    const char* geometryTypeKey,
    IGeometryDescriptor_out geometryDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getGeometryType"
        "(const char* geometryTypeKey, IGeometryDescriptor_out geometryDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Check for valid geometry type name.
        if (!geometryTypeMap_.found(geometryTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid geometry type name '" + word(geometryTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the GeometryDescriptor object.
        geometryDesc = geometryTypeMap_[geometryTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::patchTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::patchTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(patchTypes_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IPropertiesImpl::getPatchType
(
    const char* patchTypeKey,
    IPatchDescriptor_out patchDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getPatchType"
        "(const char* patchTypeKey, IPatchDescriptor_out patchDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Check for valid patch name.
        if (!patchMap_.found(patchTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch type name '" + word(patchTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the PatchDescriptor object.
        patchDesc = patchMap_[patchTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IPropertiesImpl::findPatchType
(
    const char* patchTypeName,
    IPatchDescriptor_out patchDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::findPatchType"
        "(const char* patchTypeName, IPatchDescriptor_out patchDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        IPatchDescriptorImpl* pPatchDescriptor = NULL;

        // Loop over all known patch types and find the one with the
        // specified name.
        for (unsigned int i = 0; i <patchTypes_.length(); i++)
        {
            IPatchDescriptor_var patchDesc;
            getPatchType(patchTypes_[i], patchDesc.out());

            // Check for matching name.
            CORBA::String_var typeName = patchDesc->name();
            if (strcmp(typeName, patchTypeName) == 0)
            {
                pPatchDescriptor = patchMap_[word(patchTypes_[i])];
                break;
            }
        }

        // Return a reference to the PatchFieldDescriptor object.
        if (pPatchDescriptor != NULL)
        {
            patchDesc = pPatchDescriptor->_this();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::patchFieldTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::patchFieldTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(patchFieldTypes_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getPatchFieldType
(
    const char* patchFieldTypeKey,
    ITypeDescriptor_out patchFieldDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getPatchFieldType"
        "(const char* patchFieldTypeKey,"
        "ITypeDescriptor_out patchFieldDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check for valid patch field name.
        if (!patchFieldMap_.found(patchFieldTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch field type name '"
              + word(patchFieldTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the PatchFieldDescriptor object.
        patchFieldDesc = patchFieldMap_[patchFieldTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::findPatchFieldType
(
    const char* patchFieldTypeName,
    ITypeDescriptor_out patchFieldDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::findPatchFieldType"
        "(const char* patchFieldTypeName, "
        "ITypeDescriptor_out patchFieldDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        ITypeDescriptorImpl* pPatchFieldDescriptor = NULL;

        // Loop over all knwon patch field types and find the one with
        // the specified name.
        for (unsigned int i = 0; i <patchFieldTypes_.length(); i++)
        {
            ITypeDescriptor_var patchFieldDesc;
            getPatchFieldType(patchFieldTypes_[i], patchFieldDesc.out());

            // Check for matching name.
            CORBA::String_var typeName = patchFieldDesc->name();
            if (strcmp(typeName, patchFieldTypeName) == 0)
            {
                pPatchFieldDescriptor =
                    patchFieldMap_[word(patchFieldTypes_[i])];
                break;
            }
        }

        // Return a reference to the PatchFieldDescriptor object.
        if (pPatchFieldDescriptor != NULL)
        {
            patchFieldDesc = pPatchFieldDescriptor->_this();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::ApplicationClassDescriptorList*
FoamX::IPropertiesImpl::applicationClasses()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::applicationClasses()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Construct an application class list and return.
    ApplicationClassDescriptorList* pAppClassList =
        new ApplicationClassDescriptorList();
    pAppClassList->length(appClassDescriptorMap_.size());

    label i = 0;

    for
    (
        Foam::HashPtrTable<ApplicationClassDescriptor>::iterator iter = 
            appClassDescriptorMap_.begin();
        iter != appClassDescriptorMap_.end();
        ++iter
    )
    {
        (*pAppClassList)[i++] = *iter();
    }

    return pAppClassList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getApplicationClass
(
    const char* appClassKey,
    IApplicationClass_out appClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getApplicationClass"
        "(const char* appClassKey, IApplicationClass_out appClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        IApplicationClassImpl* pAppClass = NULL;

        // See if the specified application class name is valid.
        if (!appClassDescriptorMap_.found(appClassKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "getApplicationClass::Invalid application class name '"
              + word(appClassKey) + "'.",
                "IPropertiesImpl",
                __FILE__, __LINE__
            );
        }

        // See if we have this application class object cached.
        if (appClassMap_.found(appClassKey))
        {
            log << "Existing." << endl;
            pAppClass = appClassMap_[appClassKey];
        }
        else
        {
            log << "New application class " << appClass << endl;

            // Create and initialise a new ApplicationClass object.
            pAppClass = new IApplicationClassImpl
            (
                *appClassDescriptorMap_[appClassKey],
                *this
            );

            if (pAppClass == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create ApplicationClass object for '"
                  + word(appClassKey) + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Load the application class from the definition dictionary.
            // Allow user defined application classes to override the system
            // classes.
            pAppClass->load();

            // Add application class object to map.
            appClassMap_.insert(appClassKey, pAppClass);
        }

        // Return a reference to the ApplicationClass object.
        appClass = pAppClass->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::addApplicationClass
(
    const char* appClassKey,
    IApplicationClass_out appClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::addApplicationClass"
        "(const char* appClassKey, IApplicationClass_out appClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to addApplicationClass for '" + word(appClassKey)
              + "'. Object is read only.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified application class name is valid.
        if (appClassDescriptorMap_.found(appClassKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid application class name '" + word(appClassKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new ApplicationClassDescriptor object.
        ApplicationClassDescriptor* pDesc = new ApplicationClassDescriptor();
        pDesc->name        = appClassKey;
        pDesc->category    = "";
        pDesc->systemClass = false;

        // Create and initialise a new ApplicationClass object.
        // User defined application class.
        IApplicationClassImpl* pAppClass = new IApplicationClassImpl
        (
            *pDesc,
            *this
        );

        if (pAppClass == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create ApplicationClass object for '"
              + word(appClassKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add application class object to map.
        appClassMap_.insert(appClassKey, pAppClass);

        // Add application class descriptor to map.
        appClassDescriptorMap_.insert(appClassKey, pDesc);

        // Return a reference to the application class object.
        appClass = pAppClass->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::deleteApplicationClass(const char* appClassKey)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::deleteApplicationClass"
        "(const char* appClassKey)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a wobbler.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to deleteApplicationClass for " 
              + word(appClassKey) + "'. Object is read only.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified application class name is valid.
        if (!appClassDescriptorMap_.found(appClassKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid application class name '" + word(appClassKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified application class name is user defined.
        if (appClassDescriptorMap_[appClassKey]->systemClass)
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Unable to delete system application class '"
              + word(appClassKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if we have this application class object cached.
        if (appClassMap_.found(appClassKey))
        {
            ObjRefHashTable<IApplicationClassImpl*>::iterator iter
            (
                appClassMap_.find(appClassKey)
            );

            appClassMap_.erase(iter);   // Releases implementation reference.
        }

        // Remove application class descriptor.
        Foam::HashPtrTable<ApplicationClassDescriptor>::iterator iter
        (
            appClassDescriptorMap_.find(appClassKey)
        );
        appClassDescriptorMap_.erase(iter);

        // Remove the complete aplicationClass directory structure
        if (!rmDir(Paths::user/"applications"/appClassKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Unable to delete application class "
              + Paths::user/"applications"/appClassKey,
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::cloneApplicationClass
(
    const char* appClassKeySrc,
    const char* appClassKeyDest,
    IApplicationClass_out appClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::cloneApplicationClass"
        "(const char* appClassKeySrc, const char* appClassKeyDest, "
        "IApplicationClass_out appClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to cloneApplicationClass for '"
              + word(appClassKeySrc) + "'. Object is read only.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified source application class name is valid.
        if (!appClassDescriptorMap_.found(appClassKeySrc))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid source application class name '"
              + word(appClassKeySrc) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified destination application class name is valid.
        if (appClassDescriptorMap_.found(appClassKeyDest))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid destination application class name '"
              + word(appClassKeyDest) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Copy application class definition files.
        ApplicationClassDescriptor* pDescriptor =
            appClassDescriptorMap_[appClassKeySrc];

        fileName sourcePath =
            (pDescriptor->systemClass ? Paths::system : Paths::user)
           /"applications"/appClassKeySrc;

        fileName destPath = Paths::user/"applications"/appClassKeyDest;

        // Make application class directory.
        if (!dir(destPath) && !mkDir(destPath))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to create directory '" + destPath + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Copy main configuration file.
        fileName srcConfig(fileName(appClassKeySrc) + ".cfg");
        fileName destConfig(fileName(appClassKeyDest) + ".cfg");
        cp(sourcePath/srcConfig, destPath/destConfig); // File to file copy.

        // Copy rest of files.
        fileNameList contents = readDir(sourcePath, fileName::FILE);
        forAll(contents, i)
        {
             // Do not copy original app class config file.
            if (contents[i] != srcConfig)
            {
                cp(sourcePath/contents[i], destPath/contents[i]);
            }
        }

        // Create and initialise a new ApplicationClassDescriptor object.
        // User defined application class.
        ApplicationClassDescriptor* pDesc = new ApplicationClassDescriptor();
        pDesc->name        = appClassKeyDest;
        pDesc->category    = "";
        pDesc->systemClass = false;

        // Create and initialise a new ApplicationClass object.
        IApplicationClassImpl* pAppClass = new IApplicationClassImpl
        (
            *pDesc,
            *this
        );

        if (pAppClass == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create ApplicationClass object '"
              + word(appClassKeyDest) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Load the application class details from the definition dictionary.
        // Allow user defined application classes to override the system
        // classes.
        pAppClass->load();

        // Add to application class map.
        appClassMap_.insert(appClassKeyDest, pAppClass);

        // Add application class descriptor to map.
        appClassDescriptorMap_.insert(appClassKeyDest, pDesc);

        // Return a reference to the application class object.
        appClass = pAppClass->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::UtilityDescriptorList* FoamX::IPropertiesImpl::utilities()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::utilities()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Construct an application class list and return.
    UtilityDescriptorList* pUtilityList = new UtilityDescriptorList();
    pUtilityList->length(utilityDescriptorMap_.size());

    label i = 0;

    for
    (
        Foam::HashPtrTable<UtilityDescriptor>::iterator iter = 
            utilityDescriptorMap_.begin();
        iter != utilityDescriptorMap_.end();
        ++iter
    )
    {
        (*pUtilityList)[i++] = *iter();
    }

    return pUtilityList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getUtilityControlDict
(
    const char* utilityName,
    const char* rootDir,
    const char* caseName,
    FoamXServer::IDictionaryEntry_out controlDict
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getUtilityDescriptor"
        "(const char* utilityName, FoamXServer::UtilityDescriptor_out utilityDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check for valid utility name.
        if (!utilityDescriptorMap_.found(utilityName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid foam utility name '" + word(utilityName) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        UtilityDescriptor* utilityDescriptor = utilityDescriptorMap_[utilityName];

        if (!CORBA::is_nil(utilityDescriptor->controlDict))
        {
            // Create an appropriate sub entry object and store its reference.
            IDictionaryEntryImpl* controlDictPtr = new RootDictionary
            (
                utilityDescriptor->controlDict,
                rootDir,
                caseName
            );

            if (controlDictPtr == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Couldn't create IDictionaryEntryImpl object for utility "
                  + word(utilityName),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            controlDict = controlDictPtr->_this();

            // Get dictionary full file name
            // (using the dictionary path information in the type descriptor).
            fileName dictName =
                fileName(rootDir)/fileName(caseName)
               /fileName(controlDict->typeDescriptor()->dictionaryPath())
               /controlDict->typeDescriptor()->name();

            // See if the dictionary file exists.
            if (exists(dictName))
            {
                // Load the values from the file.
                controlDictPtr->load((IFstream(dictName)()));
            }
        }
        else
        {
            controlDict = FoamXServer::IDictionaryEntry::_nil();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::validate()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::saveSystemProperties()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::saveSystemProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to save. IPropertiesImpl object is read only",
                functionName,
                __FILE__, __LINE__
            );
        }

        {
        fileName configDictFileName = Paths::system/"FoamX.cfg";
        DictionaryWriter dict(configDictFileName);

        dict.writeHeader
        (
            "FoamX System Properties.",
            "dictionary"
        );

        dict.writeEntry("availableModules", availableModules_);

        dict.startSubDict("processControl");
        dict.writeEntry("remoteShell", Foam::string("rsh"));
        dict.endSubDict();
        dict.writeBar();

        // Geometry type definitions.
        dict.writeSectionHeader("Geometry type definitions.");
        dict.startSubDict("geometryTypes");
        for
        (
            ObjRefHashTable<IGeometryDescriptorImpl*>::iterator iter =
                geometryTypeMap_.begin();
            iter != geometryTypeMap_.end();
            ++iter
        )
        {
            iter()->save(dict);
            dict.writeEndl();
        }
        dict.endSubDict();
        dict.writeBar();

        // Patch type definitions.
        dict.writeSectionHeader("Patch type definitions.");
        dict.startSubDict("patchTypes");
        for
        (
            ObjRefHashTable<IPatchDescriptorImpl*>::iterator iter =
                patchMap_.begin();
            iter != patchMap_.end();
            ++iter
        )
        {
            iter()->save(dict);
            dict.writeEndl();
        }
        dict.endSubDict();
        dict.writeBar();

        // Patch field type definitions.
        dict.writeSectionHeader("Patch field type definitions.");
        dict.startSubDict("patchFieldTypes");
        for
        (
            ObjRefHashTable<ITypeDescriptorImpl*>::iterator iter =
                patchFieldMap_.begin();
            iter != patchFieldMap_.end();
            ++iter
        )
        {
            iter()->save(dict, true);
            dict.writeEndl();
        }
        dict.endSubDict();
        dict.writeEndl();
        dict.writeEndBar();
        }



        {
        fileName configDictFileName = Paths::system/"types/types.cfg";
        DictionaryWriter dict(configDictFileName);

        dict.writeHeader
        (
            "Primitive types.",
            "dictionary"
        );

        for
        (
            ObjRefHashTable<ITypeDescriptorImpl*>::iterator iter =
                foamTypeMap_.begin();
            iter != foamTypeMap_.end();
            ++iter
        )
        {
            iter()->save(dict, false);
            dict.writeEndl();
        }

        dict.writeEndl();
        dict.writeEndBar();
        }


        // Foam utilities.
        for
        (
            Foam::HashPtrTable<UtilityDescriptor>::iterator iter = 
                utilityDescriptorMap_.begin();
            iter != utilityDescriptorMap_.end();
            ++iter
        )
        {
            fileName utilityDictFileName = 
                Paths::system/"utilities"/(word(iter()->name) + ".cfg");
            DictionaryWriter utilityDict(utilityDictFileName);

            utilityDict.writeHeader
            (
                "FoamX Utility Configuration File " + word(iter()->name),
                "dictionary"
            );

            utilityDict.writeEntry("changesMesh", bool(iter()->changesMesh));
            utilityDict.writeEntry("changesFields", bool(iter()->changesFields));
            utilityDict.writeEntry("description", iter()->description);
            utilityDict.writeEntry("category", iter()->category);

            if (!CORBA::is_nil(iter()->controlDict))
            {
                utilityDict.writeEntry
                (
                    iter()->controlDict->name(),
                    iter()->controlDict
                );
            }

            utilityDict.writeEndl();
            utilityDict.writeEndBar();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::saveUserProperties()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::saveUserProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to save. IPropertiesImpl object is read only",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName controlDictFileName = Foam::dotFoam("controlDict");
        DictionaryWriter dict(controlDictFileName);

        dict.writeHeader
        (
            "FoamX User Properties.",
            "dictionary"
        );

        dict.writeEndl();
        dict.writeEntry("caseRoots", rootDirectories_);
        dict.writeEndl();
        dict.writeEndBar();
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
