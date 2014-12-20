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

Class
    IPropertiesImpl

Description

SourceFiles
    IPropertiesImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IPropertiesImpl_H
#define IPropertiesImpl_H

// Foam header files.
#include "HashPtrTable.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "ObjRefHashTable.H"
#include "FoamXWordList.H"
#include "FoamXFileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Class forward declarations.
class IApplicationClassImpl;
class ITypeDescriptorImpl;
class IGeometryDescriptorImpl;
class IPatchDescriptorImpl;
class ITypeDescriptorImpl;

/*---------------------------------------------------------------------------*\
                           Class className Declaration
\*---------------------------------------------------------------------------*/

class IPropertiesImpl
:
    public POA_FoamXServer::CaseServer::IFoamProperties,
    public PortableServer::RefCountServantBase
{
    // Private data

        //- Allow or disallow editing.
        bool readOnly_;

        //- Available FoamX modules.
        FoamXWordList availableModules_;

        //- Root directories for cases (as read in)
        FoamXFileNameList rawRootDirectories_;

        //- Root directories for cases (after expansion)
        FoamXFileNameList rootDirectories_;

        //- Available foam types dictionary.
        Foam::dictionary foamTypesDict_;

        //- Available foam type names.
        FoamXWordList foamTypes_;

        //- Available geometry types.
        FoamXWordList geometryTypes_;

        //- Available patch types.
        FoamXWordList patchTypes_;

        //- Available patch field types.
        FoamXWordList patchFieldTypes_;

        //- Defined application class descriptors.
        Foam::HashPtrTable<FoamXServer::ApplicationClassDescriptor>
            appClassDescriptorMap_;

        //- Foam utility descriptor map.
        Foam::HashPtrTable<FoamXServer::UtilityDescriptor>
            utilityDescriptorMap_;


        //- Map of (created on demand) Application Class objects.
        ObjRefHashTable<IApplicationClassImpl*> appClassMap_;

        //- Map of Foam Type Descriptor objects.
        ObjRefHashTable<ITypeDescriptorImpl*> foamTypeMap_;

        //- Map of Geometry Type Descriptor objects.
        ObjRefHashTable<IGeometryDescriptorImpl*> geometryTypeMap_;

        //- Map of Patch Descriptor objects.
        ObjRefHashTable<IPatchDescriptorImpl*> patchMap_;

        //- Map of Patch Field Descriptor objects.
        ObjRefHashTable<ITypeDescriptorImpl*> patchFieldMap_;


    // Private member functions

        FoamXServer::ApplicationClassDescriptor* readAppClassDescriptor
        (
            const Foam::word& name,
            const Foam::fileName& category,
            const bool systemClass
        );

        void readAppClassDescriptors
        (
            const Foam::fileName& dir,
            const Foam::fileName& category,
            const bool systemClass
        );

        FoamXServer::UtilityDescriptor* readUtilityDescriptor
        (
            const Foam::word& name,
            const Foam::fileName& category,
            const bool systemClass
        );

        void readUtilityDescriptors
        (
            const Foam::fileName& dir,
            const Foam::fileName& category,
            const bool systemClass
        );

        IPropertiesImpl(const IPropertiesImpl&);


public:

    // Constructors

        IPropertiesImpl(bool readOnly = true);


    // Destructor

        virtual ~IPropertiesImpl();


    // Member Functions

        const Foam::dictionary& foamTypesDict() const
        {
            return foamTypesDict_;
        }

        const ObjRefHashTable<ITypeDescriptorImpl*>& foamTypes() const
        {
            return foamTypeMap_;
        }

        const ObjRefHashTable<IGeometryDescriptorImpl*>&
        geometryDescriptors() const
        {
            return geometryTypeMap_;
        }


        //- Installed FoamX modules.
        virtual FoamXServer::StringList* availableModules();


        // Case root directories.

            virtual FoamXServer::StringList* rootDirectories();
            virtual FoamXServer::StringList* rawRootDirectories();

            virtual void addRootDirectory(const char* rootDir);

            virtual void deleteRootDirectory(const char* rootDir);


        // Foam types.

            virtual FoamXServer::StringList* foamTypes();

            virtual void getFoamType
            (
                const char* foamTypeKey,
                FoamXServer::ITypeDescriptor_out typeDesc
            );


        // Geometry types.

            virtual FoamXServer::StringList* geometryTypes();

            virtual void getGeometryType
            (
                const char* geometryTypeKey,
                FoamXServer::CaseServer::IGeometryDescriptor_out
                    geometryDesc
            );


        // Patch types.

            virtual FoamXServer::StringList* patchTypes();

            virtual void getPatchType
            (
                const char* patchTypeKey,
                FoamXServer::CaseServer::IPatchDescriptor_out
                    patchDesc
            );

            virtual void findPatchType
            (
                const char* patchTypeName,
                FoamXServer::CaseServer::IPatchDescriptor_out
                    patchDesc
            );


        // Patch field types.

            virtual FoamXServer::StringList* patchFieldTypes();

            virtual void getPatchFieldType
            (
                const char* patchFieldTypeKey,
                FoamXServer::ITypeDescriptor_out patchFieldDesc
            );

            virtual void findPatchFieldType
            (
                const char* patchFieldTypeName,
                FoamXServer::ITypeDescriptor_out patchFieldDesc
            );


        // Application classes.

            virtual FoamXServer::ApplicationClassDescriptorList* 
                applicationClasses();

            virtual void getApplicationClass
            (
                const char* appClassKey,
                FoamXServer::CaseServer::IApplicationClass_out appClass
            );

            virtual void addApplicationClass
            (
                const char* appClassKey,
                FoamXServer::CaseServer::IApplicationClass_out appClass
            );

            virtual void deleteApplicationClass(const char* appClassKey);

            virtual void cloneApplicationClass
            (
                const char* appClassKeySrc,
                const char* appClassKeyDest,
                FoamXServer::CaseServer::IApplicationClass_out appClass
            );


        // Foam utilities.

            virtual FoamXServer::UtilityDescriptorList* utilities();

            virtual void getUtilityControlDict
            (
                const char* utilityName,
                const char* rootDir,
                const char* caseName,
                FoamXServer::IDictionaryEntry_out controlDict
            );


        //- Validation.
        virtual void validate();


        // Persistence.

            virtual void saveSystemProperties();

            virtual void saveUserProperties();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
