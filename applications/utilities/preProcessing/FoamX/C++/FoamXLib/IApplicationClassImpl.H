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
    ITypeDescriptorImpl

Description

SourceFiles
    IApplicationClassImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IApplicationClassImpl_H
#define IApplicationClassImpl_H

// Foam header files.
#include "fileName.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXString.H"
#include "FoamXStringList.H"
#include "FoamXWordList.H"
#include "ObjRefHashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Class forward declarations.
class IGeometricFieldDescriptorImpl;
class IBoundaryTypeDescriptorImpl;
class ITypeDescriptorImpl;
class IPropertiesImpl;

/*---------------------------------------------------------------------------*\
                 Class IApplicationClassImpl Declaration
\*---------------------------------------------------------------------------*/

class IApplicationClassImpl
:
    public POA_FoamXServer::CaseServer::IApplicationClass,
    public PortableServer::RefCountServantBase
{
    // Private data

        Foam::word name_;
        Foam::string description_;
        Foam::fileName category_;

        //- System or user class flag.
        bool systemClass_;

        //- Reference to the Foam properties class
        const IPropertiesImpl& foamSystemProperties_;

        ITypeDescriptorImpl* arguments_;
        FoamXWordList moduleList_;
        FoamXWordList dictionaryList_;

        //- Full application class definition dictionary name.
        Foam::fileName classDictName_;

        ObjRefHashTable<IGeometricFieldDescriptorImpl*> fieldDescriptorMap_;

        ObjRefHashTable<IBoundaryTypeDescriptorImpl*>
            boundaryTypeDescriptorMap_;

        ObjRefHashTable<ITypeDescriptorImpl*> dictTypeDescriptorMap_;


public:

    // Constructors

        IApplicationClassImpl
        (
            const ApplicationClassDescriptor& appClassDesc,
            const IPropertiesImpl& foamSystemProperties
        );


    // Destructor

        virtual ~IApplicationClassImpl();


    // Member Functions

            virtual char* name();
            virtual void name(const char*);

            virtual char* description();
            virtual void description(const char*);

            virtual char* category();
            virtual void category(const char*);

            virtual FoamXServer::StringList* modules();
            virtual void modules(const FoamXServer::StringList&);

            virtual CORBA::Boolean systemClass();

            virtual FoamXServer::StringList* fields();

            virtual void getField
            (
                const char* fieldKey,
                FoamXServer::CaseServer::IGeometricFieldDescriptor_out
                    fieldDescriptor
            );

            virtual void findField
            (
                const char* fieldName,
                FoamXServer::CaseServer::IGeometricFieldDescriptor_out
                    fieldDescriptor
            );

            virtual void addField
            (
                const char* fieldKey,
                FoamXServer::CaseServer::IGeometricFieldDescriptor_out
                    fieldDescriptor
            );

            virtual void deleteField(const char* fieldKey);

            virtual FoamXServer::StringList* boundaryTypes();

            virtual void getBoundaryType
            (
                const char* boundaryTypeKey,
                FoamXServer::CaseServer::IBoundaryTypeDescriptor_out
                    boundaryDescriptor
            );

            virtual void findBoundaryType
            (
                const char* boundaryTypeName,
                FoamXServer::CaseServer::IBoundaryTypeDescriptor_out
                    boundaryDescriptor
            );

            virtual void addBoundaryType
            (
                const char* boundaryTypeKey,
                FoamXServer::CaseServer::IBoundaryTypeDescriptor_out
                    boundaryDescriptor
            );

            virtual void deleteBoundaryType(const char* boundaryTypeKey);

            virtual FoamXServer::StringList* dictionaries();

            virtual void getDictionary
            (
                const char* dictName,
                FoamXServer::ITypeDescriptor_out
                    dictTypeDescriptor
            );
    
            virtual void addDictionary
            (
                const char* dictName,
                FoamXServer::ITypeDescriptor_out
                    dictTypeDescriptor
            );

            virtual void deleteDictionary(const char* dictName);


        //- Validation.
        virtual void validate();

        //- Persistence.
        virtual void save();
 
        //- Load from application class definition dictionary.
        void load();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
