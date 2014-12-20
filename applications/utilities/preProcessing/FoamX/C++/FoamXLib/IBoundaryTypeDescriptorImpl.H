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
    IBoundaryTypeDescriptorImpl

Description

SourceFiles
    IBoundaryTypeDescriptorImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IBoundaryTypeDescriptorImpl_H
#define IBoundaryTypeDescriptorImpl_H

// Foam header files.
#include "string.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Class forward declarations.
class ITypeDescriptorImpl;

/*---------------------------------------------------------------------------*\
              Class IBoundaryTypeDescriptorImpl Declaration
\*---------------------------------------------------------------------------*/

class IBoundaryTypeDescriptorImpl
:
    public POA_FoamXServer::CaseServer::IBoundaryTypeDescriptor,
    public PortableServer::RefCountServantBase
{
    // Private data

        Foam::word name_;
        Foam::string displayName_;
        Foam::string description_;
        Foam::word patchType_;
        Foam::word superType_;
        FoamXServer::StringPairList patchFieldTypes_;


public:

    // Constructors

        IBoundaryTypeDescriptorImpl
        (
            const word& boundaryTypeName
        );

        IBoundaryTypeDescriptorImpl
        (
            const word& boundaryTypeName,
            const Foam::wordList& fieldList
        );

        IBoundaryTypeDescriptorImpl
        (
            const word& boundaryTypeName,
            const Foam::dictionary& boundaryTypeDict,
            const Foam::wordList& fieldList
        );


    // Destructor

        virtual ~IBoundaryTypeDescriptorImpl();


    // Member Functions

            virtual char* name();
            virtual void  name(const char*);
            virtual char* displayName();
            virtual void  displayName(const char*);
            virtual char* description();
            virtual void description(const char*);
            virtual char* patchType();
            virtual void patchType(const char*);
            virtual char* superType();
            virtual void superType(const char*);
            virtual FoamXServer::StringPairList* patchFieldTypes();
            virtual void patchFieldTypes(const FoamXServer::StringPairList&);


        //- Persist to file.
        void save(DictionaryWriter& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
