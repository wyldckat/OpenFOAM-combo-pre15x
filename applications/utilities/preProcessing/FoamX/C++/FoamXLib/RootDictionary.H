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
    IDictionaryEntryImpl

Description

SourceFiles
    RootDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef RootDictionary_H
#define RootDictionary_H

// FoamX header files.
#include "IDictionaryEntryImpl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class RootDictionary Declaration
\*---------------------------------------------------------------------------*/

class RootDictionary
:
    public IDictionaryEntryImpl
{
    // Private data

        Foam::fileName caseRoot_;
        Foam::fileName caseName_;

public:

    // Constructors

        RootDictionary
        (
            FoamXServer::ITypeDescriptor_ptr typeDesc,
            const Foam::fileName& caseRoot,
            const Foam::fileName& caseName
        );


    // Destructor

        virtual ~RootDictionary();


    // Member Functions

        // Overridden IDictionaryEntry load
        virtual void load(Foam::Istream&);

        // Overridden IDictionaryEntry save
        virtual void save();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
