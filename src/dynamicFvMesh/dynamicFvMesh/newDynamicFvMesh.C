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

\*---------------------------------------------------------------------------*/

#include "dynamicFvMesh.H"
#include "Time.H"

#include <dlfcn.h>

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicFvMesh> Foam::dynamicFvMesh::New(const IOobject& io)
{
    // Enclose the creation of the dynamicMesh to ensure it is
    // deleted before the dynamicFvMesh is created otherwise the dictionary
    // is entered in the database twice
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            io.time().constant(),
            io.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    word dynamicFvMeshTypeName(dynamicMeshDict.lookup("dynamicFvMesh"));

    Info<< "Selecting dynamicFvMesh " << dynamicFvMeshTypeName << endl;

    if (dynamicMeshDict.found("dynamicFvMeshLib"))
    {
        string dynamicFvMeshLibName =
            dynamicMeshDict.lookup("dynamicFvMeshLib");

        if (dynamicFvMeshLibName.size())
        {
            label nSolvers = 0;

            if (IOobjectConstructorTablePtr_)
            {
                nSolvers = IOobjectConstructorTablePtr_->size();
            }

            void* dynamicFvMeshLibPtr =
                dlopen(dynamicFvMeshLibName.c_str(), RTLD_LAZY|RTLD_GLOBAL);

            if (!dynamicFvMeshLibPtr)
            {
                FatalErrorIn
                (
                    "dynamicFvMesh::New(const IOobject&)"
                )   << "could not load " << dlerror()
                    << exit(FatalError);
            }
            else if 
            (
               !IOobjectConstructorTablePtr_
             || IOobjectConstructorTablePtr_->size() <= nSolvers
            )
            {
                WarningIn
                (
                    "dynamicFvMesh::New(const IOobject&)"
                )   << "library " << dynamicFvMeshLibName
                    << " did not introduce any new solvers"
                    << endl << endl;
            }
        }
    }


    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(dynamicFvMeshTypeName);

    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dynamicFvMesh::New(const IOobject&)"
        )   << "Unknown dynamicFvMesh type " << dynamicFvMeshTypeName
            << endl << endl
            << "Valid dynamicFvMesh types are :" << endl
            << IOobjectConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<dynamicFvMesh>(cstrIter()(io));
}


// ************************************************************************* //
