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

#include "topoFvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoFvMesh> Foam::topoFvMesh::New
(
    const Foam::IOobject& io
)
{
    word topoFvMeshTypeName;

    // Enclose the creation of the meshMotionDict to ensure it is
    // deleted before the topoFvMesh is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary meshMotionDictDict
        (
            IOobject
            (
                "meshMotionDict",
                io.time().constant(),
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        meshMotionDictDict.lookup("topoFvMesh") >> topoFvMeshTypeName;
    }

    Info<< "Selecting topoFvMesh " << topoFvMeshTypeName << endl;

    IOobjectConstructorTable::iterator cstrIter =
        IOobjectConstructorTablePtr_->find(topoFvMeshTypeName);

    if (cstrIter == IOobjectConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "topoFvMesh::New(const IOobject&)"
        )   << "Unknown topoFvMesh type " << topoFvMeshTypeName
            << endl << endl
            << "Valid topoFvMesh types are :" << endl
            << IOobjectConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<topoFvMesh>(cstrIter()(io));
}


// ************************************************************************* //
