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

Description

\*---------------------------------------------------------------------------*/

#include "IOField.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<IOField<Type> > reconstructLagrangianField
(
    const Time& runTime,
    PtrList<Time>& databases,
    const IOobject& fieldIoObject
)
{
    tmp<IOField<Type> > tfield
    (
        new IOField<Type>
        (
            IOobject
            (
                fieldIoObject.name(),
                runTime.timeName(),
                "lagrangian",
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Field<Type>(0)
        )
    );
    Field<Type>& field = tfield();

    forAll(databases, i)
    {
        IOobject fieldIOobject
        ( 
            fieldIoObject.name(),
            databases[i].timeName(),
            "lagrangian",
            databases[i],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldIOobject.headerOk())
        {
            IOField<Type> fieldi(fieldIOobject);

            label offset = field.size();
            field.setSize(offset + fieldi.size());
            
            forAll(fieldi, j)
            {
                field[offset + j] = fieldi[j];
            }
        }
    }

    return tfield;
}


template<class Type>
void reconstructLagrangianFields
(
    const Time& runTime,
    PtrList<Time>& databases,
    const IOobjectList& objects
)
{
    word fieldClassName(IOField<Type>::typeName);

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing lagrangian "
            << fieldClassName << "s\n" << endl;

        for
        (
            IOobjectList::iterator fieldIter = fields.begin();
            fieldIter != fields.end();
            ++fieldIter
        )
        {
            Info<< "        " << fieldIter()->name() << endl;
            reconstructLagrangianField<Type>
            (
                runTime,
                databases,
                *fieldIter()
            )().write();
        }

        Info<< endl;
    }
}


// ************************************************************************* //
