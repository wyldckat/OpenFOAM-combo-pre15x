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
    runTimeSelectionTables

Description

\*---------------------------------------------------------------------------*/

#include "token.H"

#ifndef runTimeSelectionTables_H
#define runTimeSelectionTables_H

#include "autoPtr.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define declareRunTimeSelectionTable(autoPtr,baseType,argNames,argList,parList)\
                                                                        \
    /* Construct from argList function pointer type */                  \
    typedef autoPtr<baseType> (*argNames##ConstructorPtr)argList;       \
                                                                        \
    /* Construct from argList function table type */                    \
    typedef HashTable<argNames##ConstructorPtr, word, Hash<string> >    \
        argNames##ConstructorTable;                                     \
                                                                        \
    /* Construct from argList function pointer table pointer */         \
    static argNames##ConstructorTable* argNames##ConstructorTablePtr_;  \
                                                                        \
    /* Class to add constructor from argList to table */                \
    template<class baseType##Type>                                      \
    class add##argNames##ConstructorToTable                             \
    {                                                                   \
    public:                                                             \
                                                                        \
        static autoPtr<baseType> New argList                            \
        {                                                               \
            return autoPtr<baseType>(new baseType##Type parList);       \
        }                                                               \
                                                                        \
        add##argNames##ConstructorToTable()                             \
        {                                                               \
            construct##argNames##ConstructorTables();                   \
                                                                        \
            argNames##ConstructorTablePtr_                              \
                ->insert(baseType##Type::typeName, New);                \
        }                                                               \
    };                                                                  \
                                                                        \
    /* Table Constructor called from the table add function */          \
    static void construct##argNames##ConstructorTables()


#define defineRunTimeSelectionTableConstructor(baseType,argNames)       \
                                                                        \
    /* Table Constructor called from the table add function */          \
    void baseType::construct##argNames##ConstructorTables()             \
    {                                                                   \
        static bool constructed = false;                                \
                                                                        \
        if (!constructed)                                               \
        {                                                               \
            baseType::argNames##ConstructorTablePtr_                    \
                = new baseType::argNames##ConstructorTable;             \
                                                                        \
            constructed = true;                                         \
        }                                                               \
    }

#define defineRunTimeSelectionTablePtr(baseType,argNames)               \
                                                                        \
    /* Define the constructor function table */                         \
    baseType::argNames##ConstructorTable*                               \
        baseType::argNames##ConstructorTablePtr_ = NULL

#define defineTemplateRunTimeSelectionTablePtr(baseType,argNames)       \
                                                                        \
    /* Define the constructor function table */                         \
    typename baseType::argNames##ConstructorTable*                      \
        baseType::argNames##ConstructorTablePtr_ = NULL

#define defineRunTimeSelectionTable(baseType,argNames)                  \
                                                                        \
    defineRunTimeSelectionTablePtr(baseType,argNames);                  \
    defineRunTimeSelectionTableConstructor(baseType,argNames)

#define defineTemplateRunTimeSelectionTable(baseType,argNames)          \
                                                                        \
    template<>                                                          \
    defineRunTimeSelectionTablePtr(baseType,argNames);                  \
    template<>                                                          \
    defineRunTimeSelectionTableConstructor(baseType,argNames)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
