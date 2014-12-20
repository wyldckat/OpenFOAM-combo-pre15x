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
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace debug
{

dictionary* debugSwitchesPtr_ = NULL;
dictionary* infoSwitchesPtr_ = NULL;
dictionary* optimisationSwitchesPtr_ = NULL;

dictionary* readSwitchSet
(
    const fileName& dictFileName,
    const char* switchSet
)
{
    IFstream dictFile(dictFileName);

    if (!dictFile.good())
    {
        cerr<< "readSwitchSet(const fileName& dictFileName, "
               "const char* switchSet) : "
            << "cannot open essential file " << dictFileName.c_str()
            << std::endl << std::endl;
        ::exit(1);
    }

    dictionary controlDict(dictFile);

    if (!controlDict.found(switchSet))
    {
        cerr<< "readSwitchSet(const fileName& dictFileName, "
               "const char* switchSet) : "
            << "cannot find " <<  switchSet
            << " entry in file " << dictFileName.c_str()
            << std::endl << std::endl;
        ::exit(1);
    }

    return new dictionary(controlDict.subDict(switchSet));
}


typedef dictionary* dictionaryPtr;

dictionary& switchSet
(
    const char* switchSetName,
    dictionaryPtr& switchSetPtr
)
{
    if (!switchSetPtr)
    {
        switchSetPtr = readSwitchSet(dotFoam("controlDict"), switchSetName);
    }

    return *switchSetPtr;
}


int debugSwitch
(
    dictionary& switchSet,
    const char* switchName,
    const int defaultValue
)
{
    if (switchSet.found(switchName))
    {
        return readInt(switchSet.lookup(switchName));
    }
    else
    {
        switchSet.add(switchName, defaultValue);
        return defaultValue;
    }
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary& debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


int debug::debugSwitch(const char* switchName, const int defaultValue)
{
    return debugSwitch
    (
        debugSwitches(),
        switchName,
        defaultValue
    );
}


dictionary& debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


int debug::infoSwitch(const char* switchName, const int defaultValue)
{
    return debugSwitch
    (
        infoSwitches(),
        switchName,
        defaultValue
    );
}


dictionary& debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


int debug::optimisationSwitch(const char* switchName, const int defaultValue)
{
    return debugSwitch
    (
        optimisationSwitches(),
        switchName,
        defaultValue
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
