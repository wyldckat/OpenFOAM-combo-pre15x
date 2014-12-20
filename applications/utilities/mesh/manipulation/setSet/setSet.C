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
    Manipulate a cell/face/point set interactively.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "topoSetSource.H"
#include "topoSet.H"
#include "IStringStream.H"
#include "topoSet.H"
#include "cellSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "demandDrivenData.H"
#include "writeFaceSet.H"
#include "writePointSet.H"

#include <stdio.h>

#if defined(READLINE)
#include <readline/readline.h>
#include <readline/history.h>
#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#if defined(READLINE)
static const char* historyFile = ".setSet";
#endif

Istream& selectStream(Istream* is0Ptr, Istream* is1Ptr)
{
    if (is0Ptr)
    {
        return *is0Ptr;
    }
    else if (is1Ptr)
    {
        return *is1Ptr;
    }
    else
    {
        FatalErrorIn("selectStream(Istream*, Istream*)")
            << "No valid stream opened" << abort(FatalError);

        return *is0Ptr;
    }
}

// Copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const topoSet& fromSet,
    const word& toName
)
{
    if (fromSet.size() > 0)
    {
        Info<< "    Backing up " << fromName << " into " << toName << endl;

        topoSet backupSet(mesh, toName, fromSet);

        backupSet.write();
    }
}


// Read and copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const word& toName
)
{
    topoSet fromSet(mesh, fromName, IOobject::READ_IF_PRESENT);

    backup(mesh, fromName, fromSet, toName);
}


// Write set to VTK readable files
void writeVTK
(
    const polyMesh& mesh,
    const topoSet& currentSet,
    const fileName& vtkName
)
{
    if (typeid(currentSet) == typeid(faceSet))
    {
        writeFaceSet
        (
            true,
            mesh,
            currentSet,
            mesh.time().path()/vtkName
        );
    }
    else if (typeid(currentSet) == typeid(cellSet))
    {
        faceSet cellFaces(mesh, "cellFaces", currentSet.size());

        forAllConstIter(cellSet, currentSet, iter)
        {
            const cell& cFaces = mesh.cells()[iter.key()];

            forAll(cFaces, i)
            {
                cellFaces.insert(cFaces[i]);
            }
        }

        writeFaceSet
        (
            true,
            mesh,
            cellFaces,
            mesh.time().path()/vtkName
        );
    }
    else if (typeid(currentSet) == typeid(pointSet))
    {
        writePointSet
        (
            true,
            mesh,
            currentSet,
            mesh.time().path()/vtkName
        );        
    }
    else
    {
        Warning<< "Don't know how to handle set of type " << currentSet.type()
            << endl;
    }
}


void printHelp(Ostream& os)
{
    os  << "Please type 'help', 'quit' or a set command after prompt." << endl
        << endl
        << "A set command should be of the following form" << endl
        << endl
        << "    cellSet|faceSet|pointSet <setName> <action> <source>"
        << endl << endl
        << "The <action> is one of" << endl
        << "    list            - prints the contents of the set" << endl
        << "    clear           - clears the set" << endl
        << "    invert          - inverts the set" << endl
        << "    new <source>    - sets to set to the source set" << endl
        << "    add <source>    - adds all elements from the source set" << endl
        << "    delete <source> - deletes      ,," << endl
        << "    subset <source> - combines current set with the source set"
        << endl
        << endl
        << "The sources come in various forms. Type a wrong source"
        << " to see all the types available." << endl
        << endl
        << "Example: pick up all cells connected by point or face to patch"
        << " movingWall" << endl
        << endl
        << "Pick up all faces of patch:" << endl
        << "    faceSet f0 new patchToFace movingWall" << endl
        << "Add faces 0,1,2:" << endl
        << "    faceSet f0 add labelToFace (0 1 2)" << endl
        << "Pick up all points used by faces in faceSet f0:" << endl
        << "    pointSet p0 new faceToPoint f0 all" << endl
        << "Pick up cell which has any face in f0:" << endl
        << "    cellSet c0 new faceToCell f0 any" << endl
        << "Add cells which have any point in p0:" << endl
        << "    cellSet c0 add pointToCell p0 any" << endl
        << "List set:" << endl
        << "    cellSet c0 list" << endl
        << endl;
}



// Read command and execute. Return true if ok, false otherwise.
bool doCommand
(
    const polyMesh& mesh,
    const word& setType,
    const word& setName,
    const word& actionName,
    const bool writeVTKFile,
    Istream& is
)
{

    // Get some size estimate for set.
    label typSize = max(mesh.nCells(), max(mesh.nFaces(), mesh.nPoints()))/10;


    bool error = false;

    // Set to work on
    topoSet* currentSetPtr(NULL);

    word sourceType;

    try
    {
        topoSetSource::setAction action =
            topoSetSource::toAction(actionName);


        IOobject::readOption r;

        if
        (
            (action == topoSetSource::NEW)
         || (action == topoSetSource::CLEAR)
        )
        {
            r = IOobject::NO_READ;

            backup(mesh, setName, setName + "_old");

            currentSetPtr =
                topoSet::New
                (
                    setType,
                    mesh,
                    setName,
                    typSize
                ).ptr();
        }
        else
        {
            r = IOobject::MUST_READ;

            if (exists(mesh.time().path()/topoSet::localPath(mesh, setName)))
            {
                currentSetPtr =
                    topoSet::New
                    (
                        setType,
                        mesh,
                        setName,
                        r
                    ).ptr();

                topoSet& currentSet = *currentSetPtr;

                currentSet.resize(max(currentSet.size(), typSize));               
            }
        }

        if (!currentSetPtr)
        {
            Info<< "    Cannot construct/load set "
                << topoSet::localPath(mesh, setName) << endl;

            error = true;
        }
        else
        {
            topoSet& currentSet = *currentSetPtr;

            Info<< "    Set:" << currentSet.name()
                << "  Size:" << currentSet.size()
                << "  Action:" << actionName
                << endl;

            if ((r == IOobject::MUST_READ) && (action != topoSetSource::LIST))
            {
                // currentSet has been read so can make copy.
                backup(mesh, setName, currentSet, setName + "_old");
            }

            if (action == topoSetSource::CLEAR)
            {
                // Already handled above by not reading
            }
            else if (action == topoSetSource::INVERT)
            {
                currentSet.invert(currentSet.maxSize(mesh));
            }
            else if (action == topoSetSource::LIST)
            {
                currentSet.writeDebug(Info, mesh, 100);
                Info<< endl;
            }
            else if (action == topoSetSource::SUBSET)
            {
                if (is >> sourceType)
                {
                    autoPtr<topoSetSource> setSource
                    (
                        topoSetSource::New
                        (
                            sourceType,
                            mesh,
                            is
                        )
                    );

                    // Backup current set.
                    topoSet oldSet
                    (
                        mesh,
                        currentSet.name() + "_old2",
                        currentSet
                    );

                    currentSet.clear();
                    currentSet.resize(oldSet.size());
                    setSource().applyToSet(topoSetSource::NEW, currentSet);

                    // Combine new value of currentSet with old one.
                    currentSet.subset(oldSet);
                }
            }
            else
            {        
                if (is >> sourceType)
                {
                    autoPtr<topoSetSource> setSource
                    (
                        topoSetSource::New
                        (
                            sourceType,
                            mesh,
                            is
                        )
                    );

                    setSource().applyToSet(action, currentSet);
                }
            }


            if (action != topoSetSource::LIST)
            {
                if (writeVTKFile)
                {
                    mkDir(mesh.time().path()/"VTK"/currentSet.name());

                    fileName vtkName
                    (
                        "VTK"/currentSet.name()/currentSet.name()
                      + "_"
                      + name(mesh.time().timeIndex())
                      + ".vtk"
                    );

                    Info<< "    Writing " << currentSet.name()
                        << " (size " << currentSet.size() << ") to "
                        << currentSet.instance()/currentSet.local()
                           /currentSet.name()
                        << " and to vtk file " << vtkName << endl << endl;

                    currentSet.write();

                    writeVTK(mesh, currentSet, vtkName);
                }
                else
                {
                    Info<< "    Writing " << currentSet.name()
                        << " (size " << currentSet.size() << ") to "
                        << currentSet.instance()/currentSet.local()
                           /currentSet.name() << endl << endl;

                    currentSet.write();
                }
            }
        }
    }
    catch (Foam::IOerror& fIOErr)
    {
        error = true;

        Info<< fIOErr.message().c_str() << endl;

        if (sourceType.size() != 0)
        {
            Info << topoSetSource::usage(sourceType).c_str();
        }
    }
    catch (Foam::error& fErr)
    {
        error = true;

        Info<< fErr.message().c_str() << endl;

        if (sourceType.size() != 0)
        {
            Info << topoSetSource::usage(sourceType).c_str();
        }
    }

    deleteDemandDrivenData(currentSetPtr);

    return error;
}


enum commandStatus
{
    QUIT,
    INVALID,
    VALID
};


commandStatus parseType(const word& setType)
{
    if (setType.size() == 0)
    {
        Info<< "Type 'help' for usage information" << endl;

        return INVALID;
    }
    else if (setType == "help")
    {
        printHelp(Info);

        return INVALID;
    }
    else if (setType == "quit")
    {
        Info<< "Quitting ..." << endl;

        return QUIT;
    }
    else if
    (
        setType == "cellSet"
     || setType == "faceSet"
     || setType == "pointSet"
    )
    {
        return VALID;
    }
    else
    {
        SeriousError
            << "Illegal set type " << setType << endl
            << "Should be one of 'cellSet' 'faceSet' 'pointSet'"
            << endl;

        return INVALID;
    }
}


commandStatus parseAction(const word& actionName)
{
    commandStatus stat = INVALID;

    if (actionName.size() != 0)
    {
        try
        {
            (void)topoSetSource::toAction(actionName);

            stat = VALID;
        }
        catch (Foam::IOerror& fIOErr)
        {
            stat = INVALID;
        }
        catch (Foam::error& fErr)
        {
            stat = INVALID;
        }
    }
    return stat;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("noVTK", "");
    argList::validOptions.insert("batch", "file");
#   include "addTimeOptions.H"

#   include "setRootCase.H"
#   include "createTime.H"

    bool writeVTK = !args.options().found("noVTK");

    // Get times list
    instantList Times = runTime.times();

    label startTime = Times.size()-1;
    label endTime = Times.size();

#   include "checkTimeOption.H"
#   include "checkLatestTimeOption.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createPolyMesh.H"

    std::ifstream* fileStreamPtr(NULL);

    if (args.options().found("batch"))
    {
        fileName batchFile(args.options()["batch"]);

        Info<< "Reading commands from file " << batchFile << endl;

        if (!exists(batchFile))
        {
            FatalErrorIn(args.executable())
                << "Cannot open file " << batchFile << exit(FatalError);
        }

        fileStreamPtr = new std::ifstream(batchFile.c_str());
    }
#if defined(READLINE)
    else if (!read_history(historyFile))
    {
        Info<< "Successfully read history from " << historyFile << endl;

        string msg =
          "# Time:" + runTime.timeName()
        + "  cells:" + Foam::name(mesh.nCells())
        + "  faces:" + Foam::name(mesh.nFaces())
        + "  points:" + Foam::name(mesh.nPoints());
    
        add_history(msg.c_str());
    }        
#endif

    Info<< "Please type 'help', 'quit' or a set command after prompt." << endl;

    bool ok = false;

    FatalError.throwExceptions();
    FatalIOError.throwExceptions();

    while (1)
    {
        string rawLine;

        // Type: cellSet, faceSet, pointSet
        word setType;
        // Name of destination set.
        word setName;
        // Action (new, invert etc.)
        word actionName;

        commandStatus stat = INVALID;

        if (fileStreamPtr)
        {
            if (!fileStreamPtr->good())
            {
                Info<< "End of batch file" << endl;
                break;
            }

            std::getline(*fileStreamPtr, rawLine);

            if (rawLine.size() > 0)
            {
                Info<< "Doing:" << rawLine << endl;
            }
        }
        else
        {
#           if defined(READLINE)
            {
                char* linePtr = readline("readline>");

                rawLine = string(linePtr);

                if (*linePtr)
                {
                    add_history(linePtr);
                    write_history(historyFile);
                }

                free(linePtr);   // readline uses malloc, not new.
            }
#           else
            {
                Info<< "Command>" << flush;
                std::getline(std::cin, rawLine);
            }
#           endif
        }

        if (rawLine.size() == 0 || rawLine[0] == '#')
        {
            continue;
        }

        IStringStream is(rawLine + ' ');

        // Type: cellSet, faceSet, pointSet
        is >> setType;

        stat = parseType(setType);

        if (stat == VALID)
        {
            if (is >> setName)
            {
                if (is >> actionName)
                {
                    stat = parseAction(actionName);
                }
            }
        }

        ok = false;

        if (stat == QUIT)
        {
            break;
        }
        else if (stat == VALID)
        {
            ok = doCommand(mesh, setType, setName, actionName, writeVTK, is);
        }

    } while (ok);


    if (fileStreamPtr)
    {
        delete fileStreamPtr;
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
