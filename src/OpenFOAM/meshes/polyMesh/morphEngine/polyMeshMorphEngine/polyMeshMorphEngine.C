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

\*---------------------------------------------------------------------------*/

#include "polyMeshMorphEngine.H"
#include "polyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshMorphEngine, 1);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Read constructor given IOobject and a polyMesh reference
Foam::polyMeshMorphEngine::polyMeshMorphEngine
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    ptrList<polyMeshModifier>(),
    regIOobject(io),
    mesh_(mesh)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        ptrList<polyMeshModifier>& modifiers = *this;

        // Read modifiers
        Istream& is = readStream(typeName);

        ptrList<entry> patchEntries(is);
        modifiers.setSize(patchEntries.size());

        forAll(modifiers, modifierI)
        {
            modifiers.hook
            (
                polyMeshModifier::New
                (
                    patchEntries[modifierI].keyword(),
                    patchEntries[modifierI].dict(),
                    modifierI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "polyMeshMorphEngine::polyMeshMorphEngine"
            "(const IOobject&, const polyMesh&)"
        );

        close();
    }
}


// Construct given size. Modifiers will be hooked later
Foam::polyMeshMorphEngine::polyMeshMorphEngine
(
    const IOobject& io,
    const polyMesh& pm,
    const label size
)
:
    ptrList<polyMeshModifier>(size),
    regIOobject(io),
    mesh_(pm)
{}


// Return a list of modifier types
Foam::wordList Foam::polyMeshMorphEngine::types() const
{
    const ptrList<polyMeshModifier>& modifiers = *this;

    wordList t(modifiers.size());

    forAll (modifiers, modifierI)
    {
        t[modifierI] = modifiers[modifierI].type();
    }

    return t;
}


// Return a list of modifier names
Foam::wordList Foam::polyMeshMorphEngine::names() const
{
    const ptrList<polyMeshModifier>& modifiers = *this;

    wordList t(modifiers.size());

    forAll (modifiers, modifierI)
    {
        t[modifierI] = modifiers[modifierI].name();
    }

    return t;
}


// Is topology change required
bool Foam::polyMeshMorphEngine::changeTopology() const
{
    // Go through all mesh modifiers and accumulate the morphing information
    const ptrList<polyMeshModifier>& morphs = *this;

    bool triggerChange = false;

    forAll (morphs, morphI)
    {
        if (morphs[morphI].active())
        {
            bool curTriggerChange = morphs[morphI].changeTopology();

            if (polyMesh::morphDebug)
            {
                Info<< "Modifier " << morphI << " named "
                    << morphs[morphI].name();
                
                if (curTriggerChange)
                {
                    Info << " morphing" << endl;
                }
                else
                {
                    Info << " unchanged" << endl;
                }
            }

            triggerChange = triggerChange || curTriggerChange;
        }
        else
        {
            if (polyMesh::morphDebug)
            {
                Info<< "Modifier " << morphI  << " named "
                    << morphs[morphI].name() << " inactive" << endl;
            }
        }
            
    }

    return triggerChange;
}


// Return topology change request
Foam::tmp<Foam::polyTopoChange>
Foam::polyMeshMorphEngine::topoChangeRequest() const
{
    // Collect changes from all modifiers
    const ptrList<polyMeshModifier>& morphs = *this;

    tmp<polyTopoChange> tref(new polyTopoChange(mesh()));
    polyTopoChange& ref = tref();

    forAll (morphs, morphI)
    {
        if (morphs[morphI].active())
        {
            morphs[morphI].setRefinement(ref);
        }
    }

    return tref;
}


// Correct polyMeshMorphEngine after moving points
void Foam::polyMeshMorphEngine::modifyMotionPoints(pointField& p) const
{
    const ptrList<polyMeshModifier>& morphs = *this;

    forAll (morphs, morphI)
    {
        if (morphs[morphI].active())
        {
            morphs[morphI].modifyMotionPoints(p);
        }
    }
}


// Force recalculation of locally stored data on topological change
void Foam::polyMeshMorphEngine::updateTopology(const mapPolyMesh& m)
{
    // Go through all mesh modifiers and accumulate the morphing information
    ptrList<polyMeshModifier>& morphs = *this;

    forAll (morphs, morphI)
    {
        morphs[morphI].updateTopology(m);
    }
}


Foam::label Foam::polyMeshMorphEngine::findModifierID
(
    const word& modName
) const
{
    const ptrList<polyMeshModifier>& morphs = *this;

    forAll (morphs, morphI)
    {
        if (morphs[morphI].name() == modName)
        {
            return morphI;
        }
    }

    // Modifier not found
    if (debug)
    {
        Info<< "label polyMeshMorphEngine::::findModifierID(const word& "
            << "modName) const"
            << "Modifier named " << modName << " not found.  "
            << "List of available modifier names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


// writeData member function required by regIOobject
bool Foam::polyMeshMorphEngine::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


Foam::Ostream& Foam::operator<<(Ostream& os, const polyMeshMorphEngine& mme)
{
    os  << mme.size() << nl << token::BEGIN_LIST;

    forAll(mme, mmeI)
    {
        mme[mmeI].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// ************************************************************************* //
