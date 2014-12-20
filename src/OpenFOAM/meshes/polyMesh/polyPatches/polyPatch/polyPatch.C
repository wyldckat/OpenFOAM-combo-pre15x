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

#include "polyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "SubField.H"
#include "entry.H"
#include "dictionary.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyPatch, 0);

defineRunTimeSelectionTable(polyPatch, word);
defineRunTimeSelectionTable(polyPatch, Istream);
defineRunTimeSelectionTable(polyPatch, dictionary);

addToRunTimeSelectionTable(polyPatch, polyPatch, word);
addToRunTimeSelectionTable(polyPatch, polyPatch, Istream);
addToRunTimeSelectionTable(polyPatch, polyPatch, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void polyPatch::movePoints(const pointField& p)
{
    primitivePatch::movePoints(p);
}

void polyPatch::updateTopology()
{
    deleteDemandDrivenData(mePtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(name, index),
    primitivePatch
    (
        faceSubList(bm.mesh().allFaces(), size, start),
        bm.mesh().allPoints()
    ),
    start_(start),
    boundaryMesh_(bm),
    mePtr_(NULL)
{}


// Construct from Istream
polyPatch::polyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(is, index),
    primitivePatch
    (
        faceSubList(bm.mesh().allFaces(), 0, 0),
        bm.mesh().allPoints()
    ),
    start_(-1),
    boundaryMesh_(bm),
    mePtr_(NULL)
{
    label size = readLabel(is);
    start_ = readLabel(is);
    UList<face>::operator=
    (
        faceSubList(bm.mesh().allFaces(), size, start_)
    );
}


// Construct from dictionary
polyPatch::polyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(name, dict, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().allFaces(),
            readLabel(dict.lookup("nFaces")),
            readLabel(dict.lookup("startFace"))
        ),
        bm.mesh().allPoints()
    ),
    start_(readLabel(dict.lookup("startFace"))),
    boundaryMesh_(bm),
    mePtr_(NULL)
{}


// Construct as copy. Resets the reference to polyBoundaryMesh but the start
// and size of the patch remain the same
// Used in clone function
polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(pp),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().allFaces(),
            pp.size(),
            pp.start()
        ),
        bm.mesh().allPoints()
    ),
    start_(pp.start()),
    boundaryMesh_(bm),
    mePtr_(NULL)
{}


// Construct as copy. Resets the reference to face list and polyBoundaryMesh
// Start and size of the patch can be changed as well
// Used in clone function
polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    patchIdentifier(pp, index),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().allFaces(),
            newSize,
            newStart
        ),
        bm.mesh().allPoints()
    ),
    start_(newStart),
    boundaryMesh_(bm),
    mePtr_(NULL)
{}


polyPatch::polyPatch(const polyPatch& p)
:
    patchIdentifier(p),
    primitivePatch(p),
    start_(p.start_),
    boundaryMesh_(p.boundaryMesh_),
    mePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPatch::~polyPatch()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const polyBoundaryMesh& polyPatch::boundaryMesh() const
{
    return boundaryMesh_;
}


const pointField& polyPatch::allPoints() const
{
    return boundaryMesh_.mesh().allPoints();
}


bool polyPatch::constraintType(const word& pt)
{
    return fvPatchField<scalar>::patchConstructorTablePtr_->found(pt);
}


wordList polyPatch::constraintTypes()
{
    wordList cTypes(dictionaryConstructorTablePtr_->size());

    label i = 0;

    for
    (
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->begin();
        cstrIter != dictionaryConstructorTablePtr_->end();
        ++cstrIter
    )
    {
        if (constraintType(cstrIter.key()))
        {
            cTypes[i++] = cstrIter.key();
        }
    }

    cTypes.setSize(i);

    return cTypes;
}


const vectorField::subField polyPatch::faceCentres() const
{
    return patchSlice(boundaryMesh().mesh().faceCentres());
}


const vectorField::subField polyPatch::faceAreas() const
{
    return patchSlice(boundaryMesh().mesh().faceAreas());
}


// Return the patch face neighbour cell centres
tmp<vectorField> polyPatch::faceCellCentres() const
{
    tmp<vectorField> tcc(new vectorField(size()));
    vectorField& cc = tcc();

    // get reference to global cell centres
    const vectorField& gcc = boundaryMesh_.mesh().cellCentres();
    const labelList::subList cellLabels = faceCells();

    forAll (cellLabels, faceI)
    {
        cc[faceI] = gcc[cellLabels[faceI]];
    }

    /* This should not be used because it's inconsistent with the
    // mesh().cellCentres() used elsewhere
    const pointField& points = allPoints();
    const cellList& cells = boundaryMesh_.mesh().cells();
    const faceList& faces = boundaryMesh_.mesh().faces();
    const labelList::subList cellLabels = faceCells();

    forAll (cellLabels, faceI)
    {
        cc[faceI] = cells[cellLabels[faceI]].centre(points, faces);
    }
    */

    return tcc;
}


const labelList::subList polyPatch::faceCells() const
{
    return patchSlice(boundaryMesh().mesh().faceOwner());
}


const labelList& polyPatch::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_ =
            new labelList
            (
                primitivePatch::meshEdges
                (
                    boundaryMesh().mesh().edges(),
                    boundaryMesh().mesh().cellEdges(),
                    faceCells()
                )
            );
    }

    return *mePtr_;
}


void polyPatch::clearAddressing()
{
    deleteDemandDrivenData(mePtr_);
}


void polyPatch::write(Ostream& os) const
{
    os  << nl << type()
        << static_cast<const patchIdentifier&>(*this) << endl
        << this->size() << tab << start_;
}


void polyPatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    patchIdentifier::writeDict(os);

    os  << "    nFaces " << this->size() << token::END_STATEMENT << nl
        << "    startFace " << start() << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void polyPatch::operator=(const polyPatch& p)
{
    clearAddressing();

    patchIdentifier::operator=(p);
    primitivePatch::operator=(p);
    start_ = p.start_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const polyPatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const polyPatch& p");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
