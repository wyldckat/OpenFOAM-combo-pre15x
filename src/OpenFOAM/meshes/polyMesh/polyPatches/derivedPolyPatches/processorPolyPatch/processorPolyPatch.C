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

#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorPolyPatch, 0);

addToRunTimeSelectionTable(polyPatch, processorPolyPatch, Istream);
addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Construct from components
processorPolyPatch::processorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo
)
:
    coupledPolyPatch(name, size, start, index, bm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo)
{}


// Construct from Istream
processorPolyPatch::processorPolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(is, index, bm),
    myProcNo_(readLabel(is)),
    neighbProcNo_(readLabel(is))
{}


// Construct from dictionary
processorPolyPatch::processorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo")))
{}


//- Construct as copy, resetting the boundary mesh
processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_)
{}


//- Construct as copy, resetting the face list and boundary mesh data
processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorPolyPatch::~processorPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void processorPolyPatch::initGeometry()
{
    if (Pstream::parRun())
    {
        OPstream toNeighbProc(neighbProcNo(), 3*size()*sizeof(vector));

        toNeighbProc
            << faceCentres()
            << faceAreas()
            << faceCellCentres();
    }
}


void processorPolyPatch::calcGeometry()
{
    if (Pstream::parRun())
    {
        {
            IPstream fromNeighbProc(neighbProcNo(), 3*size()*sizeof(vector));
            fromNeighbProc
                >> neighbFaceCentres_
                >> neighbFaceAreas_
                >> neighbFaceCellCentres_;
        }

        scalarField magSf = mag(faceAreas());

        forAll(magSf, facei)
        {
            scalar nmagSf = mag(neighbFaceAreas_[facei]);
            scalar avSf = (magSf[facei] + nmagSf)/2.0;

            if (avSf > VSMALL && mag(magSf[facei] - nmagSf)/avSf > 1e-6)
            {
                FatalErrorIn
                (
                    "processorFvPatch::makeWeights(scalarField& w) const"
                )   << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf[facei] - nmagSf)/avSf
                    << "% -- possible face ordering problem"
                    << exit(FatalError);
            }
        }

        calcTransformTensors
        (
            faceCentres(),
            neighbFaceCentres_,
            faceNormals(),
            neighbFaceAreas_/(mag(neighbFaceAreas_)+VSMALL)
        );
    }
}


void processorPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    initGeometry();
}


void processorPolyPatch::movePoints(const pointField&)
{
    calcGeometry();
}


// Write
void processorPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);

    os  << nl << myProcNo_ << token::SPACE << neighbProcNo_;
}

void processorPolyPatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;
    patchIdentifier::writeDict(os);
    os  << "    nFaces " << size() << token::END_STATEMENT << nl
        << "    startFace " << start() << token::END_STATEMENT << nl
        << "    myProcNo " << myProcNo_ << token::END_STATEMENT << nl
        << "    neighbProcNo " << neighbProcNo_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
