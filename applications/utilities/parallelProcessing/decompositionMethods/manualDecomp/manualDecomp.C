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
    Decomposition given a cell-to-processor association in a file

\*---------------------------------------------------------------------------*/

#include "manualDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        manualDecomp,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        manualDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualDecomp::manualDecomp(const dictionary& decompositionDict)
:
    decompositionMethod(decompositionDict),
    decompDataFile_
    (
        decompositionDict.subDict(word(decompositionDict.lookup("method"))
      + "Coeffs").lookup("dataFile")
    )
{}


Foam::manualDecomp::manualDecomp
(
    const dictionary& decompositionDict, 
    const primitiveMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    decompDataFile_
    (
        decompositionDict.subDict(word(decompositionDict.lookup("method"))
      + "Coeffs").lookup("dataFile")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::manualDecomp::decompose(const pointField& points)
{
    IFstream decompStream(decompDataFile_);

    if (!decompStream)
    {
        FatalIOErrorIn("manualDecomp::decompose()", decompStream)
            << "Cannot read manual decomposition data file "
            << decompDataFile_ << "." << endl
            << exit(FatalIOError);
    }

    labelList finalDecomp(decompStream);

    // check if the final decomposition is OK

    if (finalDecomp.size() != points.size())
    {
        FatalErrorIn("manualDecomp::decompose()")
            << "Size of decomposition list does not correspond "
            << "to the number of points.  Size: "
            << finalDecomp.size() << " Number of points: "
            << points.size()
            << ".\n" << "Manual decomposition data read from file "
            << decompDataFile_ << "." << endl
            << exit(FatalError);
    }

    if (min(finalDecomp) < 0 || max(finalDecomp) > nProcessors_ - 1)
    {
        FatalErrorIn("labelList manualDecomp::decompose()")
            << "According to the decomposition, cells assigned to "
            << "impossible processor numbers.  Min processor = "
            << min(finalDecomp) << " Max processor = " << max(finalDecomp)
            << ".\n" << "Manual decomposition data read from file "
            << decompDataFile_ << "." << endl
            << exit(FatalError);
    }

    return finalDecomp;
}


// ************************************************************************* //
