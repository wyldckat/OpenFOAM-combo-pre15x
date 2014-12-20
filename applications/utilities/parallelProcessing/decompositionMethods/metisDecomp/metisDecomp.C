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

#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "IFstream.H"

extern "C"
{
#   include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisDecomp,
        dictionaryMesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisDecomp::metisDecomp
(
    const dictionary& decompositionDict,
    const primitiveMesh& mesh
)
:
    decompositionMethod(decompositionDict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisDecomp::decompose(const pointField& points)
{

    //
    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    //

    labelList xadj(mesh_.nCells()+1);

    labelList adjncy(2*mesh_.nInternalFaces());

    // Fill in xadj

    label freeAdj = 0;

    for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
    {
        xadj[cellI] = freeAdj;

        const labelList& cFaces = mesh_.cells()[cellI];

        forAll(cFaces, i)
        {
            label faceI = cFaces[i];

            if (mesh_.isInternalFace(faceI))
            {
                freeAdj++;
            }
        }
    }
    xadj[mesh_.nCells()] = freeAdj;


    // Fill in adjncy

    labelList nFacesPerCell(mesh_.nCells(), 0);

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        adjncy[xadj[own] + nFacesPerCell[own]++] = nei;
        adjncy[xadj[nei] + nFacesPerCell[nei]++] = own;
    }


    // C style numbering
    int numFlag = 0;

    // decomposition options. 0 = use defaults
    labelList options(5, 0);

    // processor weights. Use even weighting
    Field<floatScalar> processorWeights(nProcessors_, 1.0/nProcessors_);

    // cell weights (so on the vertices of the dual)
    labelList cellWeights;

    // face weights (so on the edges of the dual)
    labelList faceWeights;

    Info<< endl;

    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        dictionary metisDecompCoeffs
        (
            decompositionDict_.subDict("metisCoeffs")
        );

        if (metisDecompCoeffs.found("processorWeights"))
        {
            metisDecompCoeffs.lookup("processorWeights") >> processorWeights;

            if (processorWeights.size() != nProcessors_)
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of processor weights "
                    << processorWeights.size()
                    << " does not equal number of domains " << nProcessors_
                    << exit(FatalError);
            }
        }

        if (metisDecompCoeffs.found("cellWeightsFile"))
        {
            Info<< "metisDecomp : Using cell-based weights." << endl;

            fileName cellWeightsFile
            (
                metisDecompCoeffs.lookup("cellWeightsFile")
            );

            IFstream decompStream(cellWeightsFile);

            if (!decompStream)
            {
                FatalIOErrorIn("manualDecomp::decompose()", decompStream)
                    << "Cannot read cell weights data file "
                    << cellWeightsFile << "." << endl
                    << exit(FatalIOError);
            }

            decompStream >> cellWeights;

            if (cellWeights.size() != mesh_.nCells())
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of cell weights " << cellWeights.size()
                    << " does not equal number of cells " << mesh_.nCells()
                    << exit(FatalError);
            }
        }

        if (metisDecompCoeffs.found("faceWeightsFile"))
        {
            Info<< "metisDecomp : Using face-based weights." << endl;

            fileName faceWeightsFile
            (
                metisDecompCoeffs.lookup("faceWeightsFile")
            );

            IFstream decompStream(faceWeightsFile);

            if (!decompStream)
            {
                FatalIOErrorIn("manualDecomp::decompose()", decompStream)
                    << "Cannot read face weights data file "
                    << faceWeightsFile << "." << endl
                    << exit(FatalIOError);
            }

            labelList weights(decompStream);

            if (weights.size() != mesh_.nInternalFaces())
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of face weights " << weights.size()
                    << " does not equal number of internal faces "
                    << mesh_.nInternalFaces()
                    << exit(FatalError);
            }

            // Assume symmetric weights. Keep same ordering as adjncy.
            faceWeights.setSize(2*mesh_.nInternalFaces());

            labelList nFacesPerCell(mesh_.nCells(), 0);

            for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
            {
                label w = weights[faceI];

                label own = mesh_.faceOwner()[faceI];
                label nei = mesh_.faceNeighbour()[faceI];

                faceWeights[xadj[own] + nFacesPerCell[own]++] = w;
                faceWeights[xadj[nei] + nFacesPerCell[nei]++] = w;
            }
        }

        if (metisDecompCoeffs.found("options"))
        {
            metisDecompCoeffs.lookup("options") >> options;

            if (options.size() != 5)
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 5"
                    << abort(FatalError);
            }

            Info<< "Using Metis options     " << options
                << endl << endl;
        }
    }

    processorWeights /= sum(processorWeights);


    int numCells = mesh_.nCells();

    // output: cell -> processor addressing
    labelList finalDecomp(mesh_.nCells());

    // output: number of cut edges
    int edgeCut = 0;

    // Vertex weight info
    int wgtFlag = 0;
    label* vwgtPtr = NULL;
    label* adjwgtPtr = NULL;

    if (cellWeights.size() > 0)
    {
        vwgtPtr = cellWeights.begin();
        wgtFlag += 2;       // Weights on vertices
    }
    if (faceWeights.size() > 0)
    {
        adjwgtPtr = faceWeights.begin();
        wgtFlag += 1;       // Weights on edges
    }

        
    METIS_WPartGraphKway
    (
        &numCells,         // num vertices in graph
        xadj.begin(),      // indexing into adjncy
        adjncy.begin(),    // neighbour info
        vwgtPtr,           // vertexweights
        adjwgtPtr,         // no edgeweights
        &wgtFlag,
        &numFlag,
        &nProcessors_,
        processorWeights.begin(),
        options.begin(),
        &edgeCut,
        finalDecomp.begin()
    );

    return finalDecomp;
}


// ************************************************************************* //
