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
    label numCells = points.size();

    //
    // Make Metis CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph)
    //   xadj(celli) : start of information in adjncy for celli
    //

    labelList xadj(numCells+1);

    labelList adjncy(2*mesh_.nInternalFaces());

    label freeAdj = 0;

    const labelListList& cellsList = mesh_.cellCells();

    forAll(cellsList, celli)
    {
        xadj[celli] = freeAdj;

        const labelList& neighbours = cellsList[celli];

        forAll(neighbours, neighbouri)
        {
            adjncy[freeAdj++] = neighbours[neighbouri];
        }
    }
    xadj[numCells] = freeAdj;

    // no weight info provided
    int wgtFlag = 0;

    // C style numbering
    int numFlag = 0;

    // decomposition options. 0 = use defaults
    int options[5] = {0};  

    // processor weights. Use even weighting
    scalarField processorWeights(nProcessors_, 1.0/nProcessors_);

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
            processorWeights = 
                scalarField(metisDecompCoeffs.lookup("processorWeights"));
        }

        if (metisDecompCoeffs.found("options"))
        {
            scalarField decompOptions(metisDecompCoeffs.lookup("options"));

            if (decompOptions.size() != 5)
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 5"
                    << abort(FatalError);
            }

            Info<< "Using Metis options     " << decompOptions
                << endl << endl;

            // convert decompOptions to int
            for (int i=0; i<5; i++)
            {
                options[i] = int(decompOptions[i]);
            }
        }
    }

    processorWeights /= sum(processorWeights);

    // convert processorWeights to float
    float *weights = new float[nProcessors_];
    for (int i=0; i<nProcessors_; i++)
    {
        weights[i] = processorWeights[i];
    }


    // output: cell -> processor addressing
    labelList finalDecomp(numCells);

    // output: number of cut edges
    int edgeCut = 0;

    METIS_WPartGraphKway
    (
        &numCells,         // num vertices in graph
        xadj.begin(),      // indexing into adjncy
        adjncy.begin(),    // neighbour info
        NULL,              // no vertexweights
        NULL,              // no edgeweights
        &wgtFlag,
        &numFlag,
        &nProcessors_,
        weights,
        options,
        &edgeCut,
        finalDecomp.begin()
    );

    delete[] weights;

    return finalDecomp;
}


// ************************************************************************* //
