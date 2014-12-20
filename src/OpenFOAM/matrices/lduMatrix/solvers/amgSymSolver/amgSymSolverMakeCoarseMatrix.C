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
    Given a level index calculate the next coarser matrix.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "amgSymSolver.H"
#include "amgCoupledInterface.H"
#include "SLList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void amgSymSolver::makeCoarseMatrix(const label fineLevelIndex)
{
    // Construct the coarse matrix and ldu addressing for the next level
    // Algorithm:
    // 1) Loop through all fine faces. If the cluster labels on two sides are
    //    different, this creates a coarse face. Define owner and neighbour
    //    for this face based on cluster IDs.
    // 2) Check if the face has been seen before. If yes, add the coefficient
    //    to the appropriate field (stored with the cell). If no, create a new
    //    face with neighbour ID and add the coefficient
    // 3) Once all the faces have been created, loop through all clusters and
    //    insert the faces in the upper order. At the same time, collect the
    //    owner and neighbour addressing.
    // 4) Agglomerate the diagonal by summing up the fine diagonal

    // Get addressing
    const lduMatrix& fineMatrix = matrixLevel(fineLevelIndex);

    const unallocLabelList& upperAddr = fineMatrix.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = fineMatrix.lduAddr().lowerAddr();

    // Get off-diagonal matrix coefficients
    const scalarField& upper = fineMatrix.upper();

    // Get restriction map for current level
    const labelField& restrictMap = restrictAddressing_[fineLevelIndex];

#   ifdef FULLDEBUG
    if (restrictMap.size() != fineMatrix.lduAddr().size())
    {
        FatalErrorIn
        (
            "amgSymSolver::makeCoarseMatrix(const label fineLevelIndex)"
        )   << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << restrictMap.size()
            << " nEqns: " << fineMatrix.lduAddr().size()
            << abort(FatalError);
    }
#   endif


    // Count the number of coarse cells
    const label nCoarseEqns = max(restrictMap) + 1;

    // Storage for block neighbours and coefficients
    List<SLList<label> > blockNbrs(nCoarseEqns);
    List<SLList<scalar> > blockFaceCoeffs(nCoarseEqns);

    // Counter for coarse faces
    label nCoarseFaces = 0;

    // Calculate the coarse diagonal. For each fine face inside of the
    // coarse cell, add twice the fine face coefficient into the diagonal
    scalarField coarseDiag(nCoarseEqns, 0.0);

    // Loop through all fine faces

    forAll (upperAddr, fineFaceI)
    {
        if
        (
            restrictMap[upperAddr[fineFaceI]]
         == restrictMap[lowerAddr[fineFaceI]]
        )
        {
            // face internal to cell. Add twice the face coefficient
            // into the diagonal
            coarseDiag[restrictMap[upperAddr[fineFaceI]]] +=
                2.0*upper[fineFaceI];
        }
        else
        {
            // this face is a part of a coarse face

            // get coarse owner and neighbour
            const label coarseOwner =
                min
                (
                    restrictMap[upperAddr[fineFaceI]],
                    restrictMap[lowerAddr[fineFaceI]]
                );

            const label coarseNeighbour =
                max
                (
                    restrictMap[upperAddr[fineFaceI]],
                    restrictMap[lowerAddr[fineFaceI]]
                );

            // check the owner block to see if this face has already been found
            SLList<label>& curBlockNbrs = blockNbrs[coarseOwner];

            SLList<scalar>& curBlockFaceCoeffs = blockFaceCoeffs[coarseOwner];

            bool nbrFound = false;

            // Warning. Synchronous iterators
            SLList<label>::iterator curBlockNbrsIter = curBlockNbrs.begin();

            SLList<scalar>::iterator curBlockFaceCoeffsIter =
                curBlockFaceCoeffs.begin();

            for
            (
                ;

                curBlockNbrsIter != curBlockNbrs.end(),
                curBlockFaceCoeffsIter != curBlockFaceCoeffs.end();

                ++curBlockNbrsIter,
                ++curBlockFaceCoeffsIter
            )
            {
                if (curBlockNbrsIter() == coarseNeighbour)
                {
                    nbrFound = true;

                    curBlockFaceCoeffsIter() += upper[fineFaceI];

                    break;
                }
            }

            if (!nbrFound)
            {
                curBlockNbrs.append(coarseNeighbour);

                curBlockFaceCoeffs.append(upper[fineFaceI]);

                // new coarse face created
                nCoarseFaces++;
            }
        }
    } // end for all fine faces

    // All coarse faces created. Make the coarse off-diagonal matrix
    labelList coarseOwner(nCoarseFaces, -1);
    labelList coarseNeighbour(nCoarseFaces, -1);

    scalarField coarseUpper(nCoarseFaces, 0.0);

    // Reorganise the storage of the coefficients into owner-neighbour
    // addressing.

    label coarseFaceI = 0;

    forAll (blockNbrs, blockI)
    {
        SLList<label>::iterator curBlockNbrsIter = blockNbrs[blockI].begin();

        SLList<scalar>::iterator curBlockFaceCoeffsIter =
            blockFaceCoeffs[blockI].begin();

        for
        (
            ;

            curBlockNbrsIter != blockNbrs[blockI].end(),
            curBlockFaceCoeffsIter != blockFaceCoeffs[blockI].end();

            ++curBlockNbrsIter,
            ++curBlockFaceCoeffsIter
        )
        {
            coarseOwner[coarseFaceI] = blockI;
            coarseNeighbour[coarseFaceI] = curBlockNbrsIter();

            coarseUpper[coarseFaceI] = curBlockFaceCoeffsIter();

            coarseFaceI++;
        }
    }

    // Create coarse-level coupled interfaces

    // Get reference to fine-level interfaces
    const lduCoupledInterfacePtrsList& fineInterfaces =
        interfaceLevel(fineLevelIndex);

    // Get reference to fine-level coefficients
    const FieldField<Field, scalar>& fineInterfaceCoeffs =
        interfaceCoeffsLevel(fineLevelIndex);

    // Create coarse interfaces, addressing and coefficients
    // Hook the coarse interfaces and coefficients
    interfaceLevels_.hook
    (
        new lduCoupledInterfacePtrsList(fineInterfaces.size())
    );

    interfaceCoeffs_.hook
    (
        new FieldField<Field, scalar>(fineInterfaces.size())
    );

    labelListList coarseInterfaceAddr(fineInterfaces.size());

    // Note: references offset by one since the fine level is stored separately
    // 
    lduCoupledInterfacePtrsList&  coarseInterfaces =
        interfaceLevels_[fineLevelIndex];
    FieldField<Field, scalar>& coarseInterfaceCoeffs =
        interfaceCoeffs_[fineLevelIndex];

    // Initialise transfer of colouring on the interface
    forAll (fineInterfaces, intI)
    {
        fineInterfaces[intI]->initNbrColour(restrictMap, true);
    }

    // Add the coarse level
    forAll (fineInterfaces, intI)
    {
        // Coarse level interface
        coarseInterfaces[intI] =
            amgCoupledInterface::New(fineInterfaces[intI], intI).ptr();

        // Get the local colour
        const unallocLabelList& localAddr =
            fineMatrix.lduAddr().patchAddr(intI);

        labelField localColour(localAddr.size());

        forAll (localColour, faceI)
        {
            localColour[faceI] = restrictMap[localAddr[faceI]];
        }

        const amgCoupledInterface& amgInterface =
            refCast<const amgCoupledInterface>(*coarseInterfaces[intI]);

        // Coefficients and addressing
        coarseInterfaceCoeffs.hook
        (
            new scalarField
            (
                amgInterface.coeffs
                (
                    localColour,
                    fineInterfaces[intI]->nbrColour(restrictMap),
                    fineInterfaceCoeffs[intI]
                )()
            )
        );

        coarseInterfaceAddr[intI] = amgInterface.addressing();
    }
        
    // Matrix restriction done!

    // Hook the coarse ldu addressing onto the list
    addrLevels_.hook
    (
        new lduAddressingStore
        (
            nCoarseEqns,
            coarseOwner,
            coarseNeighbour,
            coarseInterfaceAddr
        )
    );

    // Hook the coarse level matrix
    matrixLevels_.hook
    (
        new lduMatrix
        (
            addrLevels_[fineLevelIndex],
            fineMatrix.lduCoupledInterfaceSchedule()
        )
    );

    // Insert coarse upper
    lduMatrix& coarseMatrix = matrixLevels_[fineLevelIndex];

    coarseMatrix.upper() = coarseUpper;

    coarseMatrix.diag() =
        coarseDiag
      + restrictField(fineMatrix.diag(), fineLevelIndex);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
