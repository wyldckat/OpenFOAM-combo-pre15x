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
    Parallel cyclic interface for the use with the AMG solver.

\*---------------------------------------------------------------------------*/

#include "cyclicAmgCoupledInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicAmgCoupledInterface, 0);
addToRunTimeSelectionTable
(
    amgCoupledInterface,
    cyclicAmgCoupledInterface,
    lduInterface
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > cyclicAmgCoupledInterface::patchInternalField
(
    const Field<Type>& iF
) const
{
    const unallocLabelList& addr = addressing();

    tmp<Field<Type> > tresult(new Field<Type>(size()));
    Field<Type>& result = tresult();

    forAll (result, elemI)
    {
        result[elemI] = iF[addr[elemI]];
    }

    return tresult;
}


const labelField& cyclicAmgCoupledInterface::addressing() const
{
    if (!addrPtr_)
    {
        FatalErrorIn
        (
            "const labelField& cyclicAmgCoupledInterface::addressing() const"
        )   << "addressing not available"
            << abort(FatalError);
    }

    return *addrPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cyclicAmgCoupledInterface::cyclicAmgCoupledInterface
(
    const lduCoupledInterface* fineInterfacePtr,
    const label index
)
:
    amgCoupledInterface(fineInterfacePtr, index),
    doTransform_(false),
    forwardT_(1, tensor::zero),
    reverseT_(1, tensor::zero),
    rank_(0),
    addrPtr_(NULL)
{
    const cyclicLduCoupledInterface& p =
        refCast<const cyclicLduCoupledInterface>(*fineInterfacePtr);

    doTransform_ = p.doTransform();

    if (doTransform())
    {
        forwardT_ = p.forwardT();
        reverseT_ = p.reverseT();
    }

    rank_ = p.rank();

    // Size remains unknown until the addressing is calculated
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label cyclicAmgCoupledInterface::size() const
{
    return addressing().size();
}


// Given colouring on both sides, calculate the addressing
tmp<scalarField> cyclicAmgCoupledInterface::coeffs
(
    const labelField& localColour,
    const labelField&,
    const scalarField& fineCoeffs
) const
{
    // Make a lookup table of entries for owner/neighbour
    HashTable<SLList<label>, label, Hash<label> > neighbourLookup
        (localColour.size());

    HashTable<SLList<scalar>, label, Hash<label> > neiCoeffs
        (localColour.size());

    label nCoarseFaces = 0;

    label sizeBy2 = localColour.size()/2;

    for (label fineFaceI = 0; fineFaceI < sizeBy2; fineFaceI++)
    {
        label curMaster = localColour[fineFaceI];
        label curSlave = localColour[fineFaceI + sizeBy2];

        // Look for the master cell.  If it has already got a face,
        // add the coefficient to the face.  If not, create a new
        // face.
        if (neighbourLookup.found(curMaster))
        {
            // Check all current neighbours to see if the current
            // slave already exists.  If so, add the coefficient.

            SLList<label>& curNbrs = neighbourLookup.find(curMaster)();

            SLList<scalar>& curNeiCoeffs = neiCoeffs.find(curMaster)();

            bool nbrFound = false;

            // Warning. Synchronous iterators
            SLList<label>::iterator curNbrsIter = curNbrs.begin();

            SLList<scalar>::iterator curNeiCoeffsIter =
                curNeiCoeffs.begin();

            for
            (
                ;

                curNbrsIter != curNbrs.end(),
                curNeiCoeffsIter != curNeiCoeffs.end();

                ++curNbrsIter,
                ++curNeiCoeffsIter
            )
            {
                if (curNbrsIter() == curSlave)
                {
                    nbrFound = true;

                    curNeiCoeffsIter() += fineCoeffs[fineFaceI];

                    break;
                }
            }

            if (!nbrFound)
            {
                curNbrs.append(curSlave);

                curNeiCoeffs.append(fineCoeffs[fineFaceI]);

                // New coarse face created
                nCoarseFaces++;
            }
        }
        else
        {
            // This master has got no neighbours yet.  Add a neighbour
            // and a coefficient, thus creating a new face
            neighbourLookup.insert(curMaster, SLList<label>(curSlave));
            neiCoeffs.insert(curMaster, SLList<scalar>(fineCoeffs[fineFaceI]));

            // New coarse face created
            nCoarseFaces++;
        }
    } // end for all fine faces

    // Create the coarse patch and coarse patch coefficients
    if (addrPtr_)
    {
        FatalErrorIn
       (
           "tmp<scalarField> cyclicAmgCoupledInterface::coeffsAddressing\n"
           "(\n"
           "    const labelField& localColour,\n"
           "    const labelField& nbrColour,\n"
           "    const scalarField& fineCoeffs\n"
           ") const"
       )   << "coefficients and addressing already calculated."
           << abort(FatalError);
    }

    // Allocate addressing and coefficients

    addrPtr_ = new labelField(2*nCoarseFaces, -1);
    labelField& coarseAddr = *addrPtr_;

    tmp<scalarField> tcoarseCoeffs(new scalarField(2*nCoarseFaces));
    scalarField& coarseCoeffs = tcoarseCoeffs();

    labelList contents = neighbourLookup.toc();

    // Reset face counter for re-use
    nCoarseFaces = 0;

    // On master side, the owner addressing is stored in table of contents
    forAll (contents, masterI)
    {
        SLList<label>& curNbrs = neighbourLookup.find(contents[masterI])();
        SLList<scalar>& curNeiCoeffs = neiCoeffs.find(contents[masterI])();

        // Warning. Synchronous iterators
        SLList<label>::iterator curNbrsIter = curNbrs.begin();

        SLList<scalar>::iterator curNeiCoeffsIter = curNeiCoeffs.begin();

        for
        (
            ;

            curNbrsIter != curNbrs.end(),
            curNeiCoeffsIter != curNeiCoeffs.end();

            ++curNbrsIter,
            ++curNeiCoeffsIter
        )
        {
            coarseAddr[nCoarseFaces] = contents[masterI];
            coarseCoeffs[nCoarseFaces] = curNeiCoeffsIter();
            nCoarseFaces++;
        }
    }

    // On slave side, the owner addressing is stored in linked lists
    forAll (contents, masterI)
    {
        SLList<label>& curNbrs = neighbourLookup.find(contents[masterI])();
        SLList<scalar>& curNeiCoeffs = neiCoeffs.find(contents[masterI])();

        // Warning. Synchronous iterators
        SLList<label>::iterator curNbrsIter = curNbrs.begin();

        SLList<scalar>::iterator curNeiCoeffsIter = curNeiCoeffs.begin();

        for
        (
            ;

            curNbrsIter != curNbrs.end(),
            curNeiCoeffsIter != curNeiCoeffs.end();

            ++curNbrsIter,
            ++curNeiCoeffsIter
        )
        {
            coarseAddr[nCoarseFaces] = curNbrsIter();
            coarseCoeffs[nCoarseFaces] = curNeiCoeffsIter();
            nCoarseFaces++;
        }
    }

    return tcoarseCoeffs;
}


// Return neighbour colouring
tmp<labelField> cyclicAmgCoupledInterface::nbrColour
(
    const labelField& cField
) const
{
    const unallocLabelList& addr = addressing();

    tmp<labelField> tpnf(new labelField(size()));
    labelField& pnf = tpnf();

    label sizeby2 = size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = cField[addr[facei + sizeby2]];
        pnf[facei + sizeby2] = cField[addr[facei]];
    }

    return tpnf;
}


// Transfer internal cell data to neighbour cyclic
// and receive and return neighbour cyclic internal cell data
void cyclicAmgCoupledInterface::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction
) const
{
    const unallocLabelList& addr = addressing();

    scalarField pnf(size());

    label sizeby2 = size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = psiInternal[addr[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[addr[facei]];
    }

    // Multiply the field by coefficients and add into the result

    forAll(addr, elemI)
    {
        result[addr[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
