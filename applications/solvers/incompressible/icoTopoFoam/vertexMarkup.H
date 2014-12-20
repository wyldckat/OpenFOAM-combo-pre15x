#ifndef VERTEXMARKUP
#define VERTEXMARKUP

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    tmp<scalarField> vertexMarkup
    (
        const pointField& p,
        const scalar& curLeft,
        const scalar& curRight
    );
}

#endif
