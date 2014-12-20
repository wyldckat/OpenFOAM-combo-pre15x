if (algorithm == "flux")
{

#   include "algorithms/flux/pEqn.H"

}
else if (algorithm == "Up")
{

#   include "algorithms/Up/pEqn.H"

}
else
{
    FatalErrorIn(args.executable())
        << "unknown algorithm " << algorithm
        << " expected 'flux' of 'Up'"
        << exit(FatalError);
}
