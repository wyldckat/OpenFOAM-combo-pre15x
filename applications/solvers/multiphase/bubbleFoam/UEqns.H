fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);

if (algorithm == "flux")
{

#   include "algorithms/flux/UEqns.H"

}
else if (algorithm == "Up")
{

#   include "algorithms/Up/UEqns.H"

}
else
{
    FatalErrorIn(args.executable())
        << "unknown algorithm " << algorithm
        << " expected 'flux' of 'Up'"
        << exit(FatalError);
}
