bool closedVolume = false;

pd.boundaryField() ==
    p.boundaryField() - rho.boundaryField()*gh.boundaryField() - pRef.value();

rho = thermo->rho();

volScalarField A = UEqn.A();
U = UEqn.H()/A;

phi =
(
    fvc::interpolate
    (
        rho*(U + UphiCoeff*fvc::ddt0(rho, U)/A),
        "interpolate((H(U)|A(U)))"
    ) & mesh.Sf()
) - UphiCoeff*fvc::interpolate(rho/A)*fvc::ddt0(phi)
  - fvc::interpolate(rho*gh/A)*fvc::snGrad(rho)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pdEqn
    (
        fvm::ddt(psi, pd)
      + fvc::ddt(psi)*pRef
      + fvc::ddt(psi, rho)*gh
      + fvc::div(phi)
      - fvm::laplacian(rho/A, pd)
    );

    closedVolume = pd.needReference();
    pdEqn.solve();

    if (nonOrth == nNonOrthCorr)
    {
        phi += pdEqn.flux();
    }
}

p = pd + rho*gh + pRef;
dpdt = fvc::ddt(p);

#include "rhoEqn.H"
#include "continuityErrs.H"

U -= (fvc::grad(pd) + fvc::grad(rho)*gh)/A;
U.correctBoundaryConditions();


// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
        /fvc::domainIntegrate(thermo->psi());
    rho = thermo->rho();
}
