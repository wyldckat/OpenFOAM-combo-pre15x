pd.boundaryField() ==
    p.boundaryField() - rho.boundaryField()*gh.boundaryField() - pRef.value();

volScalarField AU = UEqn().A();
U = UEqn().H()/AU;
UEqn.clear();
phi = (fvc::interpolate(rho*U) & mesh.Sf())
    - fvc::interpolate(rho*gh/AU)*fvc::snGrad(rho)*mesh.magSf();
bool closedVolume = adjustPhi(phi, U, p);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pdEqn
    (
        fvm::laplacian(rho/AU, pd) == fvc::div(phi)
    );

    fvScalarMatrix::reference pRef = pdEqn.setReference(pRefCell, pRefValue);
    pdEqn.solve();
    pdEqn.unsetReference(pRef);

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pdEqn.flux();
    }
}

#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
pd.relax();

p = pd + rho*gh + pRef;

U -= (fvc::grad(pd) + fvc::grad(rho)*gh)/AU;
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
        /fvc::domainIntegrate(thermo->psi());
}

rho = thermo->rho();
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
