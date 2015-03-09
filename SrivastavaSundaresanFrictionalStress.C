/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SrivastavaSundaresanFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(SrivastavaSundaresan, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SrivastavaSundaresan,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::SrivastavaSundaresan::
SrivastavaSundaresan
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_.lookup("Fr")),
    eta_("eta", dimless, coeffDict_.lookup("eta")),
    p_("p", dimless, coeffDict_.lookup("p")),
    phi_("phi", dimless, coeffDict_.lookup("phi"))
{
    phi_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::SrivastavaSundaresan::
~SrivastavaSundaresan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::SrivastavaSundaresan::
frictionalPressure
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{

    return
        Fr_*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta_)
       /pow(max(alphaMax - alpha1, scalar(5.0e-2)), p_);
} //critical state pressure using Johnson-Jackson frictional stress model


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::SrivastavaSundaresan::
frictionalPressurePrime
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax
) const
{
    return Fr_*
    (
        eta_*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta_ - 1.0)
       *(alphaMax-alpha1)
      + p_*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta_)
    )/pow(max(alphaMax - alpha1, scalar(5.0e-2)), p_ + 1.0);
} //derivative of critical state pressure using Johnson-Jackson frictional stress model


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::SrivastavaSundaresan::nu
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    const scalar I2Dsmall = 1.0e-15;

    // Creating nu assuming it should be 0 on the boundary which may not be
    // true
    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "SrivastavaSundaresan:nu",
                alpha1.mesh().time().timeName(),
                alpha1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            alpha1.mesh(),
            dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0.0)
        )
    );

    volScalarField& nuf = tnu();

    forAll (D, celli)
    {
        if (alpha1[celli] > alphaMax.value() - 5e-2)
        {
            nuf[celli] =
                0.5*pf[celli]*sin(phi_.value())
               /(
                    sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                  + sqr(D[celli].yy() - D[celli].zz())
                  + sqr(D[celli].zz() - D[celli].xx()))
                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                  + sqr(D[celli].yz())) + I2Dsmall
                );
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


bool Foam::kineticTheoryModels::frictionalStressModels::SrivastavaSundaresan::read()
{
    coeffDict_ <<= dict_.subDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;

    return true;
}


// ************************************************************************* //
