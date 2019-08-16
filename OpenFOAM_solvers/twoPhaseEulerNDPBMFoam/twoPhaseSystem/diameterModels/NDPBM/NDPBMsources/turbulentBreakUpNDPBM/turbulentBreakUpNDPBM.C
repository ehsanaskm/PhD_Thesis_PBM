/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "turbulentBreakUpNDPBM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace NDPBMsources
{
    defineTypeNameAndDebug(turbulentBreakUpNDPBM, 0);
    addToRunTimeSelectionTable(NDPBMsource, turbulentBreakUpNDPBM, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::NDPBMsources::turbulentBreakUpNDPBM::
turbulentBreakUpNDPBM
(
    const NDPBM& ndpbm,
    const dictionary& dict
)
:
    NDPBMsource(ndpbm),
    Cti_("Cti", dimless, dict.lookup("Cti")),
    WeCr_("WeCr", dimless, dict.lookup("WeCr"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::NDPBMsources::turbulentBreakUpNDPBM::S() const
{
    tmp<volScalarField> tS_
    (
        new volScalarField
        (
            IOobject
            (
                "S_",
                ndpbm_.phase().U().time().timeName(),
                ndpbm_.phase().mesh()
            ),
            ndpbm_.phase().U().mesh(),
            dimensionedScalar("S_", dimless/dimTime, 0)
        )
    );

    volScalarField S_ = tS_();
    scalar Cti = Cti_.value();
    scalar WeCr = WeCr_.value();
    volScalarField Ut(this->Ut());
    volScalarField We(this->We());
    const volScalarField& d(ndpbm_.d()());

    forAll(S_, celli)
    {
        if (We[celli] > WeCr)
        {
            S_[celli] =
                Cti/d[celli]
               *Ut[celli]
               *sqrt(1 - WeCr/We[celli])
               *exp(-WeCr/We[celli]);
        }
    }

    return tS_;
}


