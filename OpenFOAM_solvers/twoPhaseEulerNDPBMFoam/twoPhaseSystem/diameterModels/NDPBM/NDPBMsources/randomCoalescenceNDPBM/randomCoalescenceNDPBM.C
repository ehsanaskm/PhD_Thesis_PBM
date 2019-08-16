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

#include "randomCoalescenceNDPBM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace NDPBMsources
{
    defineTypeNameAndDebug(randomCoalescenceNDPBM, 0);
    addToRunTimeSelectionTable(NDPBMsource, randomCoalescenceNDPBM, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::NDPBMsources::randomCoalescenceNDPBM::
randomCoalescenceNDPBM
(
    const NDPBM& ndpbm,
    const dictionary& dict
)
:
    NDPBMsource(ndpbm),
    Crc_("Crc", dimless, dict.lookup("Crc")),
    C_("C", dimless, dict.lookup("C")),
    alphaMax_("alphaMax", dimless, dict.lookup("alphaMax"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::NDPBMsources::randomCoalescenceNDPBM::S() const
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

    scalar Crc = Crc_.value();
    scalar C = C_.value();
    scalar alphaMax = alphaMax_.value();
    volScalarField Ut(this->Ut());
    const volScalarField& alpha = phase();
    const volScalarField& n = ndpbm_.n();
    scalar cbrtAlphaMax = cbrt(alphaMax);
    const volScalarField& d(ndpbm_.d()());

    forAll(S_, celli)
    {

           

            S_[celli] =
                -1*n[celli]*d[celli]*d[celli]
               *0.05
               *Ut[celli]
               /(1 - cbrt(alpha[celli])+1.0e-15);
               
        
    }

    return tS_;
}


// ************************************************************************* //
