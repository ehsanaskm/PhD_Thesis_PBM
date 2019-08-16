/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "NDPBMsource.H"
#include "twoPhaseSystem.H"
#include "fvMatrix.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(NDPBMsource, 0);
    defineRunTimeSelectionTable(NDPBMsource, dictionary);
}
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::NDPBMsource>
Foam::diameterModels::NDPBMsource::New
(
    const word& type,
    const NDPBM& ndpbm,
    const dictionary& dict
)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "NDPBMsource::New"
            "(const word& type, const NDPBM&, const dictionary&)"
        )   << "Unknown NDPBM source type "
            << type << nl << nl
            << "Valid NDPBM source types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<NDPBMsource>(cstrIter()(ndpbm, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::Ur() const
{
    const uniformDimensionedVectorField& g =
        phase().U().db().lookupObject<uniformDimensionedVectorField>("g");

    return
        sqrt(2.0)
       *pow025
        (
            fluid().sigma()*mag(g)
           *(otherPhase().rho() - phase().rho())
           /sqr(otherPhase().rho())
        )
       *pow(max(1 - phase(), scalar(0)), 1.75);
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::Ut() const
{
    return sqrt(2*otherPhase().turbulence().k());
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::Re() const
{
    return max(Ur()*phase().d()/otherPhase().nu(), scalar(1.0e-3));
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::CD() const
{
    const volScalarField Eo(this->Eo());
    const volScalarField Re(this->Re());

    return
        max
        (
            min
            (
                (16/Re)*(1 + 0.15*pow(Re, 0.687)),
                48/Re
            ),
            8*Eo/(3*(Eo + 4))
        );
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::Mo() const
{
    const uniformDimensionedVectorField& g =
        phase().U().db().lookupObject<uniformDimensionedVectorField>("g");

    return
        mag(g)*pow4(otherPhase().nu())*sqr(otherPhase().rho())
       *(otherPhase().rho() - phase().rho())
       /pow3(fluid().sigma());
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::Eo() const
{
    const uniformDimensionedVectorField& g =
        phase().U().db().lookupObject<uniformDimensionedVectorField>("g");

    return
        mag(g)*sqr(phase().d())
       *(otherPhase().rho() - phase().rho())
       /fluid().sigma();
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::NDPBMsource::We() const
{
    return otherPhase().rho()*sqr(Ur())*phase().d()/fluid().sigma();
  // return otherPhase().rho()*sqr(Ur())*phase().d32()/fluid().sigma();
}


