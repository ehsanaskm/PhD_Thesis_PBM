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

#include "PBM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(PBMmethod, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        PBMmethod,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::PBMmethod::PBMmethod
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase)//,
 //   d0_("d0", dimLength, diameterProperties_.lookup("d0")),
 //   p0_("p0", dimPressure, diameterProperties_.lookup("p0"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::PBMmethod::~PBMmethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::PBMmethod::d() const
{
    const volScalarField& d32 = phase_.U().db().lookupObject<volScalarField>
    (
        "d32"
    );


      return d32;
}  



bool Foam::diameterModels::PBMmethod::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

   // diameterProperties_.lookup("d0") >> d0_;
   // diameterProperties_.lookup("p0") >> p0_;

    return true;
}


// ************************************************************************* //
