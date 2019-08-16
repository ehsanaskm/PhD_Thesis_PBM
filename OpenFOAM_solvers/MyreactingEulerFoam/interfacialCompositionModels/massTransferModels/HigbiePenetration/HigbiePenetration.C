/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "HigbiePenetration.H"
#include "phasePair.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(HigbiePenetration, 0);
    addToRunTimeSelectionTable(massTransferModel, HigbiePenetration, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::HigbiePenetration::HigbiePenetration
(
    const dictionary& dict,
    const phasePair& pair
)
:
    massTransferModel(dict, pair),
    Le_("Le", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::HigbiePenetration::~HigbiePenetration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::HigbiePenetration::K() const
{


    const fvMesh& mesh = pair_.dispersed().mesh();
dimensionedScalar diffCoef ("diffCoef",dimensionSet(0,2,-1,0,0,0,0),2.01*1e-9); 

 volVectorField Uslip =
   mesh.objectRegistry::lookupObject<volVectorField>("Uslip");

  volScalarField KLa=1.12866*(sqrt(diffCoef*mag(Uslip)/pair_.dispersed().d()))*6*pair_.dispersed()/pair_.dispersed().d();



      return KLa/diffCoef;
}


// ************************************************************************* //
