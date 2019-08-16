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

Class
    Foam::massTransferModels::HigbiePenetration

Description
    HigbiePenetration correlation for turbulent mass transfer from the surface of a
    sphere to the surrounding fluid.

SourceFiles
    HigbiePenetration.C

\*---------------------------------------------------------------------------*/

#ifndef HigbiePenetration_H
#define HigbiePenetration_H

#include "massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class HigbiePenetration Declaration
\*---------------------------------------------------------------------------*/

class HigbiePenetration
:
    public massTransferModel
{
private:

    // Private data

        //- Lewis number
        const dimensionedScalar Le_;


public:

    //- Runtime type information
    TypeName("HigbiePenetration");


    // Constructors

        //- Construct from components
        HigbiePenetration
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~HigbiePenetration();


    // Member Functions

        //- The implicit mass transfer coefficient
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
