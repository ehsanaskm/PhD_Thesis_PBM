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

Class
    Foam::diameterModels::NDPBMsources::randomCoalescenceNDPBM

Description
    Random coalescence NDPBM source as defined in paper:

    \verbatim
        "CFD modeling of gas dispersion and bubble size in a double turbine stirred tank"
        F. Kerdouss,
        A. Bannari,
        P. Proulx,
        Chemical Engineering Science, Vol.61 December 2006
    \endverbatim



SourceFiles
    randomCoalescenceNDPBM.C

\*---------------------------------------------------------------------------*/

#ifndef randomCoalescenceNDPBM_H
#define randomCoalescenceNDPBM_H

#include "NDPBMsource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace NDPBMsources
{

/*---------------------------------------------------------------------------*\
                           Class randomCoalescenceNDPBM Declaration
\*---------------------------------------------------------------------------*/

class randomCoalescenceNDPBM
:
    public NDPBMsource
{
    // Private data

        dimensionedScalar Crc_;
        dimensionedScalar C_;
        dimensionedScalar alphaMax_;


public:

    //- Runtime type information
    TypeName("randomCoalescenceNDPBM");

    // Constructors

        randomCoalescenceNDPBM
        (
            const NDPBM& ndpbm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~randomCoalescenceNDPBM()
    {}


    // Member Functions

        virtual tmp<volScalarField> S() const;
       
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace NDPBMsources
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
