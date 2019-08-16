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
    Foam::diameterModels::NDPBMsource

Description
    NDPBM (Number Density Approach Population Balance Model) bubble diameter model
    run-time selectable sources.

SourceFiles
    NDPBMsource.C

\*---------------------------------------------------------------------------*/

#ifndef NDPBMsource_H
#define NDPBMsource_H

#include "NDPBM.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class NDPBMsource Declaration
\*---------------------------------------------------------------------------*/

class NDPBMsource
{

protected:

    // Protected data

        //- Reference to the NDPBM this source applies to
        const NDPBM& ndpbm_;


public:

    //- Runtime type information
    TypeName("NDPBMsource");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            NDPBMsource,
            dictionary,
            (
                const NDPBM& ndpbm,
                const dictionary& dict
            ),
            (ndpbm, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of NDPBM sources
    class iNew
    {
        const NDPBM& ndpbm_;

    public:

        iNew(const NDPBM& ndpbm)
        :
            ndpbm_(ndpbm)
        {}

        autoPtr<NDPBMsource> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return NDPBMsource::New(type, ndpbm_, dict);
        }
    };


    // Constructors

        NDPBMsource(const NDPBM& ndpbm)
        :
            ndpbm_(ndpbm)
        {}

        autoPtr<NDPBMsource> clone() const
        {
            notImplemented("autoPtr<NDPBMsource> clone() const");
            return autoPtr<NDPBMsource>(NULL);
        }


    // Selectors

        static autoPtr<NDPBMsource> New
        (
            const word& type,
            const NDPBM& ndpbm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~NDPBMsource()
    {}


    // Member Functions

        const phaseModel& phase() const
        {
            return ndpbm_.phase();
        }

        const twoPhaseSystem& fluid() const
        {
            return ndpbm_.phase().fluid();
        }

        const phaseModel& otherPhase() const
        {
            return phase().otherPhase();
        }

        scalar phi() const
        {
            return 1.0/(36*constant::mathematical::pi);
        }

        //- Return the bubble relative velocity
        tmp<volScalarField> Ur() const;



        //- Return the bubble turbulent velocity
        tmp<volScalarField> Ut() const;

        //- Return the bubble Reynolds number
        tmp<volScalarField> Re() const;

        //- Return the bubble drag coefficient
        tmp<volScalarField> CD() const;

        //- Return the bubble Morton number
        tmp<volScalarField> Mo() const;

        //- Return the bubble Eotvos number
        tmp<volScalarField> Eo() const;

        //- Return the bubble Webber number
        tmp<volScalarField> We() const;

        //


        virtual tmp<volScalarField> S() const = 0;
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
