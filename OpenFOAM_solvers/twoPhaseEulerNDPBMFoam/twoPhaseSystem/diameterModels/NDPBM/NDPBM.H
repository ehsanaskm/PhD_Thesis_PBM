/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::diameterModels::NDPBM

Description
    NDPBM (Number Density Approach Population Balance Model) bubble diameter model.

    Solves for the interfacial curvature per unit volume of the phase rather
    than interfacial area per unit volume to avoid stability issues relating to
    the consistency requirements between the phase fraction and interfacial area
    per unit volume.  In every other respect this model is as presented in the
    paper:

    \verbatim
        "CFD modeling of gas dispersion and bubble size in a double turbine stirred tank"
        F. Kerdouss,
        A. Bannari,
        P. Proulx,
        Chemical Engineering Science, Vol.61 December 2006
    \endverbatim

SourceFiles
    NDPBM.C

\*---------------------------------------------------------------------------*/

#ifndef NDPBM_H
#define NDPBM_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

// Forward declaration of classes
class NDPBMsource;

/*---------------------------------------------------------------------------*\
                           Class NDPBM Declaration
\*---------------------------------------------------------------------------*/

class NDPBM
:
    public diameterModel
{
    // Private data

        //- Interfacial curvature (alpha*interfacial area)
 //       volScalarField kappai_;
        //add 
        volScalarField n_;
        //add
        volScalarField S_;

      
        //- Maximum diameter used for stabilisation in the limit kappai->0
        dimensionedScalar dMax_;

        //- Minimum diameter used for stabilisation in the limit kappai->inf
        dimensionedScalar dMin_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- The Sauter-mean diameter of the phase
        volScalarField d_;

 //       volScalarField d32_;   //add

        //- NDPBM sources
        PtrList<NDPBMsource> sources_;


    // Private member functions

        tmp<volScalarField> dsm() const;

//        tmp<volScalarField> d32sm() const;   //add
public:

    friend class NDPBMsource;

    //- Runtime type information
    TypeName("NDPBM");


    // Constructors

        //- Construct from components
        NDPBM
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~NDPBM();


    // Member Functions

        //- Return the interfacial curvature
 /*       const volScalarField& kappai() const
        {
            return kappai_;
        }
*/
        const volScalarField& n() const   //add
        {
            return n_;
        }

/*        //- Return the interfacial area
        tmp<volScalarField> a() const
        {
            return phase_*kappai_;
        }
*/
        //- Return the Sauter-mean diameter
        virtual tmp<volScalarField> d() const
        {
            return d_;
        }
/*
        virtual tmp<volScalarField> d32() const  //add
        {
            return d32_;
        }
*/
        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
