/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::tkeBudget

Group
    grpFieldFunctionObjects

Description
    Computes turbulence kinetic energy budget term in incompressible flows.

    Operands:
    \table
      Operand      | Type                          | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
               --> |$FOAM_CASE/\<time\>/\<inpField\>s
      output file  | -                             | -
      output field | {vol,surface}\<Type\>Field(s) <!--
               --> | $FOAM_CASE/\<time\>/\<outField\>s
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    libs("libtkeBudget.so");

    tkeBudget
    {
        type tkeBudget;
        libs (tkeBudget);
        writeInterval   1;
        timeStart       1;
        timeEnd         1;
        rho             1000;

        convectionTerm                              on;
        productionTerm                              on;
        turbulenceTransportTerm                     on;
        viscousTransportTerm                        on;
        velocity-pressureGradient-CorelationTerm    on;
        viscousDissipation                          on;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: tkeBudget               | word |  yes  | -
      libs         | Library name: tkeBudget            | word |  yes  | -
      rho          | Fluid Density                      | scalar| yes  | -
      fields       | Names of the operand fields and averaging options <!--
               --> | dict |  yes  | -
      restartOnRestart| Restart the averaging on restart | bool | no     | false
      restartOnOutput | Restart the averaging on output  | bool | no     | false
      periodicRestart | Periodically restart the averaging | bool | no   | false
      restartPeriod   | Periodic restart period   | scalar | conditional | -
      restartTime     | One-shot reset of the averaging  | scalar | no   | GREAT
      subRegion    | Name for alternative objectRegistry | word   | no   | ""
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link tkeBudget.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldAverageItem

SourceFiles
    tkeBudget.C
  
\*---------------------------------------------------------------------------*/

#ifndef functionObjects_tkeBuget_H
#define functionObjects_tkeBuget_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{   

/*---------------------------------------------------------------------------*\
                    Class tkeBudget Declaration
\*---------------------------------------------------------------------------*/

class tkeBudget
:
    public fvMeshFunctionObject
{

protected:

    // Protected Data

        //- Fluid Density
        scalar rho_;

        //- Calculate the tkeBuget convention term
        Switch convectionTerm_;

        //- Calculate the tkeBuget production term
        Switch productionTerm_;

        //- Calculate the tkeBuget turbulence transport term
        Switch turTransportTerm_;

        //- Calculate the tkeBuget viscousity transport term
        Switch visTransportTerm_;

        //- Calculate the tkeBuget velocity-pressure_gradient correlation term
        Switch vpgCorelationTerm_;

        //- Calculate the tkeBuget viscousity dissipation term
        Switch visDissipationTerm_;

        //- Name of velocity field
        word UName_;

        //- Name of mean velocity field
        word UMeanName_;

        //- Name of prime2mean velocity field
        word UPrime2MeanName_;

        //- Name of pressure field
        word pName_;

        //- Name of mean pressure field
        word pMeanName_;

        /*
        //- Velocity field
        const volVectorField U_;

        //- Mean velocity field
        const volVectorField UMean_;

        //- Prime2Mean velocity field
        const volSymmTensorField UPrime2Mean_;

        //- Kinematic viscosity field
        const volScalarField nu_;

        //- Pressure field
        const volScalarField p_;

        //- Mean pressure field
        const volScalarField pMean_;
        */

/*
        //- Temporary field for conventionTerm field
        tmp<volScalarField> convectionTermField_;

        //- Temporary field for productionTerm field
        tmp<volScalarField> productionTermField_;

        //- Temporary field for turbulenceDiffusionTerm field
        tmp<volScalarField> turbulenceDiffusionTermField_;

        //- Temporary field for viscousityDiffusionTerm field
        tmp<volScalarField> viscousDiffusionTermField_;

        //- Temporary field for velocity-pressure_gradient correlation term field
        tmp<volScalarField> vpgCorelationTermField_;
*/
        //- Fields name-switch mapping
        HashTable<bool> FieldsList;

    //protected member functions

        //- Calculate the tkeBuget Fields
        void calculateTKEBudget();

        //- Calculate the tkeBuget ConvectionTerm Fields
        void calConvectionTerm
        (
            const volVectorField& UMean,
            const volScalarField& k
        );

        //- Calculate the tkeBuget ProductionTerm Fields
        void calProductionTerm
        (
            const volSymmTensorField& UPrime2Mean,
            const volTensorField& gradUMean
        );

        //- Calculate the tkeBuget TurbulenceTransportTerm Fields
        void calTurbulenceTransportTerm
        (
            const volScalarField& k,
            const volVectorField& UPrime
        );

        //- Calculate the tkeBuget ViscousTransportTerm Fields
        void calViscousTransportTerm
        (
            const volSymmTensorField& Sij,
            const volVectorField& UPrime,
            const volScalarField& nu
        );

        //- Calculate the tkeBuget Velocity-PressureGradient CorelationTerm Fields
        void calVPGCorelationTerm
        (
            const volVectorField& UPrime,
            const volScalarField& p,
            const volScalarField& pMean
        );

        //- Calculate the tkeBuget ViscousDissipationTerm Fields
        void calViscousDissipationTerm
        (
            const volSymmTensorField& Sij,
            const volScalarField& nu
        );

        //- Insert the defined field name-switch Mapping
        void checkInsert
        (
            const word& fieldName,
            const Switch& fieldSwitch,
            HashTable<bool>& FieldsList
        );

        template<class ftp>
        bool isExist(const word& fieldName)
        {
            if(!(obr().foundObject<ftp>(fieldName)))
            {
                Info<< "Field " << fieldName << " not found" <<nl
                    << "Calculate " << fieldName << " field first!" << endl;
                
                    return false;
            }
            return true;
        }

public:

    //- Runtime type information
    TypeName("tkeBudget");

    // Constructors

        //- Construct from Time and dictionary
        tkeBudget
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Destructor
        virtual ~tkeBudget() = default;

    // Member Functions

        //- Read the input field data
        virtual bool read();

        //- Calculate the tkeBudget Fields
        virtual bool execute();

        //- Write the tkeBudget Fields
        virtual bool write();

};

} // End namespace functionObjects

} // End namespace Foam

#endif
