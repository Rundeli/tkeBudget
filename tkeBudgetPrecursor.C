
#include "tkeBudgetPrecursor.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(tkeBudgetPrecursor,0);
    addToRunTimeSelectionTable(functionObject, tkeBudgetPrecursor, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::tkeBudgetPrecursor::UPCTermPrecursorField
(
    const volVectorField& UPrime,
    const volScalarField& pPrime
)
{
    Info<< "Calculating Velocity-Pressure-Gradient Corelation TermPrecursorField" << endl;

    volVectorField upPrime = pPrime * UPrime;

    tmp<volVectorField> uppPtr_
        (
            new volVectorField
            (
                IOobject
                (
                    "UPCTermPrecursor",
                    UPrime.mesh().time().timeName(),
                    UPrime.mesh()
                ),
                upPrime
            )
        );

    if(obr().foundObject<volVectorField>("tkeBudget_UPCTermPrecursor"))
    {
        const_cast<volVectorField&>(obr().lookupObject<volVectorField>("tkeBudget_UPCTermPrecursor")) = uppPtr_;
    }
    else
    {
        obr_.store
        (
            new volVectorField
            (
                IOobject
                (
                    "tkeBudget_UPCTermPrecursor",
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                uppPtr_
            )
        );
    }
}

void Foam::functionObjects::tkeBudgetPrecursor::TurTransTermPrecursorField
(
    const volScalarField& k,
    const volVectorField& UPrime
)
{
    Info<< "Calculating Turbulence Transport Term PrecursorField" << endl;

    volScalarField k2 = k * k;
    volVectorField k2UPrime = k2 * UPrime;

    tmp<volVectorField> ttpPtr_
        (
            new volVectorField
            (
                IOobject
                (
                    "TurTransTermPrecursor",
                    UPrime.mesh().time().timeName(),
                    UPrime.mesh()
                ),
                k2UPrime
            )
        );

    if(obr().foundObject<volVectorField>("tkeBudget_TurTransTermPrecursor"))
    {
        const_cast<volVectorField&>(obr().lookupObject<volVectorField>("tkeBudget_TurTransTermPrecursor")) = ttpPtr_;
    }
    else
    {
        obr_.store
        (
            new volVectorField
            (
                IOobject
                (
                    "tkeBudget_TurTransTermPrecursor",
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                ttpPtr_
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::functionObjects::tkeBudgetPrecursor::tkeBudgetPrecursor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_(dict.getOrDefault<Foam::word>("U", "U")),
    UMeanName_(dict.getOrDefault<Foam::word>("UMean", "UMean")),
    UPrime2MeanName_(dict.getOrDefault<Foam::word>("UPrime2Mean", "UPrime2Mean")),
    pName_(dict.getOrDefault<Foam::word>("p","p")),
    pMeanName_(dict.getOrDefault<Foam::word>("pMean","pMean"))
{
   
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::tkeBudgetPrecursor::execute()
{
    Info<< "Excute tkeBudgetPrecursor" << endl;

    volVectorField U_ = obr().lookupObject<volVectorField>(UName_);
    volVectorField UMean_ = obr().lookupObject<volVectorField>(UMeanName_);
    volSymmTensorField UPrime2Mean_ = obr().lookupObject<volSymmTensorField>(UPrime2MeanName_);

    volScalarField p_ = obr().lookupObject<volScalarField>(pName_);
    volScalarField pMean_ = obr().lookupObject<volScalarField>(pMeanName_);

    volScalarField k = 0.5*tr(UPrime2Mean_);

    volVectorField UPrime = U_ - UMean_;
    volScalarField pPrime = p_ - pMean_;

    UPCTermPrecursorField(UPrime,pPrime);
    TurTransTermPrecursorField(k,UPrime);

    return true;
}

bool Foam::functionObjects::tkeBudgetPrecursor::write()
{
    Info<< "Writing TKE Budget precursor Fields:" << endl;

    Info<< "Writing UPCTermPrecursor field" << endl;
    writeObject("tkeBudget_UPCTermPrecursor");

    Info<< "Writing Turbulence Transport Term Precursor field" << endl;
    writeObject("tkeBudget_TurTransTermPrecursor");

    return true;
}


