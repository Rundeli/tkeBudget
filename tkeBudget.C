
#include "tkeBudget.H"
#include "fvCFD.H"
#include "fvcGrad.H"
#include "turbulentTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(tkeBudget,0);
    addToRunTimeSelectionTable(functionObject, tkeBudget, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::tkeBudget::calculateTKEBudget()
{
    volVectorField U_ = obr().lookupObject<volVectorField>(UName_);
    volVectorField UMean_ = obr().lookupObject<volVectorField>(UMeanName_);
    volSymmTensorField UPrime2Mean_ = obr().lookupObject<volSymmTensorField>(UPrime2MeanName_);

    volScalarField p_ = obr().lookupObject<volScalarField>(pName_);
    volScalarField pMean_ = obr().lookupObject<volScalarField>(pMeanName_);

    const incompressible::turbulenceModel& turModel_ = lookupObject<incompressible::turbulenceModel>(turbulenceModel::propertiesName);
    volScalarField nu_ = turModel_.nu();

    volVectorField UPrime = U_ - UMean_;
    tmp<volTensorField> gradUPrimePtr_ = fvc::grad(UPrime);
    volTensorField& gradUPrime_ = gradUPrimePtr_.ref();

    volSymmTensorField Sij = symm(gradUPrime_);

    tmp<volTensorField> gradUMeanPtr_ = fvc::grad(UMean_);
    volTensorField& gradUMean_ = gradUMeanPtr_.ref();

    volScalarField trUP2Mean = tr(UPrime2Mean_);
    tmp<volVectorField> gradUPrime2Mean = fvc::grad(trUP2Mean);
    volVectorField& gradUPrime2Mean_ = gradUPrime2Mean.ref();

    calConvectionTerm(UMean_,gradUPrime2Mean_);
    calProductionTerm(UPrime2Mean_,gradUMean_);
    calTurbulenceDiffusionTerm(UPrime2Mean_,UPrime);
    calViscousDiffusionTerm(Sij,UPrime,nu_);
    calVPGCorelationTerm(UPrime,p_,pMean_);
    calViscousDissipationTerm(Sij,nu_);
}

void Foam::functionObjects::tkeBudget::calConvectionTerm
(
    const volVectorField& UMean,
    const volVectorField& gradUPrime2Mean
)
{
    if(convectionTerm_)
    {
        Info<< "Calculating TKE convection term" << endl;
        tmp<volScalarField> cvPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "convectionTerm",
                    UMean.mesh().time().timeName(),
                    UMean.mesh()
                ),
                -0.5*UMean & gradUPrime2Mean
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_convectionTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_convectionTerm")) = cvPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_convectionTerm",
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    cvPtr_
                )
            );
        }
    }
}

void Foam::functionObjects::tkeBudget::calProductionTerm
(
    const volSymmTensorField& UPrime2Mean,
    const volTensorField& gradUMean
)
{
    if(productionTerm_)
    {
        Info<< "Calculating TKE production term" << endl;
        tmp<volScalarField> pdPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "productionTerm",
                    UPrime2Mean.mesh().time().timeName(),
                    UPrime2Mean.mesh()
                ),
                -1.0*UPrime2Mean && gradUMean
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_productionTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_productionTerm")) = pdPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_productionTerm",
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    pdPtr_
                )
            );
        }
    }
}

void Foam::functionObjects::tkeBudget::calTurbulenceDiffusionTerm
(
    const volSymmTensorField& UPrime2Mean,
    const volVectorField& UPrime
)
{
    if(turDiffusionTerm_)
    {
        Info<< "Calculating TKE turbulence diffusion term" << endl;
        volScalarField trUP2M = tr(UPrime2Mean);
        volVectorField q2UPrime = trUP2M * UPrime;
        volScalarField divq2UPM = fvc::div(q2UPrime);

        tmp<volScalarField> tdfPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "turbulenceDiffusionTerm",
                    UPrime.mesh().time().timeName(),
                    UPrime.mesh()
                ),
                divq2UPM
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_turbulenceDiffusionTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_turbulenceDiffusionTerm")) = tdfPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_turbulenceDiffusionTerm",
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    tdfPtr_
                )
            );
        }
    }
}

void Foam::functionObjects::tkeBudget::calViscousDiffusionTerm
(
    const volSymmTensorField& Sij,
    const volVectorField& UPrime,
    const volScalarField& nu
)
{
    if(visDiffusionTerm_)
    {
        Info<< "Calculating TKE viscous diffusion term" << endl;
        volTensorField gradUP = fvc::grad(UPrime);
        volScalarField SgradUP = Sij && gradUP;
        tmp<volScalarField> vdfPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "viscousDiffusionTerm",
                    UPrime.mesh().time().timeName(),
                    UPrime.mesh()
                ),
                2.0*nu*SgradUP
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_viscousDiffusionTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_viscousDiffusionTerm")) = vdfPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_viscousDiffusionTerm",
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    vdfPtr_
                )
            );
        }
    }
}

void Foam::functionObjects::tkeBudget::calVPGCorelationTerm
(
    const volVectorField& UPrime,
    const volScalarField& p,
    const volScalarField& pMean
)
{
    if(vpgCorelationTerm_)
    {
        Info<< "Calculating TKE velocity-pressureGradient-Corelation term" << endl;
        volScalarField PPrime = p - pMean;
    
        tmp<volVectorField> gradPPrimePtr_ = fvc::grad(PPrime);
        volVectorField& gradPPrime_ = gradPPrimePtr_.ref(); 
    
        tmp<volScalarField> vpgPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "VPGCorelationTerm",
                    UPrime.mesh().time().timeName(),
                    UPrime.mesh()
                ),
                UPrime & gradPPrime_
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_VPGCorelationTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_VPGCorelationTerm")) = vpgPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_VPGCorelationTerm",
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    vpgPtr_
                )
            );
        }
    }
}

void Foam::functionObjects::tkeBudget::calViscousDissipationTerm
(
    const volSymmTensorField& Sij,
    const volScalarField& nu
)
{
    if(visDissipationTerm_)
    {
        Info<< "Calculating TKE viscous dissipation term" << endl;
        volScalarField Sij2 = Sij && Sij;
        tmp<volScalarField> vdpPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "viscousDissipationTerm",
                    Sij.mesh().time().timeName(),
                    Sij.mesh()
                ),
                2.0*nu*Sij2
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_viscousDissipationTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_viscousDissipationTerm")) = vdpPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_viscousDissipationTerm",
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    vdpPtr_
                )
            );
        }
    }
}

void Foam::functionObjects::tkeBudget::checkInsert
(
    const word &fieldName, 
    const Switch &fieldSwitch,
    HashTable<bool>& FieldsList
)
{
    if(fieldSwitch)
    {
        FieldsList.insert(fieldName,true);
    }
    else
    {
        Info<< fieldName << " not calculated!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::functionObjects::tkeBudget::tkeBudget
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    convectionTerm_(dict.get<Switch>("convectionTerm")),
    productionTerm_(dict.get<Switch>("productionTerm")),
    turDiffusionTerm_(dict.get<Switch>("turbulenceDiffusionTerm")),
    visDiffusionTerm_(dict.get<Switch>("viscousDiffusionTerm")),
    vpgCorelationTerm_(dict.get<Switch>("velocity-pressureGradient-CorelationTerm")),
    visDissipationTerm_(dict.get<Switch>("viscousDissipationTerm")),
    UName_(dict.getOrDefault<Foam::word>("U", "U")),
    UMeanName_(dict.getOrDefault<Foam::word>("UMean", "UMean")),
    UPrime2MeanName_(dict.getOrDefault<Foam::word>("UPrime2Mean", "UPrime2Mean")),
    pName_(dict.getOrDefault<Foam::word>("p","p")),
    pMeanName_(dict.getOrDefault<Foam::word>("pMean","pMean")),
    FieldsList()
{
    checkInsert("tkeBudget_convectionTerm",convectionTerm_,FieldsList);
    checkInsert("tkeBudget_productionTerm",productionTerm_,FieldsList);
    checkInsert("tkeBudget_turbulenceDiffusionTerm",turDiffusionTerm_,FieldsList);
    checkInsert("tkeBudget_viscousDiffusionTerm",visDiffusionTerm_,FieldsList);
    checkInsert("tkeBudget_VPGCorelationTerm",vpgCorelationTerm_,FieldsList);  
    checkInsert("tkeBudget_viscousDissipationTerm",visDissipationTerm_,FieldsList); 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::tkeBudget::read()
{
    Info<< "Reading relative fields" << endl;

    return 
    (
    isExist<volVectorField>(UMeanName_)                 ||
    isExist<volSymmTensorField>(UPrime2MeanName_)       ||
    isExist<volScalarField>(pMeanName_));
}

bool Foam::functionObjects::tkeBudget::execute()
{
    read();

    Info<< "Calculating TKE Budget" << endl;
    calculateTKEBudget();

    return true;
}

bool Foam::functionObjects::tkeBudget::write()
{
    Info<< "Writing TKE Budget:" << endl;

    for(auto iter = FieldsList.begin(); iter != FieldsList.end(); iter++)
    {
        if(iter.val())
        {
            Info<< "Writing " << iter.key() << " field" << endl;
            writeObject(iter.key());
        }
    }
    return true;
}

// ************************************************************************* //