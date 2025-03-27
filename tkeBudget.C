
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

    volVectorField vpgPrecursorMean_ = obr().lookupObject<volVectorField>("tkeBudget_UPCTermPrecursor");
    volVectorField turTransPrecursorMean_ = obr().lookupObject<volVectorField>("tkeBudget_TurTransTermPrecursor");

    const incompressible::turbulenceModel& turModel_ = lookupObject<incompressible::turbulenceModel>(turbulenceModel::propertiesName);
    volScalarField nu_ = turModel_.nu();

    volVectorField UPrime = U_ - UMean_;
    tmp<volTensorField> gradUPrimePtr_ = fvc::grad(UPrime);
    volTensorField& gradUPrime_ = gradUPrimePtr_.ref();

    volSymmTensorField Sij = symm(gradUPrime_);

    tmp<volTensorField> gradUMeanPtr_ = fvc::grad(UMean_);
    volTensorField& gradUMean_ = gradUMeanPtr_.ref();

    volScalarField k_ = 0.5*tr(UPrime2Mean_);

    calConvectionTerm(UMean_,k_);
    calProductionTerm(UPrime2Mean_,gradUMean_);
    calTurbulenceTransportTerm(turTransPrecursorMean_);
    calViscousTransportTerm(Sij,UPrime,nu_);
    calVPGCorelationTerm(vpgPrecursorMean_);
    calViscousDissipationTerm(Sij,nu_);
}

void Foam::functionObjects::tkeBudget::calConvectionTerm
(
    const volVectorField& UMean,
    const volScalarField& k
)
{
    if(convectionTerm_)
    {
        Info<< "Calculating TKE convection term" << endl;

        volVectorField gradk = fvc::grad(k);

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
                -UMean & gradk
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

void Foam::functionObjects::tkeBudget::calTurbulenceTransportTerm
(
    const volVectorField& turTransPrecursorMean
)
{
    if(turTransportTerm_)
    {
        Info<< "Calculating TKE turbulence transport term" << endl;

        volScalarField divq2UPM = fvc::div(turTransPrecursorMean);

        tmp<volScalarField> tdfPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "turbulenceTransportTerm",
                    turTransPrecursorMean.mesh().time().timeName(),
                    turTransPrecursorMean.mesh()
                ),
                divq2UPM
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_turbulenceTransportTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_turbulenceTransportTerm")) = tdfPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_turbulenceTransportTerm",
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

void Foam::functionObjects::tkeBudget::calViscousTransportTerm
(
    const volSymmTensorField& Sij,
    const volVectorField& UPrime,
    const volScalarField& nu
)
{
    if(visTransportTerm_)
    {
        Info<< "Calculating TKE viscous transport term" << endl;
        volVectorField SUP = Sij & UPrime;
        volScalarField divSUP = fvc::div(SUP);
        tmp<volScalarField> vdfPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "viscousTransportTerm",
                    UPrime.mesh().time().timeName(),
                    UPrime.mesh()
                ),
                2.0*nu*divSUP
            )
        );

        if(obr().foundObject<volScalarField>("tkeBudget_viscousTransportTerm"))
        {
            const_cast<volScalarField&>(obr().lookupObject<volScalarField>("tkeBudget_viscousTransportTerm")) = vdfPtr_;
        }
        else
        {
            obr_.store
            (
                new volScalarField
                (
                    IOobject
                    (
                        "tkeBudget_viscousTransportTerm",
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
    const volVectorField& vpgPrecursorMean
)
{
    if(vpgCorelationTerm_)
    {
        Info<< "Calculating TKE velocity-pressureGradient-Corelation term" << endl;
       
        volScalarField divupPrime = fvc::div(vpgPrecursorMean);
    
        tmp<volScalarField> vpgPtr_
        (
            new volScalarField
            (
                IOobject
                (
                    "VPGCorelationTerm",
                    vpgPrecursorMean.mesh().time().timeName(),
                    vpgPrecursorMean.mesh()
                ),
                -(1.0/rho_)* divupPrime
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
    rho_(dict.get<scalar>("rho")),
    convectionTerm_(dict.get<Switch>("convectionTerm")),
    productionTerm_(dict.get<Switch>("productionTerm")),
    turTransportTerm_(dict.get<Switch>("turbulenceTransportTerm")),
    visTransportTerm_(dict.get<Switch>("viscousTransportTerm")),
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
    checkInsert("tkeBudget_turbulenceTransportTerm",turTransportTerm_,FieldsList);
    checkInsert("tkeBudget_viscousTransportTerm",visTransportTerm_,FieldsList);
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