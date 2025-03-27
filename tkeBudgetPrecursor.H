#ifndef TKE_BUDGET_PRECURSOR
#define TKE_BUDGET_PRECURSOR

#include "fvMeshFunctionObject.H"
#include "volFields.H"

namespace Foam
{
namespace functionObjects
{
    
/*---------------------------------------------------------------------------*\
                    Class tkeBudget Declaration
\*---------------------------------------------------------------------------*/

class tkeBudgetPrecursor
: 
    public fvMeshFunctionObject
{
    protected:

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

        void UPCTermPrecursorField
        (
            const volVectorField& UPrime,
            const volScalarField& pPrime
        );

        void TurTransTermPrecursorField
        (
            const volScalarField& k,
            const volVectorField& UPrime
        );

    public:
    
    //- Runtime type information
    TypeName("tkeBudgetPrecursor");

    // Constructors
    tkeBudgetPrecursor
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~tkeBudgetPrecursor()=default;

    virtual bool execute();

    virtual bool write();
};

}
}

#endif
