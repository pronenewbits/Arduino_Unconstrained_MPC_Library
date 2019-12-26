#ifndef MPC_H
#define MPC_H

#include "konfig.h"
#include "matrix.h"


#if (MPC_HP_LEN < MPC_HU_LEN)
    #error("The MPC_HP_LEN must be more than or equal MPC_HU_LEN!");
#endif
#if (((MPC_HP_LEN*SS_Z_LEN) > MATRIX_MAXIMUM_SIZE) || ((MPC_HP_LEN*SS_X_LEN) > MATRIX_MAXIMUM_SIZE) || ((MPC_HU_LEN*SS_U_LEN) > MATRIX_MAXIMUM_SIZE))
    #error("The MATRIX_MAXIMUM_SIZE is too small to do MPC calculation!");
#endif

class MPC
{
public:
    MPC(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR);
    void vReInit(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR);
    bool bUpdate(Matrix &SP, Matrix &X, Matrix &U);

protected:
    void bCalculateActiveSet(void);

private:
    Matrix CPSI     {(MPC_HP_LEN*SS_Z_LEN), SS_X_LEN};
    Matrix COMEGA   {(MPC_HP_LEN*SS_Z_LEN), SS_U_LEN};
    Matrix CTHETA   {(MPC_HP_LEN*SS_Z_LEN), (MPC_HU_LEN*SS_U_LEN)};

    Matrix DU       {(MPC_HU_LEN*SS_U_LEN), 1};

    Matrix A        {SS_X_LEN, SS_X_LEN};
    Matrix B        {SS_X_LEN, SS_U_LEN};
    Matrix C        {SS_Z_LEN, SS_X_LEN};

    Matrix Q        {(MPC_HP_LEN*SS_Z_LEN), (MPC_HP_LEN*SS_Z_LEN)};
    Matrix R        {(MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN)};
    Matrix SQ       {(MPC_HP_LEN*SS_Z_LEN), (MPC_HP_LEN*SS_Z_LEN)};
    Matrix SR       {(MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN)};
    
    Matrix XI_DU    {(SS_U_LEN), (MPC_HP_LEN*SS_Z_LEN)};
};



#endif // MPC_H
