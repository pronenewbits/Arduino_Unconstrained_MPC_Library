/**************************************************************************************************
 * Class for MPC without constraint
 *
 *  The plant to be controlled is a Linear Time-Invariant System:
 *          x(k+1)  = A*x(k) + B*u(k)   ; x = Nx1, u = Mx1
 *          z(k)    = C*x(k)            ; z = Zx1
 *
 *
 ** Calculate prediction of z(k+1..k+Hp) constants ************************************************
 *
 *      Prediction of state variable of the system:
 *        z(k+1..k+Hp) = (CPSI)*x(k) + (COMEGA)*u(k-1) + (CTHETA)*dU(k..k+Hu-1)         ...{MPC_1}
 *
 *        Constants:
 *          CPSI   = [CA C(A^2) ... C(A^Hp)]'                                       : (Hp*N)xN
 *          COMEGA = [CB C(B+A*B) ... C*Sigma(i=0->Hp-1)A^i*B]'                     : (Hp*N)xM
 *          CTHETA = [         CB                0  ....           0              ]
 *                   [       C(B+A*B)           CB   .             0              ]
 *                   [           .               .    .           CB              ] : (Hp*N)x(Hu*M)
 *                   [           .               .     .           .              ]
 *                   [C*Sigma(i=0->Hp-1)(A^i*B)  .  ....  C*Sigma(i=0->Hp-Hu)A^i*B]
 *
 *
 ** Calculate offline optimization constants ******************************************************
 * 
 *          H       = CTHETA'*Q*CTHETA + R                                              ...{MPC_2}
 * 
 *          XI_FULL = 0.5*H^-1*2*CTHETA'*Q
 *                  = H^-1 * CTHETA' * Q                                                ...{MPC_3}
 * 
 *          XI_DU   = XI_FULL(1:M, :)                                                   ...{MPC_4}
 * 
 *        Constants:
 *          Q     = Weight matrix for set-point deviation   : Hp x Hp
 *          R     = Weight matrix for control signal change : Hu x Hu
 * 
 * 
 ** MPC update algorithm **************************************************************************
 *
 *      Formulation of plant error prediction
 *              E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                ...{MPC_5}
 *
 *      Calculate the optimal control solution:
 *              dU(k)_optimal = XI_DU * E(k)                                            ...{MPC_6}
 *
 *      Integrate the du(k) to get u(k):
 *              u(k) = u(k-1) + du(k)                                                   ...{MPC_7}
 *
 *        Variables:
 *          SP(k) = Set Point vector at time-k              : (Hp*N) x 1
 *          x(k)  = State Variables at time-k               : N x 1
 *          u(k)  = Input plant at time-k                   : M x 1
 * 
 * 
 * See https://github.com/pronenewbits for more!
 *************************************************************************************************/
#include "mpc.h"


MPC::MPC(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR)
{
    vReInit(A, B, C, _bobotQ, _bobotR);
}

void MPC::vReInit(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR)
{
    this->A = A;
    this->B = B;
    this->C = C;
    Q.vSetDiag(_bobotQ);
    R.vSetDiag(_bobotR);
    
    /*  Calculate prediction of z(k+1..k+Hp) constants
     *
     *      Prediction of state variable of the system:
     *        z(k+1..k+Hp) = (CPSI)*x(k) + (COMEGA)*u(k-1) + (CTHETA)*dU(k..k+Hu-1)         ...{MPC_1}
     *
     *        Constants:
     *          CPSI   = [CA C(A^2) ... C(A^Hp)]'                                       : (Hp*N)xN
     *          COMEGA = [CB C(B+A*B) ... C*Sigma(i=0->Hp-1)A^i*B]'                     : (Hp*N)xM
     *          CTHETA = [         CB                0  ....           0              ]
     *                   [       C(B+A*B)           CB   .             0              ]
     *                   [           .               .    .           CB              ] : (Hp*N)x(Hu*M)
     *                   [           .               .     .           .              ]
     *                   [C*Sigma(i=0->Hp-1)(A^i*B)  .  ....  C*Sigma(i=0->Hp-Hu)A^i*B]
     *
     */
    Matrix _Apow(SS_X_LEN, SS_X_LEN);
    /* CPSI     : [ C *   A  ]
     *            [ C *  A^2 ]
     *            [     .    ]                                                   : (Hp*N) x N
     *            [     .    ]
     *            [ C * A^Hp ]
     */
    _Apow = A;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        CPSI = CPSI.InsertSubMatrix((C*_Apow), _i*SS_Z_LEN, 0);
        _Apow = _Apow * A;
    }
    
    /* COMEGA   : [          C * (B)         ]
     *            [        C * (B+A*B)       ]
     *            [             .            ]                                   : (Hp*N) x M
     *            [             .            ]
     *            [ C * Sigma(i=0->Hp-1)A^i*B]
     */
    Matrix _tempSigma(SS_X_LEN, SS_U_LEN);
    _Apow.vSetIdentity();
    _tempSigma = B;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        COMEGA = COMEGA.InsertSubMatrix((C*_tempSigma), _i*SS_Z_LEN, 0);
        _Apow = _Apow * A;
        _tempSigma = _tempSigma + (_Apow*B);
    }
    
    /* CTHETA   : [          C * (B)              0         ....              0             ]
     *            [       C * (B+A*B)           C * (B)      .                0             ]
     *            [            .                  .           .             C * (B)         ]: (Hp*N)x(Hu*M)
     *            [            .                  .            .              .             ]
     *            [C * Sigma(i=0->Hp-1)A^i*B      .         ....  C * Sigma(i=0->Hp-Hu)A^i*B]
     *
     *          : [COMEGA   [0 COMEGA(0:(len(COMEGA)-len(B)),:)]'  ....  [0..0 COMEGA(0:(len(COMEGA)-((Hp-Hu)*len(B))),:)]']
     */
    for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
        CTHETA = CTHETA.InsertSubMatrix(COMEGA, _i*SS_Z_LEN, _i*SS_U_LEN, (MPC_HP_LEN*SS_Z_LEN)-(_i*SS_Z_LEN), SS_U_LEN);
    }
    
    
    /* Calculate the offline optimization constants ---------------------------------------------- */
    Matrix H        {(MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN)};
    Matrix H_INV    {(MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN)};
    Matrix XI       {(MPC_HU_LEN*SS_U_LEN), (MPC_HP_LEN*SS_Z_LEN)};
    
    /*  H       = CTHETA'*Q*CTHETA + R                                                  ...{MPC_2} */
    H = ((CTHETA.Transpose()) * Q * CTHETA) + R;
    
    /*  XI_FULL = 0.5*H^-1*2*CTHETA'*Q
     *          = H^-1 * CTHETA' * Q                                                    ...{MPC_3}
     */
    H_INV = H.Invers();
    if (!H_INV.bMatrixIsValid()) {
        /* set XI as zero to signal that the offline optimization matrix calculation has failed */ 
        XI.vSetToZero();
    } else {
        XI = H_INV * (CTHETA.Transpose()) * Q;
    }
    
    /* XI_DU   = XI_FULL(1:M, :)                                                        ...{MPC_4} */
    XI_DU = XI_DU.InsertSubMatrix(XI, 0, 0, 0, 0, SS_U_LEN, (MPC_HP_LEN*SS_Z_LEN));
}

bool MPC::bUpdate(Matrix &SP, Matrix &x, Matrix &u)
{
    Matrix Err((MPC_HP_LEN*SS_Z_LEN), 1);
    
    /*  E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                        ...{MPC_5} */
    Err = SP - CPSI*x - COMEGA*u;
    
    /*  dU(k)_optimal = XI_DU * E(k)                                                    ...{MPC_6}
     * 
     * Note: If XI_DU initialization is failed in vReInit(), the DU_Out is 
     * always zero (u(k) won't change)
     */
    Matrix DU_Out(SS_U_LEN, 1);
    DU_Out = XI_DU * Err;
    
    /*  u(k) = u(k-1) + du(k)                                                           ...{MPC_7} */
    u = u + DU_Out;
    
    return true;
}
