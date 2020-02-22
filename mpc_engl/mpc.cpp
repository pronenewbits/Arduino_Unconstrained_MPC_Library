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
 ** MPC update algorithm **************************************************************************
 *
 *      Formulation of plant error prediction
 *              E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                ...{MPC_2}
 *
 * 
 *      Calculate MPC optimization variables:
 *              G = 2*CTHETA'*Q*E(k)                                                    ...{MPC_3}
 *              H = CTHETA'*Q*CTHETA + R                                                ...{MPC_4}
 *
 * 
 *      Formulation of the optimal control problem:
 *
 *              min     dU(k)'*H*dU(k) - G'*dU(k)       ; dU(k) = dU(k..k+Hu-1)
 *             dU(k)
 *
 * 
 *      MPC solution:
 *          (a) For unconstrained MPC:
 *              d [dU(k)'*H*dU(k) - G'*dU(k)]
 *              ----------------------------- = 0   -->   2*H*dU(k)-G = 0   -->   2*H*dU(k) = G
 *                      d[dU(k)]
 *
 *              --> dU(k)_optimal = 1/2 * H^-1 * G                                      ...{MPC_5a}
 *
 *          (b) For constrained MPC (quadrating programming):
 *                  min     dU(k)'*H*dU(k) - G'*dU(k)   ; subject to inequality equation
 *                 dU(k)
 *
 *              --> dU_opt(k) = ActiveSet(2H, -G, ineqLHS, ineqRHS)                     ...{MPC_5b}
 *              --> https://github.com/pronenewbits/Arduino_Constrained_MPC_Library
 *
 *      Integrate the du(k) to get u(k):
 *              u(k) = u(k-1) + du(k)                                                   ...{MPC_6}
 *
 *        Variables:
 *          SP(k) = Set Point vector at time-k              : (Hp*N) x 1
 *          x(k)  = State Variables at time-k               : N x 1
 *          u(k)  = Input plant at time-k                   : M x 1
 *          Q     = Weight matrix for set-point deviation   : Hp x Hp
 *          R     = Weight matrix for control signal change : Hu x Hu
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
}

bool MPC::bUpdate(Matrix &SP, Matrix &x, Matrix &u)
{
    Matrix Err((MPC_HP_LEN*SS_Z_LEN), 1);
    Matrix G((MPC_HU_LEN*SS_U_LEN), 1);
    Matrix H((MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN));
    
    /*  E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                        ...{MPC_2} */
    Err = SP - CPSI*x - COMEGA*u;
    
    /*  G = 2*CTHETA'*Q*E(k)                                                            ...{MPC_3} */
    G = 2.0 * (CTHETA.Transpose()) * Q * Err;
    
    /*  H = CTHETA'*Q*CTHETA + R                                                        ...{MPC_4} */
    H = ((CTHETA.Transpose()) * Q * CTHETA) + R;
    
    /*  --> dU(k)_optimal = 1/2 * H^-1 * G                                              ...{MPC_5a} */
    Matrix H_inv = H.Invers();
    if (!H_inv.bMatrixIsValid()) {
        /* return false; */
        DU.vSetToZero();
        
        return false;
    } else {
        DU = (H_inv) * G * 0.5;
    }
    
    /*  u(k) = u(k-1) + du(k)                                                           ...{MPC_6} */
    Matrix DU_Out(SS_U_LEN, 1);
    for (int32_t _i = 0; _i < SS_U_LEN; _i++) {
        DU_Out[_i][0] = DU[_i][0];
    }
    u = u + DU_Out;
    
    return true;
}
