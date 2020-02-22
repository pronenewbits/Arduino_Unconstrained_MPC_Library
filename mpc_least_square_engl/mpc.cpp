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
 *      Recreate the optimal control problem solution:
 *          [(SQ * CTHETA)] * dU(k)_optimal = [(SQ*E(k)]
 *          [      SR     ]                   [    0   ]
 * 
 *          GammaLeft * dU(k)_optimal = GammaRight
 * 
 *          Q_L * R_L = GammaLeft                                                       ...{MPC_2}
 * 
 *        Constants:
 *          SQ    = Square root of Weight matrix for set-point deviation    : Hp x Hp
 *          SR    = Square root of Weight matrix for control signal change  : Hu x Hu
 *          Q_L   = Orthogonal matrix of QR Decomposition of GammaLeft      : (Hp*Z+Hu*M) x  (Hp*Z+Hu*M)
 *          R_L   = Upper triangular matrix of QR Decomposition of GammaLeft: (Hp*Z+Hu*M) x  (Hp*Z)
 * 
 * 
 ** MPC update algorithm **************************************************************************
 *
 *      Formulation of plant error prediction
 *          E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                    ...{MPC_3}
 *
 *      Construct the optimal control solution equation:
 *          R_L * dU(k)_optimal = Qt_L * [(SQ*E(k)]                                     ...{MPC_4}
 *                                       [    0   ]
 *              
 *      Calculate the optimal control solution using back-subtitution:
 *          R_L * dU(k)_optimal = BackSubRight                                          ...{MPC_5}
 *
 *      Integrate the du(k) to get u(k):
 *          u(k) = u(k-1) + du(k)                                                       ...{MPC_6}
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
    SQ.vSetDiag(sqrt(_bobotQ));
    SR.vSetDiag(sqrt(_bobotR));
    
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
    
    
    /* Calculate offline optimization constants
     * 
     *      Recreate the optimal control problem solution:
     *          [(SQ * CTHETA)] * dU(k)_optimal = [(SQ*E(k)]
     *          [      SR     ]                   [    0   ]
     * 
     *          GammaLeft * dU(k)_optimal = GammaRight
     * 
     *          Q_L * R_L = GammaLeft                                                       ...{MPC_2}
     * 
     * NOTE: QRDec function return the transpose of Q (i.e. Q').
     */
    Matrix GammaLeft((MPC_HP_LEN*SS_Z_LEN + MPC_HU_LEN*SS_U_LEN), MPC_HP_LEN*SS_Z_LEN);
    GammaLeft = GammaLeft.InsertSubMatrix((SQ * CTHETA), 0, 0);
    GammaLeft = GammaLeft.InsertSubMatrix(SR, MPC_HP_LEN*SS_Z_LEN, 0);
    GammaLeft.QRDec(Qt_L, R_L);
}

bool MPC::bUpdate(Matrix &SP, Matrix &x, Matrix &u)
{
    Matrix Err((MPC_HP_LEN*SS_Z_LEN), 1);
    
    /*  E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                            ...{MPC_3} */
    Err = SP - CPSI*x - COMEGA*u;
    

    if (!Qt_L.bMatrixIsValid()) {
        /* The QR Decomposition in the initialization step has failed, return false */
        DU.vSetToZero();
        
        return false;
    } else {
        /*      Construct the optimal control solution equation:
         *          R_L * dU(k)_optimal = Qt_L * [(SQ*E(k)]                                 ...{MPC_4}
         *                                       [    0   ]
         * 
         * NOTE: We only need the first (Hp*Z)-th columns of Qt_L to construct the 
         *          right hand equation (encapsulated in Qt_LSQE variable).
         */
        Matrix Q1((MPC_HP_LEN*SS_Z_LEN + MPC_HU_LEN*SS_U_LEN), MPC_HP_LEN*SS_Z_LEN);
        Q1 = Q1.InsertSubMatrix(Qt_L, 0, 0, 0, 0, (MPC_HP_LEN*SS_Z_LEN + MPC_HU_LEN*SS_U_LEN), MPC_HP_LEN*SS_Z_LEN);

        Matrix Qt_LSQE(((MPC_HP_LEN*SS_Z_LEN + MPC_HU_LEN*SS_U_LEN)), 1);
        Qt_LSQE = Q1*SQ*Err;
        
        /* The linear equation is overdetermined, just need the first (Hu*M)-th row */
        Matrix BackSubRight((MPC_HU_LEN*SS_U_LEN), 1);
        BackSubRight = BackSubRight.InsertSubMatrix(Qt_LSQE, 0, 0, 0, 0, MPC_HU_LEN*SS_U_LEN, 1);
        /* The linear equation is overdetermined, just need the first (Hu*M)-th row */
        Matrix R1(MPC_HU_LEN*SS_U_LEN, MPC_HU_LEN*SS_U_LEN);
        R1 = R1.InsertSubMatrix(R_L, 0, 0, 0, 0, MPC_HU_LEN*SS_U_LEN, MPC_HU_LEN*SS_U_LEN);
        
        
        /*      Calculate the optimal control solution using back-subtitution:
         *          R_L * dU(k)_optimal = BackSubRight                                      ...{MPC_5}
         */
        DU = R1.BackSubtitution(R1, BackSubRight);
    }

    /*      Integrate the du(k) to get u(k):
     *          u(k) = u(k-1) + du(k)                                                       ...{MPC_6}
     */
    Matrix DU_Out(SS_U_LEN, 1);
    for (int32_t _i = 0; _i < SS_U_LEN; _i++) {
        DU_Out[_i][0] = DU[_i][0];
    }
    u = u + DU_Out;
    
    return true;
}
