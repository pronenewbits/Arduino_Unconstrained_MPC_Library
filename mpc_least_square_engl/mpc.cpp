/**************************************************************************************************
 * Class for MPC without constraint
 *
 *  The plant to be controlled is a Linear Time-Invariant System:
 *          x(k+1)  = A*x(k) + B*u(k)   ; x = Nx1, u = Mx1
 *          z(k)    = C*x(k)            ; z = Zx1
 *
 *
 ** Calculate prediction of x(k+1..k+Hp) constants ************************************************
 *
 *      Prediction of state variable of the system:
 *        x(k+1..k+Hp) = PSI*x(k) + OMEGA*u(k-1) + THETA*dU(k..k+Hu-1)                  ...{MPC_1}
 *
 *        Constants:
 *          PSI   = [A A^2 ... A^Hp]'                                         : (Hp*N)xN
 *          OMEGA = [B B+A*B ... Sigma(i=0->Hp-1)A^i*B]'                      : (Hp*N)xM
 *          THETA = [         B               0  ....           0            ]
 *                  [       B+A*B             B   .             0            ]
 *                  [         .               .    .            B            ]: (Hp*N)x(Hu*M)
 *                  [         .               .     .           .            ]
 *                  [Sigma(i=0->Hp-1)(A^i*B)  .  ....  Sigma(i=0->Hp-Hu)A^i*B]
 *
 *
 ** Calculate prediction of z(k+1..k+Hp) constants ************************************************
 *
 *      Prediction of output of the system:
 *        z(k+1..k+Hp) = (Cz*PSI)*x(k) + (Cz*OMEGA)*u(k-1) + (Cz*THETA)*dU(k..k+Hu-1)   ...{MPC_2}
 *
 *        Constants:
 *          Cz      : [C 0 0 .. 0]    ; 0 = zero(ZxN)       : (Hp*Z)x(Hp*N)
 *                    [0 C 0 .. 0]
 *                    [0 0 C .. 0]
 *                    [. . . .. 0]
 *                    [0 0 0 .. C]
 *
 *          CPSI   = Cz*PSI                                 : (Hp*Z)xN
 *          COMEGA = Cz*OMEGA                               : (Hp*Z)xM
 *          CTHETA = Cz*THETA                               : (Hp*Z)x(Hu*M)
 *
 *
 ** Calculate offline optimization constants ******************************************************
 * 
 *      Recreate the optimal control problem solution:
 *          [(SQ * CTHETA)] * dU(k)_optimal = [(SQ*E(k)]
 *          [      SR     ]                   [    0   ]
 * 
 *          GammaLeft * dU(k)_optimal = GammaRight
 * 
 *          Q_L * R_L = GammaLeft                                                       ...{MPC_3}
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
 *          E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                    ...{MPC_4}
 *
 *      Construct the optimal control solution equation:
 *          R_L * dU(k)_optimal = Qt_L * [(SQ*E(k)]                                     ...{MPC_5}
 *                                       [    0   ]
 *              
 *      Calculate the optimal control solution using back-subtitution:
 *          R_L * dU(k)_optimal = BackSubRight                                          ...{MPC_6}
 *
 *      Integrate the du(k) to get u(k):
 *          u(k) = u(k-1) + du(k)                                                       ...{MPC_7}
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
    SQ.vIsiDiagonal(sqrt(_bobotQ));
    SR.vIsiDiagonal(sqrt(_bobotR));
    
    /*  Calculate prediction of x(k+1..k+Hp) constants
     *
     *      Prediction of state variable of the system:
     *        x(k+1..k+Hp) = PSI*x(k) + OMEGA*u(k-1) + THETA*dU(k..k+Hu-1)                  ...{MPC_1}
     *
     *        Constants:
     *          PSI   = [A A^2 ... A^Hp]'                                         : (Hp*N)xN
     *          OMEGA = [B B+A*B ... Sigma(i=0->Hp-1)A^i*B]'                      : (Hp*N)xM
     *          THETA = [         B               0  ....           0            ]
     *                  [       B+A*B             B   .             0            ]
     *                  [         .               .    .            B            ]: (Hp*N)x(Hu*M)
     *                  [         .               .     .           .            ]
     *                  [Sigma(i=0->Hp-1)(A^i*B)  .  ....  Sigma(i=0->Hp-Hu)A^i*B]
     *
     */
    Matrix _PSI     ((MPC_HP_LEN*SS_X_LEN), SS_X_LEN);
    Matrix _OMEGA   ((MPC_HP_LEN*SS_X_LEN), SS_U_LEN);
    Matrix _THETA   ((MPC_HP_LEN*SS_X_LEN), (MPC_HU_LEN*SS_U_LEN));
    
    Matrix _Apow(SS_X_LEN, SS_X_LEN);
    /* PSI      : [  A  ]
     *            [ A^2 ]
     *            [  .  ]                                                   : (Hp*N) x N
     *            [  .  ]
     *            [A^Hp ]
     */
    _Apow = A;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        _PSI = _PSI.InsertSubMatrix(_Apow, _i*SS_X_LEN, 0);
        _Apow = _Apow * A;
    }
    
    /* OMEGA    : [          B          ]
     *            [        B+A*B        ]
     *            [          .          ]                                   : (Hp*N) x M
     *            [          .          ]
     *            [Sigma(i=0->Hp-1)A^i*B]
     */
    Matrix _tempSigma(SS_X_LEN, SS_U_LEN);
    _Apow.vSetIdentitas();
    _tempSigma = B;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        _OMEGA = _OMEGA.InsertSubMatrix(_tempSigma, _i*SS_X_LEN, 0);
        _Apow = _Apow * A;
        _tempSigma = _tempSigma + (_Apow*B);
    }
    
    /* THETA    : [         B               0  ....           0            ]
     *            [       B+A*B             B   .             0            ]
     *            [         .               .    .            B            ]: (Hp*N)x(Hu*M)
     *            [         .               .     .           .            ]
     *            [Sigma(i=0->Hp-1)A^i*B    .  ....  Sigma(i=0->Hp-Hu)A^i*B]
     *
     *          : [OMEGA   [0 OMEGA(0:(len(OMEGA)-len(B)),:)]'  ....  [0..0 OMEGA(0:(len(OMEGA)-((Hp-Hu)*len(B))),:)]']
     */
    for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
        _THETA = _THETA.InsertSubMatrix(_OMEGA, _i*SS_X_LEN, _i*SS_U_LEN, (MPC_HP_LEN*SS_X_LEN)-(_i*SS_X_LEN), SS_U_LEN);
    }
    
    
    
    /* Calculate prediction of z(k+1..k+Hp) constants
     *
     *      Prediction of output of the system:
     *        z(k+1..k+Hp) = (Cz*PSI)*x(k) + (Cz*OMEGA)*u(k-1) + (Cz*THETA)*dU(k..k+Hu-1)   ...{MPC_2}
     *
     *        Constants:
     *          Cz      : [C 0 0 .. 0]    ; 0 = zero(ZxN)       : (Hp*Z)x(Hp*N)
     *                    [0 C 0 .. 0]
     *                    [0 0 C .. 0]
     *                    [. . . .. 0]
     *                    [0 0 0 .. C]
     *
     *          CPSI   = Cz*PSI                                 : (Hp*Z)xN
     *          COMEGA = Cz*OMEGA                               : (Hp*Z)xM
     *          CTHETA = Cz*THETA                               : (Hp*Z)x(Hu*M)
     */
    Matrix Cz(MPC_HP_LEN*SS_Z_LEN, MPC_HP_LEN*SS_X_LEN);
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        Cz = Cz.InsertSubMatrix(C, _i*SS_Z_LEN, _i*SS_X_LEN);
    }
    CPSI    = Cz * _PSI;
    COMEGA  = Cz * _OMEGA;
    CTHETA  = Cz * _THETA;
    
    
    
    /* Calculate offline optimization constants
     * 
     *      Recreate the optimal control problem solution:
     *          [(SQ * CTHETA)] * dU(k)_optimal = [(SQ*E(k)]
     *          [      SR     ]                   [    0   ]
     * 
     *          GammaLeft * dU(k)_optimal = GammaRight
     * 
     *          Q_L * R_L = GammaLeft                                                       ...{MPC_3}
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
    
    /*  E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                            ...{MPC_4} */
    Err = SP - CPSI*x - COMEGA*u;
    

    if (!Qt_L.bCekMatrixValid()) {
        /* The QR Decomposition in the initialization step has failed, return false */
        DU.vIsiNol();
        
        return false;
    } else {
        /*      Construct the optimal control solution equation:
         *          R_L * dU(k)_optimal = Qt_L * [(SQ*E(k)]                                     ...{MPC_5}
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
         *          R_L * dU(k)_optimal = BackSubRight                                          ...{MPC_6}
         */
        DU = R1.BackSubtitution(R1, BackSubRight);
    }

    /*      Integrate the du(k) to get u(k):
     *          u(k) = u(k-1) + du(k)                                                       ...{MPC_7}
     */
    Matrix DU_Out(SS_U_LEN, 1);
    for (int32_t _i = 0; _i < SS_U_LEN; _i++) {
        DU_Out[_i][0] = DU[_i][0];
    }
    u = u + DU_Out;
    
    return true;
}
