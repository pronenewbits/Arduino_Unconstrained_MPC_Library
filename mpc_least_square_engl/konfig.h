#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>



/* State Space dimension */
#define SS_X_LEN    (4)
#define SS_Z_LEN    (2)
#define SS_U_LEN    (2)
#define SS_DT_MILIS (20)                            /* 20 ms */
#define SS_DT       float_prec(SS_DT_MILIS/1000.)   /* Sampling time */

#define MATRIX_MAXIMUM_SIZE     (28)        /* Change this size based on the biggest matrix you will use */

#define PAKAI_FLOAT         1
#define PAKAI_DOUBLE        2
#define PRESISI_FPU         (PAKAI_FLOAT)

#if (PRESISI_FPU == PAKAI_FLOAT)
    #define float_prec      float
    #define float_prec_ZERO (1e-8)
#elif (PRESISI_FPU == PAKAI_DOUBLE)
    #define float_prec      double
    #define float_prec_ZERO (1e-15)
#else
    #error("Kepresisian FPU belum didefinisi!");
#endif

#define MPC_HP_LEN      (7)
#define MPC_HU_LEN      (4)

/* Aktifkan definisi ini untuk implementasi kode di MCU (std lib tidak digunakan) */
#define SISTEM_PC                   1
#define SISTEM_EMBEDDED_NO_PRINT    2
#define SISTEM_EMBEDDED_ARDUINO     3

#define SISTEM_IMPLEMENTASI         SISTEM_EMBEDDED_ARDUINO


#define MATRIX_PAKAI_BOUND_CHECKING

#endif // MAIN_H
