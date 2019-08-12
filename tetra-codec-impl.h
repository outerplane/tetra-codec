#pragma once

#include <stdint.h>

//------------------------------------------------------------------------------

#ifdef __GNUC__
#pragma GCC visibility push(hidden)
#endif

//------------------------------------------------------------------------------

#define TSC_pp      (int16_t)10

#define TSC_L_next  (int16_t)40
#define TSC_L_frame  (int16_t)240
#define TSC_L_subfr  (int16_t)60
#define TSC_p        (int16_t)10
#define TSC_pp1      (int16_t)11
#define TSC_L_total  (int16_t)(TSC_L_frame+TSC_L_next+TSC_p)
#define TSC_dim_rr   (int16_t)32
#define TSC_pit_max  (int16_t)143
#define TSC_L_inter  (int16_t)15

#define  TSC_parm_size (int16_t)23

//------------------------------------------------------------------------------

struct tetra_codec {

	// SUB_DSP.C
	int16_t old_A[TSC_pp+1];

	// SCOD_TET.C & SDEC_TET.C
	int16_t old_speech[TSC_L_total];
	int16_t *speech;
	int16_t *p_window;
	int16_t *new_speech;
	int16_t old_wsp[TSC_L_frame+TSC_pit_max];
	int16_t *wsp;
	int16_t old_exc[TSC_L_frame+TSC_pit_max+TSC_L_inter];
	int16_t *exc;
	int16_t ai_zero[TSC_L_subfr+TSC_pp1];
	int16_t *zero;
	int16_t F_gamma1[TSC_p];
	int16_t F_gamma2[TSC_p];
	int16_t F_gamma3[TSC_p];
	int16_t F_gamma4[TSC_p];
	int16_t lspold[TSC_p];
	int16_t lspnew[TSC_p];
	int16_t lspnew_q[TSC_p];
	int16_t lspold_q[TSC_p];
	int16_t mem_syn[TSC_p];
	int16_t mem_w0[TSC_p];
	int16_t mem_w[TSC_p];
	int16_t rr[TSC_dim_rr][TSC_dim_rr];
	int16_t old_parm[TSC_parm_size];
	int16_t old_T0;
	int16_t last_ener_cod;
	int16_t last_ener_pit;

	// SUB_SC_D.C
	int16_t y_hi;
	int16_t y_lo;
	int16_t x0;
	int16_t lsp_old[10];

	// TETRA_OP.C
	int Overflow;
	int Carry;
};

//------------------------------------------------------------------------------

typedef struct tetra_codec tetra_codec;

typedef int16_t Word16;
typedef int32_t Word32;

//------------------------------------------------------------------------------

void init_tetra_codec(tetra_codec* st);

void Init_Pre_Process(tetra_codec* st);
void Init_Coder_Tetra(tetra_codec* st);
void Pre_Process(tetra_codec* st, Word16 signal[], Word16 lg);
void Coder_Tetra(tetra_codec* st, Word16 ana[], Word16 synth[]);
void Prm2bits_Tetra(tetra_codec* st, Word16 prm[], Word16 bits[]);

void Bits2prm_Tetra(tetra_codec* st, Word16 bits[], Word16 prm[]);
void Init_Decod_Tetra(tetra_codec* st);
void Decod_Tetra(tetra_codec* st, Word16 parm[], Word16 synth[]);
void Post_Process(tetra_codec* st, Word16 signal[], Word16 lg);

//------------------------------------------------------------------------------

#ifdef __GNUC__
#pragma GCC visibility pop
#endif

//------------------------------------------------------------------------------
