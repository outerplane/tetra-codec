
#include "tetra-codec-impl.h"

// ETSI EN 300 395-2 V1.3.1 (2005-01)
// Part 2: TETRA codec

//------------------------------------------------------------------------------

static const tetra_codec default_codec_state = {
	// SUB_DSP.C
	.old_A = { 4096,0,0,0,0,0,0,0,0,0,0 },

	// SCOD_TET.C & SDEC_TET.C
	.lspold = { 30000, 26000, 21000, 15000, 8000, 0, -8000,-15000,-21000,-26000 },
	.last_ener_cod = 0,
	.last_ener_pit = 0,

	// SUB_SC_D.C
	.lsp_old = { 30000, 26000, 21000, 15000, 8000, 0, -8000,-15000,-21000,-26000 },

	// TETRA_OP.C
	.Overflow = 0,
	.Carry = 0
};


void init_tetra_codec(tetra_codec* st)
{
	*st = default_codec_state;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- SOURCE.H ------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


/************************************************************************
*
*	FILENAME		:	source.h
*
*	DESCRIPTION		:	DEFINITIONS
*
*					FUNCTION PROTOTYPES for the TETRA
*					speech source coder and decoder
*
*					OPERATOR PROTOTYPES for TETRA codec
*
************************************************************************/

/************************************************************************
*
*      DEFINITIONS
*
************************************************************************/

//#ifndef TYPEDEF_H
//#define TYPEDEF_H


//typedef short Word16;
//typedef long  Word32;
//typedef int   Flag;

//#endif


/************************************************************************
*
*	DESCRIPTION		:	FUNCTION PROTOTYPES for the TETRA
*					speech source coder and decoder
*
************************************************************************/

//#ifndef TETRA_H
//#define TETRA_H


/* Basic functions */

static Word32 add_sh(tetra_codec* st, Word32 L_var, Word16 var1, Word16 shift);
static Word32 add_sh16(tetra_codec* st, Word32 L_var, Word16 var1);
static Word16 bin2int(tetra_codec* st, Word16 no_of_bits, Word16 *bitstream);
static void   int2bin(tetra_codec* st, Word16 value, Word16 no_of_bits, Word16 *bitstream);
static Word32 Load_sh(tetra_codec* st, Word16 var1, Word16 shift);
static Word32 Load_sh16(tetra_codec* st, Word16 var1);
static Word32 norm_v(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 *var2);
static Word16 store_hi(tetra_codec* st, Word32 L_var1, Word16 var2);
static Word32 sub_sh(tetra_codec* st, Word32 L_var, Word16 var1, Word16 shift);
static Word32 sub_sh16(tetra_codec* st, Word32 L_var, Word16 var1);

/* Extended precision functions */

static Word32 div_32(tetra_codec* st, Word32 L_num, Word16 hi, Word16 lo);
static Word32 L_comp(tetra_codec* st, Word16 hi, Word16 lo);
static void   L_extract(tetra_codec* st, Word32 L_32, Word16 *hi, Word16 *lo);
static Word32 mpy_mix(tetra_codec* st, Word16 hi1, Word16 lo1, Word16 lo2);
static Word32 mpy_32(tetra_codec* st, Word16 hi1, Word16 lo1, Word16 hi2, Word16 lo2);


/* Mathematic functions  */

static Word32 inv_sqrt(tetra_codec* st, Word32 L_x);
static void   Log2(tetra_codec* st, Word32 L_x, Word16 *exponant, Word16 *fraction);
static Word32 pow2(tetra_codec* st, Word16 exponant, Word16 fraction);

/* General signal processing */

static void   Autocorr(tetra_codec* st, Word16 x[], Word16 p, Word16 r_h[], Word16 r_l[]);
static void   Az_Lsp(tetra_codec* st, Word16 a[], Word16 lsp[], Word16 old_lsp[]);
static void   Back_Fil(tetra_codec* st, Word16 x[], Word16 h[], Word16 y[], Word16 L);
static Word16 Chebps(tetra_codec* st, Word16 x, Word16 f[], Word16 n);
static void   Convolve(tetra_codec* st, Word16 x[], Word16 h[], Word16 y[], Word16 L);
static void   Fac_Pond(tetra_codec* st, Word16 gamma, Word16 fac[]);
static void   Get_Lsp_Pol(tetra_codec* st, Word16 *lsp, Word32 *f);
static void   Int_Lpc4(tetra_codec* st, Word16 lsp_old[], Word16 lsp_new[], Word16 a_4[]);
static void   Lag_Window(tetra_codec* st, Word16 p, Word16 r_h[], Word16 r_l[]);
static void   Levin_32(tetra_codec* st, Word16 Rh[], Word16 Rl[], Word16 A[]);
static Word32 Lpc_Gain(tetra_codec* st, Word16 a[]);
static void   Lsp_Az(tetra_codec* st, Word16 lsp[], Word16 a[]);
static void   Pond_Ai(tetra_codec* st, Word16 a[], Word16 fac[], Word16 a_exp[]);
static void   Residu(tetra_codec* st, Word16 a[], Word16 x[], Word16 y[], Word16 lg);
static void   Syn_Filt(tetra_codec* st, Word16 a[], Word16 x[], Word16 y[], Word16 lg, Word16 mem[],
                Word16 update);

/* Specific coder functions */

//void   Init_Coder_Tetra(tetra_codec* st);
//void   Coder_Tetra(tetra_codec* st, Word16 ana[], Word16 synth[]);
static void   Cal_Rr2(tetra_codec* st, Word16 h[], Word16 *rr);
static void   Clsp_334(tetra_codec* st, Word16 *lsp, Word16 *lsp_q, Word16 *indice);
static Word16 D4i60_16(tetra_codec* st, Word16 dn[], Word16 f[], Word16 h[], Word16 rr[][32],
                Word16 cod[], Word16 y[], Word16 *sign, Word16 *shift_code);
static Word16 Ener_Qua(tetra_codec* st, Word16 A[], Word16 prd_lt[], Word16 code[], Word16 L_subfr,
                Word16 *gain_pit, Word16 *gain_cod);
static Word16 G_Code(tetra_codec* st, Word16 xn2[], Word16 y2[], Word16 L_subfr);
static Word16 G_Pitch(tetra_codec* st, Word16 xn[], Word16 y1[], Word16 L_subfr);
//void   Init_Pre_Process(tetra_codec* st);
static Word16 Lag_Max(tetra_codec* st, Word16 signal[], Word16 sig_dec[], Word16 L_frame,
               Word16 lag_max, Word16 lag_min, Word16 *cor_max);
static Word16 Pitch_Fr(tetra_codec* st, Word16 exc[], Word16 xn[], Word16 h[], Word16 L_subfr,
                Word16 t0_min, Word16 t0_max, Word16 i_subfr,
		    Word16 *pit_frac);
static Word16 Pitch_Ol_Dec(tetra_codec* st, Word16 signal[], Word16 L_frame);
static void   Pred_Lt(tetra_codec* st, Word16 exc[], Word16 T0, Word16 frac, Word16 L_subfr);
//void   Pre_Process(tetra_codec* st, Word16 signal[], Word16 lg);
//void   Prm2bits_Tetra(tetra_codec* st, Word16 prm[], Word16 bits[]);

/* Specific decoder functions */

//void   Init_Decod_Tetra(tetra_codec* st);
//void   Decod_Tetra(tetra_codec* st, Word16 parm[], Word16 synth[]);
//void   Bits2prm_Tetra(tetra_codec* st, Word16 bits[], Word16 prm[]);
static Word16 Dec_Ener(tetra_codec* st, Word16 index, Word16 bfi, Word16 A[], Word16 prd_lt[],
	    Word16 code[], Word16 L_subfr, Word16 *gain_pit, Word16 *gain_cod);
static void   D_D4i60(tetra_codec* st, Word16 index,Word16 sign,Word16 shift, Word16 F[],
	    Word16 cod[]);
static void   D_Lsp334(tetra_codec* st, Word16 indice[], Word16 lsp[], Word16 old_lsp[]);
//void   Post_Process(tetra_codec* st, Word16 signal[], Word16 lg);

//#endif



/************************************************************************
*
*	DESCRIPTION		:     OPERATOR PROTOTYPES for TETRA codec
*
************************************************************************/

//#ifndef TETRA_OP_H
//#define TETRA_OP_H


/*-----------------------*
 * Constants and Globals *
 *-----------------------*/


//extern Flag Overflow;
//extern Flag Carry;

#define MAX_32 (Word32)0x7fffffff
#define MIN_32 (Word32)0x80000000

#define MAX_16 (Word16)0x7fff
#define MIN_16 (Word16)0x8000


/*-----------------------*
 * Operators prototypes  *
 *-----------------------*/

static Word16 abs_s(Word16 var1);                /* Short abs,           1 */
static Word16 add(tetra_codec* st, Word16 var1, Word16 var2);     /* Short add,           1 */
static Word16 div_s(tetra_codec* st, Word16 var1, Word16 var2);   /* Short division,     18 */
static Word16 extract_h(Word32 L_var1);          /* Extract high,        1 */
static Word16 extract_l(Word32 L_var1);          /* Extract low,         1 */
static Word16 mult(tetra_codec* st, Word16 var1, Word16 var2);    /* Short mult,          1 */
static Word16 mult_r(tetra_codec* st, Word16 var1, Word16 var2);  /* Mult with round,     2 */
static Word16 negate(Word16 var1);               /* Short negate,        1 */
static Word16 norm_l(Word32 L_var1);             /* Long norm,          30 */
static Word16 norm_s(Word16 var1);               /* Short norm,         15 */
static Word16 tsc_round(tetra_codec* st, Word32 L_var1);              /* Round,               1 */
static Word16 shl(tetra_codec* st, Word16 var1, Word16 var2);     /* Short shift left,    1 */
static Word16 shr(tetra_codec* st, Word16 var1, Word16 var2);     /* Short shift right,   1 */
static Word16 sub(tetra_codec* st, Word16 var1, Word16 var2);     /* Short sub,           1 */

static Word32 L_abs(Word32 L_var1);              /* Long abs,            3 */
static Word32 L_add(tetra_codec* st, Word32 L_var1, Word32 L_var2);  /* Long add,         2 */
static Word32 L_deposit_h(Word16 var1);          /* 16 bit var1 -> MSB   2 */
static Word32 L_deposit_l(Word16 var1);          /* 16 bit var1 -> LSB,  2 */
static Word32 L_mac(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2);  /* Mac,   1 */
static Word32 L_mac0(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2); /* no shi 1 */
static Word32 L_msu(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2);  /* Msu,   1 */
static Word32 L_msu0(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2); /* no shi 1 */
static Word32 L_mult(tetra_codec* st, Word16 var1, Word16 var2);  /* Long mult,           1 */
static Word32 L_mult0(Word16 var1, Word16 var2); /* Long mult no shift,  1 */
static Word32 L_negate(Word32 L_var1);           /* Long negate,         2 */
static Word32 L_shl(tetra_codec* st, Word32 L_var1, Word16 var2); /* Long shift left,     2 */
static Word32 L_shr(tetra_codec* st, Word32 L_var1, Word16 var2); /* Long shift right,    2 */
static Word32 L_shr_r(tetra_codec* st, Word32 L_var1, Word16 var2);  /* L_shr with round, 3 */
static Word32 L_sub(tetra_codec* st, Word32 L_var1, Word32 L_var2);  /* Long sub,         2 */

//#endif




//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- FBAS_TET.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------







/************************************************************************
*
*	FILENAME		:	fbas_tet.c
*
*	DESCRIPTION		:	Library of basic functions used in the TETRA speech
*					codec, other than accepted operators
*
************************************************************************
*
*	FUNCTIONS		:	- add_sh(st, )
*					- add_sh16()
*					- bin2int(st, )
*					- int2bin(st, )
*					- Load_sh(st, )
*					- Load_sh16()
*					- norm_v()
*					- store _hi()
*					- sub_sh()
*					- sub_sh16()
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*
************************************************************************/

//#include "source.h"

static const Word16 POW2[16] = { -1, -2, -4, -8, -16, -32, -64, -128, -256, -512,
                           -1024, -2048, -4096, -8192, -16384, -32768};

/************************************************************************
*
*	Function Name : add_sh
*
*	Purpose :
*
*		Add var1 with a left shift(0-15) to L_var2.
*		Control saturation and set overflow flag
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var2
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var2 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		shift
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0 <= shift <= 15.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 add_sh(tetra_codec* st, Word32 L_var2, Word16 var1, Word16 shift)
{
	return( L_msu0(st, L_var2, var1, POW2[shift]));
}


/************************************************************************
*
*	Function Name : add_sh16
*
*	Purpose :
*
*		Add var1 with a left shift of 16 to L_var2.
*		Control saturation and set overflow flag
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var2
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var2 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 add_sh16(tetra_codec* st, Word32 L_var2, Word16 var1)
{
	return( L_msu(st, L_var2, var1, (Word16)-32768));
}


/************************************************************************
*
*	Function Name : bin2int
*
*	Purpose :
*
*		Read "no_of_bits" bits from the array bitstream[] and convert to integer
*
*	Inputs :
*
*		no_of_bits
*			16 bit
*
*		*bitstream
*			16 bit
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		value
*			32 bit
*
************************************************************************/

#define TSC_BIT_1  1

Word16 bin2int(tetra_codec* st, Word16 no_of_bits, Word16 *bitstream)
{
   Word16 value, i, bit;

   value = 0;
   for (i = 0; i < no_of_bits; i++)
   {
     value = shl(st,  value,(Word16)1 );
     bit = *bitstream++;
     if (bit == TSC_BIT_1)  value += 1;
   }
   return(value);
}

/************************************************************************
*
*	Function Name : int2bin
*
*	Purpose :
*
*		Convert integer to binary and write the bits to the array bitstream []
*
*	Inputs :
*
*		no_of_bits
*			16 bit
*
*		*bitstream
*			16 bit
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		value
*			32 bit
*
************************************************************************/

#define TSC_BIT_0     0
#define TSC_BIT_1     1
#define TSC_MASK      1

void int2bin(tetra_codec* st, Word16 value, Word16 no_of_bits, Word16 *bitstream)
{
   Word16 *pt_bitstream, i, bit;

   pt_bitstream = bitstream + no_of_bits;

   for (i = 0; i < no_of_bits; i++)
   {
     bit = value & TSC_MASK;
     if (bit == 0)
         *--pt_bitstream = TSC_BIT_0;
     else
         *--pt_bitstream = TSC_BIT_1;
     value = shr(st,  value,(Word16)1 );
   }
}


/************************************************************************
*
*	Function Name : Load_sh
*
*	Purpose :
*
*		Load the 16 bit var1 left shift(0-15) into the 32 bit output.
*		MS bits of the output are sign extended.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		shift
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0 <= var1 <= 15.
*
*	Outputs :
*
*		none.
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 Load_sh(tetra_codec* st, Word16 var1, Word16 shift)
{
	return( L_msu0(st,  (Word32)0, var1, POW2[shift]));
}


/************************************************************************
*
*	Function Name : Load_sh16
*
*	Purpose :
*
*		Load the 16 bit var1 with a left shift of 16 into the 32 bit output.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 Load_sh16(tetra_codec* st, Word16 var1)
{
	return( L_msu(st, (Word32)0, var1, (Word16)-32768));
}


/************************************************************************
*
*	Function Name : norm_v
*
*	Purpose :
*
*		Variable normalisation of a 32 bit integer (L_var3).
*		var1 gives the maximum number of left shift to be done
*		*var2 returns the actual number of left shift
*
*	Complexity Weight : 37
*
*	Inputs :
*
*		L_var3
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var3 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0 <= var1 <= 15.
*
*	Outputs :
*
*		*var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0 <= *var2 <= 15.
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 norm_v(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 *var2)
  {
   Word16 shift;

   shift = norm_l(L_var3);
   if(sub(st, shift, var1) > 0) shift = var1;
   *var2 = shift;
   return(L_shl(st, L_var3, shift));
  }


/************************************************************************
*
*	Function Name : store_hi
*
*	Purpose :
*
*		Store high part of a L_var1 with a left shift of var2.
*
*	Complexity Weight : 3
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0 <= var2 <= 7.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 store_hi(tetra_codec* st, Word32 L_var1, Word16 var2)
{
  static const Word16 SHR[8]={16, 15, 14, 13, 12, 11, 10, 9};
  return(extract_l( L_shr(st, L_var1, SHR[var2])));
}


/************************************************************************
*
*	Function Name : sub_sh
*
*	Purpose :
*
*		Subtract var1 with a left shift(0-15) to L_var2.
*		Control saturation and set overflow_flag.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var2
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var2 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		shift
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0 <= shift <= 15.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 sub_sh(tetra_codec* st, Word32 L_var2, Word16 var1, Word16 shift)
{
	return( L_mac0(st, L_var2, var1, POW2[shift]));
}


/************************************************************************
*
*	Function Name : sub_sh16
*
*	Purpose :
*
*		Subtract var1 with a left shift of 16 to L_var2.
*		Control saturation and set overflow flag.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var2
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var2 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.

************************************************************************/

Word32 sub_sh16(tetra_codec* st, Word32 L_var2, Word16 var1)
{
	return( L_mac(st, L_var2, var1, (Word16)-32768));
}





//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- FEXP_TET.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------





/************************************************************************
*
*	FILENAME		:	fexp_tet.c
*
*	DESCRIPTION		:	Library of special functions for operations in
*					extended precision used in the TETRA speech codec
*
************************************************************************
*
*	FUNCTIONS		:	- L_comp(st, )
*					- L_extract()
*					- mpy_mix(st, )
*					- mpy_32(st, )
*					- div_32()
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*
************************************************************************
*
*	COMMENTS		:	This subsection contains operations in double precision.
*					These operations are non standard double precision
*					operations.
*
*					They are used where single precision is not enough but the
*					full 32 bits precision is not necessary. For example, the
*					function div_32() has a 24 bits precision.
*
*					The double precision numbers use a special representation :
*
*					L_32 = hi<<15 + lo
*
*					L_32 is a 32 bit integer with b30 == b31.
*					hi and lo are 16 bit signed integers.
*					As the low part also contains the sign, this allows fast
*					multiplication.
*
*					0xc000 0000 <= L_32 <= 0x3fff ffff.
*
*					In general, DPF is used to specify this special
*					format.
*
************************************************************************/

//#include "source.h"


/************************************************************************
*
*	Function Name : L_comp
*
*	Purpose :
*
*		Compose from two 16 bit DPF a normal 32 bit integer.
*		L_32 = hi<<15 + lo
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		hi
*			msb
*
*		lo
*			lsb (with sign)
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_32
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0xc000 0000 <= L_32 <= 0x3fff ffff.
*
************************************************************************/

Word32 L_comp(tetra_codec* st, Word16 hi, Word16 lo)
{
  return(add_sh(st,  Load_sh(st,  lo,(Word16)0 ), hi, (Word16)15 ));
}


/************************************************************************
*
*	Function Name : L_extract
*
*	Purpose :
*
*		Extract from a 31 bit integer two 16 bit DPF.
*
*	Complexity Weight : 5
*
*	Inputs :
*
*		L_32
*			32 bit long signed integer (Word32) with b30 == b31
*			whose value falls in the range : 0xc000 0000 <= L_32 <= 0x3fff ffff.
*
*	Outputs :
*
*		hi
*			b15 to b30 of L_32
*
*		lo
*			L_32 - hi<<15
*
*	Returned Value :
*
*		none
*
************************************************************************/

void L_extract(tetra_codec* st, Word32 L_32, Word16 *hi, Word16 *lo)
{
  *hi  = extract_h( L_shl(st, L_32,(Word16)1 ) );
  *lo  = extract_l( sub_sh(st, L_32, *hi, (Word16)15 ) );
  return;
}


/************************************************************************
*
*	Function Name : mpy_mix
*
*	Purpose :
*
*		Multiply a 16 bit integer by a 32 bit (DPF).
*		The result is divided by 2**16
*		L_32 = hi1*lo2 + (lo1*lo2)>>15
*
*	Complexity Weight : 4
*
*	Inputs :
*
*		hi1
*			hi part of 32 bit number
*
*		lo1
*			lo part of 32 bit number
*
*		lo2
*			16 bit number
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 mpy_mix(tetra_codec* st, Word16 hi1, Word16 lo1, Word16 lo2)
{
  Word16 p1;
  Word32 L_32;

  p1   = extract_h(L_mult0(lo1, lo2));
  L_32 = L_mult0(hi1,lo2 );

  return(add_sh(st,  L_32, p1, (Word16)1 ));
}


/************************************************************************
*
*	Function Name : mpy_32
*
*	Purpose :
*
*		Multiply two 32 bit integers (DPF). The result is divided by 2**32
*		L_32 = hi1*hi2 + (hi1*lo2)>>15 + (lo1*hi2)>>15)
*
*	Complexity Weight : 7
*
*	Inputs :
*
*		hi1
*			hi part of first number
*
*		lo1
*			lo part of first number
*
*		hi2
*			hi part of second number
*
*		lo2
*			lo part of second number
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_var_out <= 0x7fff ffff.
*
*************************************************************************/

Word32 mpy_32(tetra_codec* st, Word16 hi1, Word16 lo1, Word16 hi2, Word16 lo2)
{
  Word16 p1, p2;
  Word32 L_32;

  p1   = extract_h(L_mult0(hi1, lo2));
  p2   = extract_h(L_mult0(lo1, hi2));
  L_32 = L_mult0(hi1, hi2);
  L_32 = add_sh(st,  L_32, p1, (Word16)1 );

  return(add_sh(st,  L_32, p2, (Word16)1 ));
}


/************************************************************************
*
*	Function Name : div_32
*
*	Purpose :
*
*		Fractionnal integer division of two 32 bit numbers.
*		L_num / L_denom
*		L_num and L_denom must be positive and L_num < L_denom
*		L_denom = denom_hi<<15 + denom_lo
*		denom_hi is a normalized number
*		The result is in Q30
*
*	Complexity Weight : 52
*
*	Inputs :
*
*		L_num
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_num <= L_denom.
*
*		(L_denom = denom_hi<<15 + denom_lo)
*
*		denom_hi
*			16 bit normalized integer whose value falls in the
*			range : 0x4000000 < hi < 0x7fff ffff.
*
*		denom_lo
*			16 bit positive integer whose value falls in the
*			range : 0 < lo < 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_div
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_div <= 0x3fff ffff.
*			L_div is a Q30 value (point between b30 and b29)
*
*	Algorithm :
*
*		 - find = 1/L_denom
*			First approximation: approx = 1 / denom_hi
*			1/L_denom = approx * (2.0 - L_denom * approx )
*		-  result = L_num * (1/L_denom)
*
************************************************************************/

Word32 div_32(tetra_codec* st, Word32 L_num, Word16 denom_hi, Word16 denom_lo)
{
  Word16 approx, hi, lo, n_hi, n_lo;
  Word32 t0;


  /* First approximation: 1 / L_denom = 1/denom_hi */

  approx = div_s(st,  (Word16)0x3fff, denom_hi);	/* result in Q15 */



/* 1/L_denom = approx * (2.0 - L_denom * approx) */

  t0 = mpy_mix(st, denom_hi, denom_lo, approx);	/* result in Q29 */

  t0 = L_sub(st,  (Word32)0x40000000, t0);		/* result in Q29 */

  L_extract(st, t0, &hi, &lo);

  t0 = mpy_mix(st, hi, lo, approx);			/* = 1/L_denom in Q28 */

  /* L_num * (1/L_denom) */

  L_extract(st, t0, &hi, &lo);
  L_extract(st, L_num, &n_hi, &n_lo);
  t0 = mpy_32(st, n_hi, n_lo, hi, lo);

  return( L_shl(st, t0,(Word16)2) );			/* From Q28 to Q30 */
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- FMAT_TET.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------





/************************************************************************
*
*	FILENAME		:	fmat_tet.c
*
*	DESCRIPTION		:	Library of mathematic functions used in the TETRA codec
*
************************************************************************
*
*	FUNCTIONS		:	- inv_sqrt()
*					- Log2()
*					- pow2()
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*
*					inv_sqrt.tab
*					log2.tab
*					pow2.tab
*
************************************************************************/

//#include "source.h"

/************************************************************************
*
*	FILENAME		:	inv_sqrt.tab
*
*	DESCRIPTION		:	Table for routine inv_sqrt().
*
************************************************************************/

static const Word16 tab_inv_sqrt[49] = {

 32767, 31790, 30894, 30070, 29309, 28602, 27945, 27330, 26755, 26214,
 25705, 25225, 24770, 24339, 23930, 23541, 23170, 22817, 22479, 22155,
 21845, 21548, 21263, 20988, 20724, 20470, 20225, 19988, 19760, 19539,
 19326, 19119, 18919, 18725, 18536, 18354, 18176, 18004, 17837, 17674,
 17515, 17361, 17211, 17064, 16921, 16782, 16646, 16514, 16384 };

//#include "inv_sqrt.tab"         /* Table for inv_sqrt() */

/************************************************************************
*
*	FILENAME		:	log2.tab
*
*	DESCRIPTION		:	Table for routine Log2().
*
************************************************************************/

static const Word16 tab_log2[33] = {
     0,  1455,  2866,  4236,  5568,  6863,  8124,  9352, 10549, 11716,
 12855, 13967, 15054, 16117, 17156, 18172, 19167, 20142, 21097, 22033,
 22951, 23852, 24735, 25603, 26455, 27291, 28113, 28922, 29716, 30497,
 31266, 32023, 32767 };

//#include "log2.tab"             /* Table for Log2() */

/************************************************************************
*
*	FILENAME		:	pow2.tab
*
*	DESCRIPTION		:	Table for routine pow2().
*
************************************************************************/

static const Word16 tab_pow2[33] = {
 16384, 16743, 17109, 17484, 17867, 18258, 18658, 19066, 19484, 19911,
 20347, 20792, 21247, 21713, 22188, 22674, 23170, 23678, 24196, 24726,
 25268, 25821, 26386, 26964, 27554, 28158, 28774, 29405, 30048, 30706,
 31379, 32066, 32767 };

//#include "pow2.tab"             /* Table for pow2() */


/************************************************************************
*
*	Function Name : inv_sqrt
*
*	Purpose :
*
*		Compute 1/sqrt(L_x).
*		L_x is positive. The result is in Q30.
*
*	Complexity Weight : 56
*
*	Inputs :
*
*		L_x
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_x <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_y
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_y <= 0x3fff ffff.
*			L_y is a Q30 value (point between b30 and b29)
*
*	Algorithm :
*
*		The function 1/sqrt(L_x) is approximated by a table (tab_inv_sqrt)
*		and linear interpolation :
*
*			1 - Normalization of L_x
*			2 - If (30-exponant) is even then shift right once
*			3 - exponant = (30-exponant)/2  +1
*			4 - i = bit25-b31 of L_x,    16 <= i <= 63  ->because of normalization
*			5 - a = bit10-b24
*			6 - i -=16
*			7 - L_y = tab_inv_sqrt[i]<<16 - (tab_inv_sqrt[i] - tab_inv_sqrt[i+1]) * a * 2
*			8 - L_y >>= exponant
*
************************************************************************/

Word32 inv_sqrt(tetra_codec* st, Word32 L_x)
{
  Word16 exp, i, a, tmp;
  Word32 L_y;

  if( L_x <= (Word32)0) return ( (Word32)0x3fffffff);


  exp = norm_l(L_x);
  L_x = L_shl(st, L_x, exp );		/* L_x is normalized */

  exp = sub(st,  (Word16)30, exp );
  if( (exp & 1) == 0 )			/* If exponant even -> shift right */
      L_x = L_shr(st, L_x, (Word16)1 );

  exp = shr(st,  exp, (Word16)1 );
  exp = add(st,  exp, (Word16)1 );

  L_x = L_shr(st, L_x, (Word16)9 );
  i   = extract_h(L_x);			/* Extract b25-b31 */
  L_x = L_shr(st, L_x, (Word16)1 );
  a   = extract_l(L_x);			/* Extract b10-b24 */
  a   = a & (Word16)0x7fff;

  i   = sub(st,  i, (Word16)16 );

  L_y = L_deposit_h(tab_inv_sqrt[i]);	/* tab_inv_sqrt[i] << 16 */
  tmp = sub(st, tab_inv_sqrt[i], tab_inv_sqrt[i+1]);
					    /* tab_inv_sqrt[i] - tab_inv_sqrt[i+1])*/
  L_y = L_msu(st, L_y, tmp, a);			/* L_y -=  tmp*a*2 */

  L_y = L_shr(st, L_y, exp);		/* denormalization */

  return(L_y);
}


/************************************************************************
*
*	Function Name : Log2
*
*	Purpose :
*
*		Compute Log2(L_x).
*		L_x is positive.
*
*	Complexity Weight : 48
*
*	Inputs :
*
*		L_x
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_x <= 0x7fff ffff.
*
*	Outputs :
*
*		exponant
*			Integer part of Log2()
*			16 bit  signed integer (Word16) whose value falls in the
*			range :   0 <= exponant <= 30
*
*		fraction
*			Fractional part of Log2()
*			16 bit signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= fraction <= 0x7fff.
*			It's a Q15 value (point between b15 and b16).
*
*	Returned Value :
*
*		none
*
*	Algorithm :
*
*		The function Log2(L_x) is approximated by a table (tab_log2)
*		and linear interpolation :
*
*			1 - Normalization of L_x
*			2 - exponant = 30-exponant
*			3 - i = bit25-b31 of L_x,    32 <= i <= 63  ->because of normalization
*			4 - a = bit10-b24
*			5 - i -=32
*			6 - fraction = tab_log2[i]<<16 - (tab_log2[i] - tab_log2[i+1]) * a * 2
*
************************************************************************/

void Log2(tetra_codec* st, Word32 L_x, Word16 *exponant, Word16 *fraction)
{
  Word16 exp, i, a, tmp;
  Word32 L_y;

  if( L_x <= (Word32)0 )
  {
    *exponant = 0;
    *fraction = 0;
    return;
  }

  exp = norm_l(L_x);
  L_x = L_shl(st, L_x, exp );			/* L_x is normalized */

  *exponant = sub(st,  (Word16)30, exp );

  L_x = L_shr(st, L_x, (Word16)9 );
  i   = extract_h(L_x);				/* Extract b25-b31 */
  L_x = L_shr(st, L_x, (Word16)1 );
  a   = extract_l(L_x);				/* Extract b10-b24 of fraction */
  a   = a & (Word16)0x7fff;

  i   = sub(st,  i, (Word16)32 );

  L_y = L_deposit_h(tab_log2[i]);		/* tab_log2[i] << 16 */
  tmp = sub(st, tab_log2[i], tab_log2[i+1]);	/* tab_log2[i] - tab_log2[i+1] */
  L_y = L_msu(st, L_y, tmp, a);			/* L_y -= tmp*a*2 */

  *fraction = extract_h( L_y);

  return;
}


/************************************************************************
*
*	Function Name : pow2
*
*	Purpose :
*
*		L_x = pow(2.0, exponant.fraction).
*
*	Complexity Weight : 17
*
*	Inputs :
*
*		exponant
*			Integer part
*			16 bit  signed integer (Word16) whose value falls in the
*			range :   0 <= exponant <= 30
*
*		fraction
*			Fractional part
*			16 bit signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= fraction <= 0x7fff.
*			It's a Q15 value (point between b15 and b16).
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_x
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_x <= 0x7fff ffff.
*
*	Algorithm :
*
*		The function pow2(L_x) is approximated by a table (tab_pow2)
*		and linear interpolation :
*
*			1 - i = bit11-b16 of fraction,   0 <= i <= 31
*			2 - a = bit0-b10  of fraction
*			3 - L_x = tab_pow2[i]<<16 - (tab_pow2[i] - tab_pow2[i+1]) * a * 2
*			4 - L_x = L_x >> (30-exponant)     (with rounding)
*
************************************************************************/

Word32 pow2(tetra_codec* st, Word16 exponant, Word16 fraction)
{
  Word16 exp, i, a, tmp;
  Word32 L_x;

  L_x = L_deposit_l(fraction);
  L_x = L_shl(st, L_x, (Word16)6 );
  i   = extract_h(L_x);				/* Extract b10-b16 of fraction */
  L_x = L_shr(st, L_x, (Word16)1 );
  a   = extract_l(L_x);				/* Extract b0-b9   of fraction */
  a   = a & (Word16)0x7fff;


  L_x = L_deposit_h(tab_pow2[i]);		/* tab_pow2[i] << 16 */
  tmp = sub(st, tab_pow2[i], tab_pow2[i+1]);	/* tab_pow2[i] - tab_pow2[i+1] */
  L_x = L_msu(st, L_x, tmp, a);			/* L_x -= tmp*a*2  */

  exp = sub(st,  (Word16)30, exponant );
  L_x = L_shr_r(st, L_x, exp);
  return(L_x);
}






//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- SUB_DSP.C -----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------





/************************************************************************
*
*	FILENAME		:	sub_dsp.c
*
*	DESCRIPTION		:	General Purpose Signal Processing Library
*
************************************************************************
*
*	SUB-ROUTINES	:	- Autocorr()
*					- Az_Lsp()
*					- Back_Fil()
*					- Chebps()
*					- Convolve()
*					- Fac_Pond()
*					- Get_Lsp_Pol()
*					- Int_Lpc4()
*					- Lag_Window()
*					- Levin_32()
*					- LPC_Gain()
*					- Lsp_Az()
*					- Pond_Ai()
*					- Residu()
*					- Syn_Filt()
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*					stdio.h
*					stdlib.h
*
*					grid .tab
*					lag_wind.tab
*					window.tab
*
************************************************************************/

//#include <stdio.h>
//#include <stdlib.h>
//#include "source.h"

/************************************************************************
*
*	FILENAME		:	grid.tab
*
*	DESCRIPTION		:	Table for routine Az_Lsp().
*
*					Vector grid[] is in Q15
*
*						grid[0] = 1.0
*						grid[TSC_grid_points+1] = -1.0
*						for (i = 1; i < TSC_grid_points; i++)
*						grid[i] = cos((6.283185307*i)/(2.0*TSC_grid_points))
*
************************************************************************/

#define TSC_grid_points 60

static const Word16 grid[TSC_grid_points+1] ={
  32760,     32723,     32588,     32364,     32051,     31651,
  31164,     30591,     29935,     29196,     28377,     27481,
  26509,     25465,     24351,     23170,     21926,     20621,
  19260,     17846,     16384,     14876,     13327,     11743,
  10125,      8480,      6812,      5126,      3425,      1714,
     -1,     -1715,     -3426,     -5127,     -6813,     -8481,
 -10126,    -11744,    -13328,    -14877,    -16385,    -17847,
 -19261,    -20622,    -21927,    -23171,    -24352,    -25466,
 -26510,    -27482,    -28378,    -29197,    -29936,    -30592,
 -31165,    -31652,    -32052,    -32365,    -32589,    -32724,
 -32760};

//#include "grid.tab"

/************************************************************************
*
*	FILENAME		:	lag_wind.tab
*
*	DESCRIPTION		:	Table of lag_window coefficents for the autocorrelation.
*
*						lag_wind[0] = 1.00000000 	(*)
*						lag_wind[1] = 0.99884027	(**)
*						lag_wind[2] = 0.99551868
*						lag_wind[3] = 0.99000722
*						lag_wind[4] = 0.98234236
*						lag_wind[5] = 0.97257471
*						lag_wind[6] = 0.96076828
*						lag_wind[7] = 0.94699991
*						lag_wind[8] = 0.93135828
*						lag_wind[9] = 0.91394323
*						lag_wind[10]= 0.89486438
*
*	COMMENTS		:	- (*) The first coefficient whose value = 1 is just
*					  mentioned for information, but not included in the code
*					- (**) All the other coefficents incorporate a scaling factor
*					  of 1,00005 corresponding to a noise floor of -43 dB
*					- This table uses a special extended precision format
*					  (see "comments" in the header of
*					   the file fexp_tet.c - Section A.5.2)
*
************************************************************************/

static const Word16 lag_h[10] = {
      32729,
      32621,
      32440,
      32189,
      31869,
      31482,
      31031,
      30518,
      29948,
      29322};

static const Word16 lag_l[10] = {
      32704,
       5120,
      18240,
      12928,
      10752,
      14912,
       9600,
      24512,
       3008,
      30016};

//#include "lag_wind.tab"

/************************************************************************
*
*	FILENAME		:	window.tab
*
*	DESCRIPTION		:	Asymetric Hamming window for LPC analysis.
*
************************************************************************/

#define TSC_L_window (Word16)256

static const Word16 window[TSC_L_window] = {
  2621,  2623,  2628,  2636,  2647,  2662,  2679,  2700,  2724,  2752,
  2782,  2816,  2852,  2892,  2936,  2982,  3031,  3084,  3140,  3199,
  3260,  3325,  3393,  3465,  3539,  3616,  3696,  3779,  3865,  3954,
  4047,  4141,  4239,  4340,  4444,  4550,  4659,  4771,  4886,  5003,
  5123,  5246,  5372,  5500,  5631,  5764,  5900,  6038,  6179,  6323,
  6468,  6617,  6767,  6920,  7075,  7233,  7392,  7554,  7718,  7884,
  8052,  8223,  8395,  8569,  8746,  8924,  9104,  9286,  9470,  9655,
  9842, 10031, 10221, 10413, 10607, 10802, 10999, 11197, 11396, 11597,
 11799, 12002, 12207, 12413, 12619, 12827, 13036, 13246, 13457, 13669,
 13882, 14095, 14309, 14524, 14740, 14956, 15173, 15391, 15608, 15827,
 16046, 16265, 16484, 16704, 16924, 17144, 17364, 17584, 17804, 18024,
 18245, 18465, 18684, 18904, 19124, 19343, 19561, 19780, 19998, 20215,
 20432, 20648, 20864, 21079, 21293, 21506, 21719, 21931, 22142, 22352,
 22561, 22769, 22976, 23181, 23386, 23589, 23791, 23992, 24191, 24389,
 24586, 24781, 24975, 25167, 25357, 25546, 25733, 25919, 26102, 26284,
 26464, 26642, 26819, 26993, 27165, 27336, 27504, 27670, 27834, 27996,
 28156, 28313, 28468, 28621, 28772, 28920, 29066, 29209, 29350, 29488,
 29624, 29757, 29888, 30016, 30142, 30265, 30385, 30502, 30617, 30729,
 30838, 30945, 31048, 31149, 31247, 31342, 31434, 31523, 31609, 31692,
 31772, 31850, 31924, 31995, 32063, 32128, 32190, 32249, 32304, 32357,
 32406, 32453, 32496, 32536, 32573, 32606, 32637, 32664, 32688, 32709,
 32727, 32741, 32753, 32761, 32765, 32767, 32767, 32718, 32572, 32329,
 31991, 31561, 31041, 30434, 29744, 28977, 28136, 27227, 26257, 25231,
 24156, 23039, 21888, 20709, 19511, 18301, 17088, 15878, 14680, 13501,
 12350, 11233, 10158,  9132,  8162,  7253,  6412,  5645,  4955,  4348,
  3828,  3397,  3059,  2817,  2670,  2621 };

//#include "window.tab"

/* TSC_pp = local LPC order, TSC_nc =TSC_pp/2 */
//#define TSC_pp (Word16)10
#define TSC_nc (Word16)5

/* Length for local impulse response */

#define  TSC_llg     (Word16)60


/**************************************************************************
*
*	ROUTINE				:	Autocorr
*
*	DESCRIPTION			:	Compute autocorrelations
*
**************************************************************************
*
*	USAGE				:	Autocorr(buffer_in,p,buffer_out1,buffer_out2)
*							(Routine_Name(input1,input2,output1,output2))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Input signal buffer
*						- Format : Word16
*
*	INPUT2			:	- Description : LPC order
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Autocorrelations buffer (msb)
*							- Format : Word16
*
*		OUTPUT2			:	- Description : Autocorrelations buffer (lsb)
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Autocorr(tetra_codec* st, Word16 x[], Word16 p, Word16 r_h[], Word16 r_l[])
{
  Word16 i, j, norm;
  Word16 y[TSC_L_window];
  Word32 sum;

  //extern Flag Overflow;

  /* Windowing of signal */

  for(i=0; i<TSC_L_window; i++)
    y[i] = mult_r(st, x[i], window[i]);

  /* Compute r[0] and test for overflow */

  do {
    st->Overflow = 0;
    sum = 1;				/* Avoid case of all zeros */
    for(i=0; i<TSC_L_window; i++)
      sum = L_mac0(st, sum, y[i], y[i]);

    /* If overflow divide y[] by 4 */

    if(st->Overflow != 0)
    {
      for(i=0; i<TSC_L_window; i++)
        y[i] = shr(st, y[i], (Word16)2);
    }

  }while(st->Overflow != 0);


 /* Normalization of r[0] */

  norm = norm_l(sum);
  sum  = L_shl(st, sum, norm);
  sum  = L_shr(st, sum, (Word16)1);		/* For special double format */
  L_extract(st, sum, &r_h[0], &r_l[0]);

 /* r[1] to r[p] */

  for (i = 1; i <= p; i++)
  {
    sum = 0;
    for(j=0; j<TSC_L_window-i; j++)
      sum = L_mac0(st, sum, y[j], y[j+i]);

    sum = L_shr(st, sum, (Word16)1);		/* Special double format */
    sum = L_shl(st, sum, norm);
    L_extract(st, sum, &r_h[i], &r_l[i]);
  }
  return;
}

/**************************************************************************
*
*	ROUTINE				:	Az_Lsp
*
*	DESCRIPTION			:	Compute the LSPs  in the cosine domain
*							from the LPC coefficients  (order=10)
*
**************************************************************************
*
*	USAGE				:	Az_Lsp(buffer_in1,buffer_out,buffer_in2)
*							(Routine_Name(input1,output1,input2))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Predictor coefficients
*						- Format : Word16 - Q12
*
*	INPUT2			:	- Description : Previous LSP values
*							          (in case not 10 roots are found)
*						- Format : Word16 - Q15
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Line spectral pairs in the
*								          cosine domain
*							- Format : Word16 - Q15
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Az_Lsp(tetra_codec* st, Word16 a[], Word16 lsp[], Word16 old_lsp[])
{
 Word16 i, j, nf, ip;
 Word16 xlow, ylow, xhigh, yhigh, xmid, ymid, xint;
 Word16 x, y, sign, exp;
 Word16 *coef;
 Word16 f1[TSC_pp/2+1], f2[TSC_pp/2+1];
 Word32 t0;

/*-------------------------------------------------------------*
 *  find the sum and diff. pol. F1(z) and F2(z)                *
 *    F1(z) <--- F1(z)/(1+z**-1) & F2(z) <--- F2(z)/(1-z**-1)  *
 *                                                             *
 * f1[0] = 1.0;                                                *
 * f2[0] = 1.0;                                                *
 *                                                             *
 * for (i = 0; i< TSC_nc; i++)                                     *
 * {                                                           *
 *   f1[i+1] = a[i+1] + a[TSC_pp-i] - f1[i] ;                       *
 *   f2[i+1] = a[i+1] - a[TSC_pp-i] + f2[i] ;                       *
 * }                                                           *
 *-------------------------------------------------------------*/

 f1[0] = 2048;				/* f1[0] = 1.0 is in Q11 */
 f2[0] = 2048;				/* f2[0] = 1.0 is in Q11 */


 for (i = 0; i< TSC_nc; i++)
 {
   t0 = Load_sh(st, a[i+1], (Word16)15);	/* a[i+1]  in Q27 */
   t0 = add_sh(st, t0, a[TSC_pp-i], (Word16)15);	/* +a[TSC_pp-i] in Q27 */
   t0 = sub_sh16(st, t0, f1[i]);			/* -f1[i]  in Q27 */
   f1[i+1] = extract_h(t0);		/* f1[i+1] = a[i+1] + a[TSC_pp-i] - f1[i] */
							/* result in Q11  */

   t0 = Load_sh(st, a[i+1], (Word16)15);	/* a[i+1]  in Q27   */
   t0 = sub_sh(st, t0, a[TSC_pp-i], (Word16)15);	/* -a[TSC_pp-i] in Q27 */
   t0 = add_sh16(st, t0, f2[i]);			/* +f2[i] in Q27  */
   f2[i+1] = extract_h(t0);		/* f2[i+1] = a[i+1] - a[TSC_pp-i] + f2[i] */
							/* result in Q11  */
 }

/*-------------------------------------------------------------*
 * find the LSPs using the Chebichev pol. evaluation           *
 *-------------------------------------------------------------*/

 nf=0;          /* number of found frequencies */
 ip=0;          /* indicator for f1 or f2      */

 coef = f1;

 xlow = grid[0];
 ylow = Chebps(st, xlow, coef, TSC_nc);

 j = 0;
 while ( (nf < TSC_pp) && (j < TSC_grid_points) )
 {
   j++;
   xhigh = xlow;
   yhigh = ylow;
   xlow  = grid[j];
   ylow  = Chebps(st, xlow,coef,TSC_nc);

   if ( L_mult0(ylow,yhigh) <= (Word32)0)
   {

     /* divide 4 times the interval */

     for (i = 0; i < 4; i++)
     {
       t0   = Load_sh(st, xlow, (Word16)15);	/* xmid = 0.5*(xlow + xhigh) */
       t0   = add_sh(st, t0, xhigh, (Word16)15);
       xmid = extract_h(t0);

       ymid = Chebps(st, xmid,coef,TSC_nc);

       if ( L_mult0(ylow,ymid) <= (Word32)0)
       {
         yhigh = ymid;
         xhigh = xmid;
       }
       else
       {
         ylow = ymid;
         xlow = xmid;
       }
     }


    /*-------------------------------------------------------------*
     * Linear interpolation                                        *
     *    xint = xlow - ylow*(xhigh-xlow)/(yhigh-ylow);            *
     *-------------------------------------------------------------*/

     x   = sub(st, xhigh, xlow);
     y   = sub(st, yhigh, ylow);

     if(y == 0)
     {
       xint = xlow;
     }
     else
     {
       sign= y;
       y   = abs_s(y);
       exp = norm_s(y);
       y   = shl(st, y, exp);
       y   = div_s(st,  (Word16)16383, y);
       t0  = L_mult0( x, y);
       t0  = L_shr(st, t0, sub(st, (Word16)19,exp) );
       y    = extract_l(t0);		/* y= (xhigh-xlow)/(yhigh-ylow) in Q10 */

       if(sign < 0) y = negate(y);

       t0   = Load_sh(st, xlow, (Word16)10);	/* xint = xlow - ylow*y */
       t0   = L_msu0(st, t0, ylow, y);
       xint = store_hi(st, t0, (Word16)6);

     }

     lsp[nf] = xint;
     xlow    = xint;
     nf++;

     if(ip == 0)
     {
       ip = 1;
       coef = f2;
     }
     else
     {
       ip = 0;
       coef = f1;
     }
     ylow = Chebps(st, xlow,coef,TSC_nc);

   }
 }

 /* Check if TSC_pp roots found */

 if( sub(st, nf, TSC_pp) < 0)
 {
    for(i=0; i<TSC_pp; i++)
       lsp[i] = old_lsp[i];
    //printf("\n !!Not 10 roots found in Az_Lsp()!!!n");
 }

 return;
}


/**************************************************************************
*
*	ROUTINE				:	Back_Fil
*
*	DESCRIPTION			:	Perform the Backward filtering of input vector
*						buffer_in1 with buffer_in2 and write the result
*						in output vector buffer_out.
*							All vectors are of length L
*
**************************************************************************
*
*	USAGE				:	Back_Fil(buffer_in1,buffer_in2,buffer_out,L)
*							(Routine_Name(input1,input2,output1,input3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Input vector
*						- Format : Word16
*
*	INPUT2			:	- Description : Impulse response
*						- Format : Word16 - Q12
*
*	INPUT3			:	- Description : Vector size
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Backward filtering result
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Back_Fil(tetra_codec* st, Word16 x[], Word16 h[], Word16 y[], Word16 L)
{
   Word16 i, j;
   Word32 s, max;
   Word32 y32[60];		/* Usually, dynamic allocation of L */

   /* first keep the result on 32 bits and find absolute maximum */

   max = 0;

   for (i = 0; i < L; i++)
   {
     s = 0;
     for (j = i; j <  L; j++)
       s = L_mac0(st, s, x[j], h[j-i]);

     y32[i] = s;

     s = L_abs(s);
     if(L_sub(st, s, max) > 0) max = s;
   }


   /* Find the number of right shifts to do on y32[]  */
   /* so that maximum is on 13 bits                   */

   j = norm_l(max);
   if( sub(st, j,(Word16)16) > 0) j = 16;

   j = sub(st, (Word16)18, j);

   for(i=0; i<L; i++)
     y[i] = extract_l( L_shr(st, y32[i], j) );

   return;
}

/**************************************************************************
*
*	ROUTINE				:	Chebps
*
*	DESCRIPTION			:	Evaluate the Chebishev polynomial series
*
**************************************************************************
*
*	USAGE				:	Chebps(x,buffer_in,n)
*							(Routine_Name(input1,input2,input3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Input value of evaluation
*							          x = cos(frequency)
*						- Format : Word16 - Q15
*
*	INPUT2			:	- Description : Coefficients of the polynomial series
*						- Format : Word16 - Q11
*
*	INPUT3			:	- Description : Polynomial order
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Value of the polynomial C(x)
*							- Format : Word16 - Q14
*							     (saturated to +/- 1.99)
*
*	COMMENTS			:	- The polynomial order is :
*						n = p/2   (p is the prediction order)
*
*						- The polynomial is given by :
*							C(x)=T_n(x)+f(1)T_n-1(x)+...+f(n-1)T_1(x)+ f(n)/2
*
**************************************************************************/

Word16 Chebps(tetra_codec* st, Word16 x, Word16 f[], Word16 n)
{
  Word16 i, cheb;
  Word16 b0_h, b0_l, b1_h, b1_l, b2_h, b2_l;
  Word32 t0;

 /* Note: All computation are done in Q24. */

  b2_h = 512;					/* b2 = 1.0 in Q24 DPF */
  b2_l = 0;

  t0 = Load_sh(st,  x,(Word16)10);		/* 2*x in Q24          */
  t0 = add_sh(st,  t0, f[1], (Word16)13);	/* + f[1] in Q24       */
  L_extract(st, t0, &b1_h, &b1_l);		/* b1 = 2*x + f[1]     */

  for (i = 2; i<n; i++)
  {
    t0 = mpy_mix(st, b1_h, b1_l, x);		/* t0 = x*b1                  */
    t0 = L_shl(st, t0,(Word16)1);			/* t0 = 2.0*x*b1              */
    t0 = sub_sh(st, t0, b2_l, (Word16)0);	/* t0 = 2.0*x*b1 - b2         */
    t0 = sub_sh(st, t0, b2_h, (Word16)15);
    t0 = add_sh(st, t0, f[i], (Word16)13);	/* t0 = 2.0*x*b1 - b2 + f[i]; */


    L_extract(st, t0, &b0_h, &b0_l);		/* b0 = 2.0*x*b1 - b2 + f[i]; */
    b2_l = b1_l;					/* b2 = b1; */
    b2_h = b1_h;
    b1_l = b0_l;					/* b1 = b0; */
    b1_h = b0_h;
  }
  t0 = mpy_mix(st, b1_h, b1_l, x);		/* t0    = x*b1;     */
  t0 = sub_sh(st, t0, b2_l, (Word16)0);		/* t0   -= b2;       */
  t0 = sub_sh(st, t0, b2_h, (Word16)15);
  t0 = add_sh(st, t0, f[n], (Word16)12);	/* t0   += 0.5*f[n]; */

  t0 = L_shl(st, t0, (Word16)6);			/* Q24 to Q30 with saturation */
  cheb = extract_h(t0);				/* Result in Q14           */

  return(cheb);
}

/**************************************************************************
*
*	ROUTINE				:	Convolve
*
*	DESCRIPTION			:	Perform the convolution between two input vectors
*						buffer_in1 and buffer_in2 and write the result in
*						output vector buffer_out.
*							All vectors are of length L
*
**************************************************************************
*
*	USAGE				:	Convolve(buffer_in1,buffer_in2,buffer_out,L)
*							(Routine_Name(input1,input2,output1,input3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Input vector
*						- Format : Word16
*
*	INPUT2			:	- Description : Impulse response
*						- Format : Word16 - Q12
*
*	INPUT3			:	- Description : Vector size
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Output vector
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Convolve(tetra_codec* st, Word16 x[], Word16 h[], Word16 y[], Word16 L)
{
   Word16 i, n;
   Word32 s;

   for (n = 0; n < L; n++)
   {
     s = 0;
     for (i = 0; i <= n; i++)
       s = L_mac0(st, s, x[i], h[n-i]);
     y[n] = store_hi(st, s, (Word16)4);		/* h is in Q12 */
   }

   return;
}


/**************************************************************************
*
*	ROUTINE				:	Fac_Pond
*
*	DESCRIPTION			:	Compute LPC spectral expansion factors (fac[])
*							with the LPC order fixed to 10
*
**************************************************************************
*
*	USAGE				:	Fac_Pond(gamma,buffer_out)
*							(Routine_Name(input1,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Spectral expansion
*						- Format : Word16 - Q15
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Computed expansion factors
*							- Format : Word16 - Q15
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	- fac[0] = gamma
*							- fac[i] = fac[i-1] * gamma	i=1,9
*
**************************************************************************/

void Fac_Pond(tetra_codec* st, Word16 gamma, Word16 fac[])
{
  Word16 i;

  fac[0] = gamma;
  for(i=1; i<TSC_pp; i++)
    fac[i] = tsc_round(st,  L_mult(st, fac[i-1], gamma) );

  return;
}


/**************************************************************************
*
*	ROUTINE				:	Get_Lsp_Pol
*
*	DESCRIPTION			:	Find the polynomial F1(z) or F2(z) from the LSPs
*
**************************************************************************
*
*	USAGE				:	Get_Lsp_Pol(buffer_in,buffer_out)
*							(Routine_Name(input1,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : line spectral pairs
*							          (cosine domaine)
*						- Format : Word16 - Q15
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Coefficients of F1 or F2
*							- Format : Word32 - Q24
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Get_Lsp_Pol(tetra_codec* st, Word16 *lsp, Word32 *f)
{
  Word16 i,j, hi, lo;
  Word32 t0;

   /* Computation in Q24 */

   *f = Load_sh(st, (Word16)4096,(Word16)12);	/* f[0] = 1.0;           in Q24  */
   f++;
   *f = 0;
   *f = sub_sh(st, *f, *lsp, (Word16)10);	/* f[1] = -2.0 * lsp[0]; in Q24  */
   f++;
   lsp += 2;					/* Advance lsp pointer           */

   for(i=2; i<=5; i++)
   {
     *f = f[-2];

     for(j=1; j<i; j++, f--)
     {
       L_extract(st, f[-1] ,&hi, &lo);
       t0 = mpy_mix(st, hi, lo, *lsp);		/* t0 = f[-1] * lsp    */
       t0 = L_shl(st, t0, (Word16)1);
       *f = L_add(st, *f, f[-2]);			/* *f += f[-2]         */
       *f = L_sub(st, *f, t0);			/* *f -= t0            */
     }
     *f   = sub_sh(st, *f, *lsp, (Word16)10);	/* *f -= lsp<<10       */
     f   += i;					/* Advance f pointer   */
     lsp += 2;					/* Advance lsp pointer */
   }

   return;
}


/**************************************************************************
*
*	ROUTINE				:	Int_Lpc4
*
*	DESCRIPTION			:	Perform the LPC interpolation for the 4 sub-frames.
*							The interpolation is done on the LSP computed
*							in the cosine domain
*
**************************************************************************
*
*	USAGE				:	Int_Lpc4(buffer_in1,buffer_in2,buffer_out)
*							(Routine_Name(input1,input2,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : LSP of previous frame
*						- Format : Word16
*
*	INPUT2			:	- Description : LSP of current frame
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : LPC coeff. vector for the 4 sub-frames
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Int_Lpc4(tetra_codec* st, Word16 lsp_old[], Word16 lsp_new[], Word16 a[])
{
  Word16 i, j, fac_new, fac_old;
  Word16 lsp[TSC_pp];
  Word32 t0;

  fac_new = 8192;       /* 1/4 in Q15 */
  fac_old = 24576;      /* 3/4 in Q15 */

  for(j=0; j<33; j+=11)
  {
    for(i=0; i<TSC_pp; i++)
    {
      t0 = L_mult(st, lsp_old[i], fac_old);
      t0 = L_mac(st, t0, lsp_new[i], fac_new);
      lsp[i] = extract_h(t0);
    }
    Lsp_Az(st, lsp, &a[j]);

    fac_old = sub(st, fac_old, (Word16)8192);
    fac_new = add(st, fac_new, (Word16)8192);
  }
  Lsp_Az(st, lsp_new, &a[33]);

  return;
}

/**************************************************************************
*
*	ROUTINE				:	Lag_Window
*
*	DESCRIPTION			:	Lag_window on autocorrelations
*
**************************************************************************
*
*	USAGE				:	Lag_Window(p,buffer1,buffer2)
*							(Routine_Name(input1,arg2,arg3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description :  LPC order
*						- Format : Word16
*
*	ARG2				:	- Description : Autocorrelations (msb)
*						- Format : Word16
*
*	ARG3				:	- Description : Autocorrelations (lsb)
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*	ARG2				:	- Description : Lag_Windowed autocorrelations (msb)
*						- Format : Word16
*
*	ARG3				:	- Description : Lag_Windowed autocorrelations (lsb)
*						- Format : Word16
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	r[i] *= lag_wind[i]
*
*						r[i] and lag_wind[i] are in special extended precision
*						(for details on this format, see "comments"
*						in the header of the file fexp_tet.c - Section A.5.2)
*
**************************************************************************/

void Lag_Window(tetra_codec* st, Word16 p, Word16 r_h[], Word16 r_l[])
{
  Word16 i;
  Word32 x;

  for(i=1; i<=p; i++)
  {
     x  = mpy_32(st, r_h[i], r_l[i], lag_h[i-1], lag_l[i-1]);
     L_extract(st, x, &r_h[i], &r_l[i]);
  }
  return;
}


/**************************************************************************
*
*	ROUTINE				:	Levin_32
*
*	DESCRIPTION			:	Computation of 10 LPC coefficients
*							based on the Levison-Durbin algorithm
*							in double precision
*
**************************************************************************
*
*	USAGE				:	Levin_32(buffer_in1,buffer_in2,buffer_out)
*							(Routine_Name(input1,input2,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Vector of autocorrelations (msb)
*						- Format : Word16 - 11 values
*
*	INPUT2			:	- Description : Vector of autocorrelations (lsb)
*						- Format : Word16 - 11 values
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : LPC coefficients
*							- Format : Word16 - Q12
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	Algorithm :
*
*							R[i]	autocorrelations
*							A[i]	filter coefficients
*							K 	reflection coefficients
*							Alpha	prediction gain.
*
*						Initialisations :
*
*							A[0] = 1
*							K    = -R[1]/R[0]
*							A[1] = K
*							Alpha = R[0] * (1-K**2]
*
*						DO for  i = 2 to TSC_pp
*
*							S =  SUM ( R[j]*A[i-j] ,j=1,i-1 ) +  R[i]
*							K = -S / Alpha
* 							An[j] = A[j] + K*A[i-j]   for j=1 to i-1
*								where   An[i] = new A[i]
*							An[i]=K
*							Alpha=Alpha * (1-K**2)
*
*						END
*
**************************************************************************


**************************************************************************
*
*	NOTES ON THE DYNAMICS OF THE COMPUTATIONS	:
*
*	The numbers used are in double precision with the following format :
*
*	A = AH <<15 + AL.  AH and AL are 16 bit signed integers. Since the LSB's also contain
*	a sign bit, this format does not correspond to standard 32 bit integers.  This format is used
*	since it allows fast execution of multiplications and divisions.
*
*	"DPF" will refer to this special format in the following text.
*	(for details on this format, see "comments" in the header of file fexp_tet.c - Section A.5.2)
*
*	The R[i] were normalized in routine Autocorr (hence, R[i] < 1.0).
*	The K[i] and Alpha are theoretically < 1.0.
*	The A[i], for a sampling frequency of 8 kHz, are in practice always inferior to 16.0.
*
*	These characteristics allow straigthforward fixed-point implementation.  The parameters are
*	represented as follows :
*
*	R[i]    Q30   +- .99..
*	K[i]    Q30   +- .99..
*	Alpha   Normalised -> mantissa in Q30 plus exponant
*	A[i]    Q26   +- 15.999..
*
*	The additions are performed in 32 bit.  For the summation used to compute the K[i],
*	numbers in Q30 are multiplied by numbers in Q26, with the results of the multiplications in
*	Q26, resulting in a dynamic of +/- 32.  This is sufficient to avoid overflow, since the final
*	result of the summation is necessarily < 1.0 as both the K[i] and Alpha are
*	theoretically < 1.0.
*
**************************************************************************/

/* Last A(z) for case of unstable filter */

//static const Word16 old_A[TSC_pp+1]={4096,0,0,0,0,0,0,0,0,0,0};

void Levin_32(tetra_codec* st, Word16 Rh[], Word16 Rl[], Word16 A[])
{
 Word16 i, j;
 Word16 hi, lo;
 Word16 Kh, Kl;              	/* reflexion coefficient; hi and lo          */
 Word16 alp_h, alp_l, alp_e;	/* Prediction gain; hi lo and exponant       */
 Word16 Ah[TSC_pp+1], Al[TSC_pp+1];	/* LPC coef. in double prec.                 */
 Word16 Anh[TSC_pp+1], Anl[TSC_pp+1];	/* LPC coef.for next iterat. in double prec. */
 Word32 t0, t1, t2;		/* temporary variable                        */


/* K = A[1] = -R[1] / R[0] */

  t1  = L_comp(st, Rh[1], Rl[1]);			/* R[1] in Q30      */
  t2  = L_abs(t1);				/* abs R[1]         */
  t0  = div_32(st, t2, Rh[0], Rl[0]);		/* R[1]/R[0] in Q30 */
  if(t1 > 0) t0= L_negate(t0);		/* -R[1]/R[0]       */
  L_extract(st, t0, &Kh, &Kl);			/* K in DPF         */
  t0 = L_shr(st, t0,(Word16)4);			/* A[1] in Q26      */
  L_extract(st, t0, &Ah[1], &Al[1]);		/* A[1] in DPF      */


/*  Alpha = R[0] * (1-K**2) */

  t0 = mpy_32(st, Kh ,Kl, Kh, Kl);		/* K*K      in Q30 */
  t0 = L_abs(t0);					/* Some case <0 !! */
  t0 = L_sub(st,  (Word32)0x3fffffff, t0 );	/* 1 - K*K  in Q30 */
  L_extract(st, t0, &hi, &lo);			/* DPF format      */
  t0 = mpy_32(st, Rh[0] ,Rl[0], hi, lo);	/* Alpha in Q30    */

/* Normalize Alpha */

  t0 = norm_v(st, t0, (Word16)12, &i);
  t0 = L_shr(st, t0, (Word16)1);
  L_extract(st, t0, &alp_h, &alp_l);		/* DPF format    */
  alp_e = i-1;					/* t0 was in Q30 */

/*--------------------------------------*
 * ITERATIONS  I=2 to TSC_pp                *
 *--------------------------------------*/

  for(i= 2; i<=TSC_pp; i++)
  {

    /* t0 = SUM ( R[j]*A[i-j] ,j=1,i-1 ) +  R[i] */

    t0 = 0;
    for(j=1; j<i; j++)
      t0 = L_add(st, t0, mpy_32(st, Rh[j], Rl[j], Ah[i-j], Al[i-j]));

    t0 = L_shl(st, t0,(Word16)4);			/* result in Q26->convert to Q30 */
							/* No overflow possible          */
    t1 = L_comp(st, Rh[i],Rl[i]);
    t0 = L_add(st, t0, t1);				/* add R[i] in Q30               */

    /* K = -t0 / Alpha */

    t1 = L_abs(t0);
    t2 = div_32(st, t1, alp_h, alp_l);		/* abs(t0)/Alpha                 */
    if(t0 > 0) t2= L_negate(t2);		/* K =-t0/Alpha                  */
    t2 = L_shl(st, t2, alp_e);			/* denormalize; compare to Alpha */
    L_extract(st, t2, &Kh, &Kl);			/* K in DPF                      */

    /* Test for unstable filter, if unstable keep old A(z) */

    if ( abs_s(Kh) > 32750)
    {
       for(j=0; j<=TSC_pp; j++) A[j] = st->old_A[j];
       return;
     }


     /*------------------------------------------*
     *  Compute new LPC coeff. -> An[i]         *
     *  An[j]= A[j] + K*A[i-j]     , j=1 to i-1 *
     *  An[i]= K                                *
     *------------------------------------------*/


    for(j=1; j<i; j++)
    {
      t0 = mpy_32(st, Kh, Kl, Ah[i-j], Al[i-j]);
      t0 = add_sh(st, t0, Ah[j], (Word16)15);
      t0 = add_sh(st, t0, Al[j], (Word16)0);
      L_extract(st, t0, &Anh[j], &Anl[j]);
    }
    t2 = L_shr(st, t2, (Word16)4);		/* t2 = K in Q30->convert to Q26 */
    L_extract(st, t2, &Anh[i], &Anl[i]);	/* An[i] in Q26                  */

    /*  Alpha = Alpha * (1-K**2) */

    t0 = mpy_32(st, Kh ,Kl, Kh, Kl);		/* K*K      in Q30 */
    t0 = L_abs(t0);				/* Some case <0 !! */
    t0 = L_sub(st,  (Word32)0x3fffffff, t0 );	/* 1 - K*K  in Q30 */
    L_extract(st, t0, &hi, &lo);			/* DPF format      */
    t0 = mpy_32(st, alp_h , alp_l, hi, lo);	/* Alpha in Q30    */

    /* Normalize Alpha */

    t0 = norm_v(st, t0, (Word16)12, &j);
    t0 = L_shr(st, t0, (Word16)1);
    L_extract(st, t0, &alp_h, &alp_l);		/* DPF format    */
    alp_e += j-1;					/* t0 was in Q30 */

    /* A[j] = An[j]   Note: can be done with pointers */

    for(j=1; j<=i; j++)
    {
      Ah[j] =Anh[j];
      Al[j] =Anl[j];
    }
  }

  /* Troncate A[i] in Q26 to Q12 with rounding */

  A[0] = 4096;
  for(i=1; i<=TSC_pp; i++)
  {
    t0   = L_comp(st, Ah[i], Al[i]);
    t0   = add_sh(st, t0, (Word16)1, (Word16)13);	/* rounding */
    st->old_A[i] = A[i] = store_hi(st, t0,(Word16)2);
  }

  return;
}


/**************************************************************************
*
*	ROUTINE				:	Lpc_Gain
*
*	DESCRIPTION			:	Compute energy of impulse response of 1/A(z)
*							on 60 points
*
**************************************************************************
*
*	USAGE				:	Lpc_Gain(buffer_in)
*							(Routine_Name(input1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : LPC coefficients
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Energy of impulse response of 1/A(z)
*							- Format : Word32 - Q20
*
**************************************************************************/

Word32 Lpc_Gain(tetra_codec* st, Word16 a[])
{
  Word16 i;
  Word32 ener;
  Word16 h[TSC_llg];

  /* Compute the impulse response of A(z) */

  h[0] = 1024;				/* 1.0 in Q10 */
  for(i=1; i<TSC_llg; i++) h[i]=0;
  Syn_Filt(st, a, h, h, TSC_llg, &h[1], (Word16)0);

  /* Compute the energy of the impulse response */

  ener = 0;
  for(i=0; i<TSC_llg; i++)
    ener = L_mac0(st, ener, h[i], h[i]);

  return(ener);
}


/**************************************************************************
*
*	ROUTINE				:	Lsp_Az
*
*	DESCRIPTION			:	Compute the LPC coefficients
*							from the LSPs in the cosine domain
*							(order=10)
*
**************************************************************************
*
*	USAGE				:	Lsp_Az(buffer_in,buffer_out)
*							(Routine_Name(input1,output1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Line spectral pairs in the
*								          cosine domain
*						- Format : Word16 - Q15
*
*	OUTPUT ARGUMENT(S)		:
*
*	OUTPUT1			:	- Description : Predictor coefficients
*							- Format : Word16 - Q12
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Lsp_Az(tetra_codec* st, Word16 lsp[], Word16 a[])
{
  Word16 i, j;
  Word32 f1[6], f2[6];
  Word32 t0;

  Get_Lsp_Pol(st, &lsp[0],f1);
  Get_Lsp_Pol(st, &lsp[1],f2);

  for (i = 5; i > 0; i--)
  {
    f1[i] = L_add(st, f1[i], f1[i-1]);			/* f1[i] += f1[i-1]; */
    f2[i] = L_sub(st, f2[i], f2[i-1]);			/* f2[i] -= f2[i-1]; */
  }

  a[0] = 4096;
  for (i = 1, j = 10; i <= 5; i++, j--)
  {
    t0   = L_add(st, f1[i], f2[i]);			/* f1[i] + f2[i] */
    a[i] = extract_l( L_shr_r(st, t0,(Word16)13) );	/*from Q24 to Q12 and * 0.5*/
    t0   = L_sub(st, f1[i], f2[i]);			/* f1[i] - f2[i] */
    a[j] = extract_l( L_shr_r(st, t0,(Word16)13) );	/*from Q24 to Q12 and * 0.5*/
  }

  return;
}


/**************************************************************************
*
*	ROUTINE				:	Pond_Ai
*
*	DESCRIPTION			:	Compute spectral expansion (a_exp[]) of LPC
*							coefficients (a[]) using spectral expansion
*							factors (fac[]) with the LPC order fixed to 10
*
**************************************************************************
*
*	USAGE				:	Pond_Ai(buffer_in1,buffer_in2,buffer_out)
*							(Routine_Name(input1,input2,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : LPC coefficients
*						- Format : Word16
*
*	INPUT2			:	- Description : Spectral expansion factors
*						- Format : Word16 - Q15
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Spectral expanded LPC coefficients
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	- a_exp[i] = a[i] * fac[i-1]	i=1,10
*
**************************************************************************/

void Pond_Ai(tetra_codec* st, Word16 a[], Word16 fac[], Word16 a_exp[])
{
  Word16 i;

  a_exp[0] = a[0];
  for(i=1; i<=TSC_pp; i++)
    a_exp[i] = tsc_round(st,  L_mult(st, a[i], fac[i-1]) );

  return;
}


/**************************************************************************
*
*	ROUTINE				:	Residu
*
*	DESCRIPTION			:	Compute the LPC residual  by filtering the input
*							speech through A(z)
*
**************************************************************************
*
*	USAGE				:	Residu(buffer_in1,buffer_in2,buffer_out,lg)
*							(Routine_Name(input1,input2,output1,input3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Prediction coefficients
*						- Format : Word16 - Q12
*
*	INPUT2			:	- Description : Input speech - values of buffer_in2[]
*							          from -p to -1 are needed
*						- Format : Word16
*
*	INPUT3			:	- Description : Size of filtering
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Residual signal
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Residu(tetra_codec* st, Word16 a[], Word16 x[], Word16 y[], Word16 lg)
{
  Word16 i, j;
  Word32 s;

  for (i = 0; i < lg; i++)
  {
    s = Load_sh(st, x[i], (Word16)12);
    for (j = 1; j <= TSC_pp; j++)
      s = L_mac0(st, s, a[j], x[i-j]);

    s = add_sh(st, s, (Word16)1, (Word16)11);		/* Rounding */
    s = L_shl(st, s, (Word16)4);				/* Saturation */
    y[i] = extract_h(s);
  }
  return;
}


/**************************************************************************
*
*	ROUTINE				:	Syn_Filt
*
*	DESCRIPTION			:	Perform the synthesis filtering 1/A(z)
*
**************************************************************************
*
*	USAGE				:	Syn_Filt(buffer_in1,buffer_in2,buffer_out,
*							lg,buffer,flag)
*							(Routine_Name(input1,input2,output1,
*							input3,arg4,input5))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Prediction coefficients
*						- Format : Word16 - Q12
*
*	INPUT2			:	- Description : Input signal
*						- Format : Word16
*
*	INPUT3			:	- Description : Size of filtering
*						- Format : Word16
*
*	ARG4				:	- Description : Memory associated with this filtering
*						- Format : Word16
*
*		INPUT5			:	- Description :	- flag = 0  -> no update of memory
*									- flag = 1  -> update of memory
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Output signal
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Syn_Filt(tetra_codec* st, Word16 a[], Word16 x[], Word16 y[], Word16 lg, Word16 mem[],
              Word16 update)
{
  Word16 i, j;
  Word32 s;
  Word16 tmp[80];	/* This is usually done by memory allocation (lg+TSC_pp) */
  Word16 *yy;

  /* Copy mem[] to yy[] */

  yy = tmp;

  for(i=0; i<TSC_pp; i++)
    *yy++ = mem[i];


  /* Do the filtering. */

  for (i = 0; i < lg; i++)
  {
    s = Load_sh(st, x[i], (Word16)12);				/* a[] is in Q12 */
    for (j = 1; j <= TSC_pp; j++)
      s = L_msu0(st, s, a[j], yy[-j]);

    s     = add_sh(st, s, (Word16)1, (Word16)11);		/* Rounding */
    *yy++ = extract_h( L_shl(st, s, (Word16)4) );
  }

  for(i=0; i<lg; i++) y[i] = tmp[i+TSC_pp];

  /* Update of memory if update==1 */

  if(update != 0)
     for (i = 0; i < TSC_pp; i++) mem[i] = y[lg-TSC_pp+i];

 return;
}





//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- SCOD_TET.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


/************************************************************************
*
*	FILENAME		:	scod_tet.c
*
*	DESCRIPTION		:	Main routines for speech source encoding
*
************************************************************************
*
*	SUB-ROUTINES	:	- Init_Coder_Tetra()
*					- Coder_Tetra()
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*
************************************************************************/

//#include "source.h"

/*----------------------------------------------------------------------*
 *         Coder constant parameters.                                   *
 *                                                                      *
 *   TSC_L_window    : LPC analysis window size.                            *
 *   TSC_L_next      : Samples of next frame needed for autocor.            *
 *   TSC_L_frame     : Frame size.                                          *
 *   TSC_L_subfr     : Sub-frame size.                                      *
 *   TSC_p           : LPC order.                                           *
 *   TSC_pp1         : LPC order+1                                          *
 *   TSC_L_total     : Total speech size.                                   *
 *   TSC_dim_rr      : Dimension of matrix rr[][].                          *
 *   TSC_pit_min     : Minimum pitch lag.                                   *
 *   TSC_pit_max     : Maximum pitch lag.                                   *
 *   TSC_L_inter     : Length of filter for interpolation                   *
 *----------------------------------------------------------------------*/

//#define  TSC_L_window (Word16)256
//#define  TSC_L_next   (Word16)40
//#define  TSC_L_frame  (Word16)240
//#define  TSC_L_subfr  (Word16)60
//#define  TSC_p        (Word16)10
//#define  TSC_pp1      (Word16)11
//#define  TSC_L_total  (Word16)(TSC_L_frame+TSC_L_next+TSC_p)
//#define  TSC_dim_rr   (Word16)32
#define  TSC_pit_min  (Word16)20
//#define  TSC_pit_max  (Word16)143
//#define  TSC_L_inter  (Word16)15

/*--------------------------------------------------------*
 *   LPC bandwidth expansion factors.                     *
 *      In Q15 = 0.95, 0.60, 0.75, 0.85                   *
 *--------------------------------------------------------*/

#define TSC_gamma1  (Word16)31130
#define TSC_gamma2  (Word16)19661
#define TSC_gamma3  (Word16)24576
#define TSC_gamma4  (Word16)27853


/*--------------------------------------------------------*
 *         Static memory allocation.                      *
 *--------------------------------------------------------*/

        /* Speech vector */

//static Word16 old_speech[TSC_L_total];
//static Word16 *speech, *p_window;
//Word16 *new_speech;                    /* Global variable */

        /* Weighted speech vector */

//static Word16 old_wsp[TSC_L_frame+TSC_pit_max];
//static Word16 *wsp;

        /* Excitation vector */

//static Word16 old_exc[TSC_L_frame+TSC_pit_max+TSC_L_inter];
//static Word16 *exc;

        /* All-zero vector */

//static Word16 ai_zero[TSC_L_subfr+TSC_pp1];
//static Word16 *zero;

        /* Spectral expansion factors */

//static Word16 F_gamma1[TSC_p];
//static Word16 F_gamma2[TSC_p];
//static Word16 F_gamma3[TSC_p];
//static Word16 F_gamma4[TSC_p];

        /* Lsp (Line Spectral Pairs in the cosine domain) */

//static Word16 lspold[TSC_p]={
//              30000, 26000, 21000, 15000, 8000, 0,
//		  -8000,-15000,-21000,-26000};
//static Word16 lspnew[TSC_p];
//static Word16 lspnew_q[TSC_p], lspold_q[TSC_p];

	  /* Initial lsp values used after each time */
        /* a reset is executed */

static const Word16 lspold_init[TSC_p]={
              30000, 26000, 21000, 15000, 8000, 0,
		  -8000,-15000,-21000,-26000};

        /* Filters memories */

//static Word16 mem_syn[TSC_p], mem_w0[TSC_p], mem_w[TSC_p];

        /* Matrix rr[TSC_dim_rr][TSC_dim_rr] */

//static Word16 rr[TSC_dim_rr][TSC_dim_rr];

       /* Global definition */

//Word16 last_ener_cod;
//Word16 last_ener_pit;


/**************************************************************************
*
*	ROUTINE				:	Init_Coder_Tetra
*
*	DESCRIPTION			:	Initialization of variables for the speech encoder
*
**************************************************************************
*
*	USAGE				:	Init_Coder_Tetra()
*
*	INPUT ARGUMENT(S)		:	None
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Init_Coder_Tetra(tetra_codec* st)
{
  Word16 i,j;


/*-----------------------------------------------------------------------*
 *      Initialize pointers to speech vector.                            *
 *                                                                       *
 *                                                                       *
 *   |--------------------|--------|--------|......|-------|-------|     *
 *     previous speech       sf1      sf2             sf4    TSC_L_next      *
 *                                                                       *
 *   <----------------  Total speech vector (TSC_L_total)   ----------->     *
 *   |           <----  LPC analysis window (TSC_L_window)  ----------->     *
 *   |           |        <---- present frame (TSC_L_frame) --->             *
 *  old_speech   |        |       <----- new speech (TSC_L_frame) ----->     *
 *            p_window    |       |                                      *
 *                     speech     |                                      *
 *                             new_speech                                *
 *-----------------------------------------------------------------------*/

  st->new_speech = st->old_speech + TSC_L_total - TSC_L_frame;	/* New speech     */
  st->speech     = st->new_speech - TSC_L_next;			/* Present frame  */
  st->p_window   = st->old_speech + TSC_L_total - TSC_L_window;	/* For LPC window */

  /* Initialize global variables */

  st->last_ener_cod = 0;
  st->last_ener_pit = 0;

  /* Initialize static pointers */

  st->wsp    = st->old_wsp + TSC_pit_max;
  st->exc    = st->old_exc + TSC_pit_max + TSC_L_inter;
  st->zero   = st->ai_zero + TSC_pp1;


  /* Static vectors to zero */

  for(i=0; i<TSC_L_total; i++)
    st->old_speech[i] = 0;

  for(i=0; i<TSC_pit_max + TSC_L_inter; i++)
    st->old_exc[i] = st->old_wsp[i] = 0;

  for(i=0; i<TSC_p; i++)
    st->mem_syn[i] = st->mem_w[i] = st->mem_w0[i] = 0;

  for(i=0; i<TSC_L_subfr; i++)
    st->zero[i] = 0;

  for(i=0; i< TSC_dim_rr; i++)
    for(j=0; j< TSC_dim_rr; j++)
    st->rr[i][j] = 0;

  /* Initialisation of lsp values for first */
  /* frame lsp interpolation */

  for(i=0; i<TSC_p; i++)
    st->lspold_q[i] = st->lspold[i] = lspold_init[i];

  /* Compute LPC spectral expansion factors */

  Fac_Pond(st, TSC_gamma1, st->F_gamma1);
  Fac_Pond(st, TSC_gamma2, st->F_gamma2);
  Fac_Pond(st, TSC_gamma3, st->F_gamma3);
  Fac_Pond(st, TSC_gamma4, st->F_gamma4);

 return;
}


/**************************************************************************
*
*	ROUTINE				:	Coder_Tetra
*
*	DESCRIPTION			:	Main speech coder function
*
**************************************************************************
*
*	USAGE				:	Coder_Tetra(ana,synth)
*							(Routine_Name(output1,output2))
*
*	INPUT ARGUMENT(S)		:	None
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Analysis parameters
*							- Format : 23 * 16 bit-samples
*
*		OUTPUT2			:	- Description : Local synthesis
*							- Format : 240 * 16 bit-samples
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	- 240 speech data should have been copied to vector
*						new_speech[].  This vector is global and is declared in
*						this function.
*						- Output2 is for debugging only
*
**************************************************************************/

void Coder_Tetra(tetra_codec* st, Word16 ana[], Word16 synth[])
{

  /* LPC coefficients */

  Word16 r_l[TSC_pp1], r_h[TSC_pp1];	/* Autocorrelations low and high        */
  Word16 A_t[(TSC_pp1)*4];		/* A(z) unquantized for the 4 subframes */
  Word16 Aq_t[(TSC_pp1)*4];		/* A(z)   quantized for the 4 subframes */
  Word16 Ap1[TSC_pp1];		/* A(z) with spectral expansion         */
  Word16 Ap2[TSC_pp1];		/* A(z) with spectral expansion         */
  Word16 Ap3[TSC_pp1];		/* A(z) with spectral expansion         */
  Word16 Ap4[TSC_pp1];		/* A(z) with spectral expansion         */
  Word16 *A, *Aq;			/* Pointer on A_t and Aq_t              */

  /* Other vectors */

  Word16 h1[TSC_L_subfr];
  Word16 zero_h2[TSC_L_subfr+64], *h2;
  Word16 zero_F[TSC_L_subfr+64],  *F;
  Word16 res[TSC_L_subfr];
  Word16 xn[TSC_L_subfr];
  Word16 xn2[TSC_L_subfr];
  Word16 dn[TSC_L_subfr+4];
  Word16 code[TSC_L_subfr+4];
  Word16 y1[TSC_L_subfr];
  Word16 y2[TSC_L_subfr];



  /* Scalars */

  Word16 i, i_subfr;
  Word16 T0, T0_min, T0_max, T0_frac;
  Word16 gain_pit, gain_code, index;
  Word16 sign_code, shift_code;
  Word16 temp;
  Word32 L_temp;

  /* Initialization of F and h2 */

  F  = &zero_F[64];
  h2 = &zero_h2[64];
  for(i=0; i<64; i++)
   zero_F[i] = zero_h2[i] = 0;

/*------------------------------------------------------------------------*
 *  - Perform LPC analysis:                                               *
 *       * autocorrelation + lag windowing                                *
 *       * Levinson-Durbin algorithm to find a[]                          *
 *       * convert a[] to lsp[]                                           *
 *       * quantize and code the LSPs                                     *
 *       * find the interpolated LSPs and convert to a[] for all          *
 *         subframes (both quantized and unquantized)                     *
 *------------------------------------------------------------------------*/

  Autocorr(st, st->p_window, TSC_p, r_h, r_l);		/* Autocorrelations */

  Lag_Window(st, TSC_p, r_h, r_l);			/* Lag windowing    */

  Levin_32(st, r_h, r_l, A_t);			/* Levinson-Durbin  */

  Az_Lsp(st, A_t, st->lspnew, st->lspold);		/* From A(z) to lsp */

  Clsp_334(st, st->lspnew, st->lspnew_q, ana);		/* Lsp quantization */

  ana += 3;				/* Increment analysis parameters pointer */

  /* Interpolation of LPC for the 4 subframes */

  Int_Lpc4(st, st->lspold,   st->lspnew,   A_t);
  Int_Lpc4(st, st->lspold_q, st->lspnew_q, Aq_t);

  /* update the LSPs for the next frame */

  for(i=0; i<TSC_p; i++)
  {
    st->lspold[i]   = st->lspnew[i];
    st->lspold_q[i] = st->lspnew_q[i];
  }


 /*----------------------------------------------------------------------*
  * - Find the weighted input speech wsp[] for the whole speech frame    *
  * - Find open-loop pitch delay                                         *
  * - Set the range for searching closed-loop pitch                      *
  *----------------------------------------------------------------------*/

  A = A_t;
  for (i = 0; i < TSC_L_frame; i += TSC_L_subfr)
  {
    Pond_Ai(st, A, st->F_gamma1, Ap1);
    Pond_Ai(st, A, st->F_gamma2, Ap2);
    Residu(st, Ap1, &st->speech[i], &st->wsp[i], TSC_L_subfr);
    Syn_Filt(st, Ap2, &st->wsp[i], &st->wsp[i], TSC_L_subfr, st->mem_w, (Word16)1);
    A += TSC_pp1;
  }

  /* Find open loop pitch delay */

  T0 = Pitch_Ol_Dec(st, st->wsp, TSC_L_frame);

  /* range for closed loop pitch search */

  T0_min = sub(st, T0, (Word16)2);
  if (T0_min < TSC_pit_min) T0_min = TSC_pit_min;
  T0_max = add(st, T0_min, (Word16)4);
  if (T0_max > TSC_pit_max)
  {
     T0_max = TSC_pit_max;
     T0_min = sub(st, T0_max, (Word16)4);
  }


 /*------------------------------------------------------------------------*
  *          Loop for every subframe in the analysis frame                 *
  *------------------------------------------------------------------------*
  *  To find the pitch and innovation parameters. The subframe size is     *
  *  TSC_L_subfr and the loop is repeated TSC_L_frame/TSC_L_subfr times.               *
  *     - find the weighted LPC coefficients                               *
  *     - find the LPC residual signal res[]                               *
  *     - compute the target signal for pitch search                       *
  *     - compute impulse response of weighted synthesis filter (h1[])     *
  *     - find the closed-loop pitch parameters                            *
  *     - encode the pitch delay                                           *
  *     - update the impulse response h1[] by including fixed-gain pitch   *
  *     - find the autocorrelations of h1[] (rr[][])                       *
  *     - find target vector for codebook search                           *
  *     - backward filtering of target vector                              *
  *     - codebook search                                                  *
  *     - encode codebook address                                          *
  *     - VQ of pitch and codebook gains                                   *
  *     - find synthesis speech                                            *
  *     - update states of weighting filter                                *
  *------------------------------------------------------------------------*/

  Aq = Aq_t;	/* pointer to interpolated quantized LPC parameters */

  for (i_subfr = 0;  i_subfr < TSC_L_frame; i_subfr += TSC_L_subfr)
  {


   /*---------------------------------------------------------------*
     * Find the weighted LPC coefficients for the weighting filter.  *
     *---------------------------------------------------------------*/

    Pond_Ai(st, Aq, st->F_gamma3, Ap3);
    Pond_Ai(st, Aq, st->F_gamma4, Ap4);


    /*---------------------------------------------------------------*
     * Compute impulse response, h1[], of weighted synthesis filter  *
     *---------------------------------------------------------------*/

    st->ai_zero[0] = 4096;				/* 1 in Q12 */
    for (i = 1; i <= TSC_p; i++) st->ai_zero[i] = 0;

    Syn_Filt(st, Ap4, st->ai_zero, h1, TSC_L_subfr, st->zero, (Word16)0);

    /*---------------------------------------------------------------*
     * Compute LPC residual and copy it to exc[i_subfr]              *
     *---------------------------------------------------------------*/

    Residu(st, Aq, &st->speech[i_subfr], res, TSC_L_subfr);

    for(i=0; i<TSC_L_subfr; i++) st->exc[i_subfr+i] = res[i];

    /*---------------------------------------------------------------*
     * Find the target vector for pitch search:  ->xn[]              *
     *---------------------------------------------------------------*/

    Syn_Filt(st, Ap4, res, xn, TSC_L_subfr, st->mem_w0, (Word16)0);

    /*----------------------------------------------------------------------*
     *                 Closed-loop fractional pitch search                  *
     *----------------------------------------------------------------------*
     * The pitch range for the first subframe is divided as follows:        *
     *   19 1/3  to   84 2/3   resolution 1/3                               *
     *   85      to   143      resolution 1                                 *
     *                                                                      *
     * The period in the first subframe is encoded with 8 bits.             *
     * For the range with fractions:                                        *
     *   code = (T0-19)*3 + frac - 1;   where T0=[19..85] and frac=[-1,0,1] *
     * and for the integer only range                                       *
     *   code = (T0 - 85) + 197;        where T0=[86..143]                  *
     *----------------------------------------------------------------------*
     * For other subframes: if t0 is the period in the first subframe then  *
     * T0_min=t0-5   and  T0_max=T0_min+9   and  the range is given by      *
     *      T0_min-1 + 1/3   to  T0_max + 2/3                               *
     *                                                                      *
     * The period in the 2nd,3rd,4th subframe is encoded with 5 bits:       *
     *  code = (T0-(T0_min-1))*3 + frac - 1;  where T0[T0_min-1 .. T0_max+1]*
     *---------------------------------------------------------------------*/

    T0 = Pitch_Fr(st, &st->exc[i_subfr], xn, h1, TSC_L_subfr, T0_min, T0_max,
                  i_subfr, &T0_frac);

    if (i_subfr == 0)
    {
      /* encode pitch delay (with fraction) */

      if (T0 <= 85)
      {
        /* index = T0*3 - 58 + T0_frac; */
        index = add(st, T0, add(st, T0, T0));
        index = sub(st, index, (Word16)58);
        index = add(st, index, T0_frac);
      }
      else
        index = add(st, T0, (Word16)112);


      /* find T0_min and T0_max for other subframes */

      T0_min = sub(st, T0, (Word16)5);
      if (T0_min < TSC_pit_min) T0_min = TSC_pit_min;
      T0_max = add(st, T0_min, (Word16)9);
      if (T0_max > TSC_pit_max)
      {
        T0_max = TSC_pit_max;
        T0_min = sub(st, T0_max, (Word16)9);
      }
    }

    else						/* other subframes */
    {
      i = sub(st, T0, T0_min);
							/* index = i*3 + 2 + T0_frac;  */
      index = add(st, i, add(st, i, i));
      index = add(st, index, (Word16)2);
      index = add(st, index, T0_frac);
    }

    *ana++ = index;


   /*-----------------------------------------------------------------*
    *   - find unity gain pitch excitation (adaptive codebook entry)  *
    *     with fractional interpolation.                              *
    *   - find filtered pitch exc. y1[]=exc[] filtered by 1/Ap4(z)    *
    *   - compute pitch gain and limit between 0 and 1.2              *
    *   - update target vector for codebook search                    *
    *-----------------------------------------------------------------*/


    Pred_Lt(st, &st->exc[i_subfr], T0, T0_frac, TSC_L_subfr);

    Syn_Filt(st, Ap4, &st->exc[i_subfr], y1, TSC_L_subfr, st->zero, (Word16)0);

    gain_pit = G_Pitch(st, xn, y1, TSC_L_subfr);

	/* xn2[i] = xn[i] - y1[i]*gain_pit */

    for (i = 0; i < TSC_L_subfr; i++)
    {
      L_temp = L_mult(st, y1[i], gain_pit);
      L_temp = L_shl(st, L_temp, (Word16)3);	/* gain_pit in Q12 */
      L_temp = L_sub(st,  Load_sh16(st, xn[i]), L_temp);
      xn2[i] = extract_h(L_temp);
    }


   /*----------------------------------------------------------------*
    * -Compute impulse response F[] and h2[] for innovation codebook *
    * -Find correlations of h2[];  rr[i][j] = sum h2[n-i]*h2[n-j]    *
    *----------------------------------------------------------------*/

    for (i = 0; i <= TSC_p; i++) st->ai_zero[i] = Ap3[i];
    Syn_Filt(st, Ap4, st->ai_zero, F, TSC_L_subfr, st->zero, (Word16)0);

    /* Introduce pitch contribution with fixe gain of 0.8 to F[] */

    for (i = T0; i < TSC_L_subfr; i++)
    {
      temp = mult(st, F[i-T0], (Word16)26216);
      F[i] = add(st, F[i], temp);
    }

    /* Compute h2[]; -> F[] filtered by 1/Ap4(z) */

    Syn_Filt(st, Ap4, F, h2, TSC_L_subfr, st->zero, (Word16)0);

    Cal_Rr2(st, h2, (Word16*)st->rr);

   /*-----------------------------------------------------------------*
    * - Backward filtering of target vector (find dn[] from xn2[])    *
    * - Innovative codebook search (find index and gain)              *
    *-----------------------------------------------------------------*/

    Back_Fil(st, xn2, h2, dn, TSC_L_subfr);	/* backward filtered target vector dn */

    *ana++ =D4i60_16(st, dn,F,h2, st->rr, code, y2, &sign_code, &shift_code);
    *ana++ = sign_code;
    *ana++ = shift_code;
    gain_code = G_Code(st, xn2, y2, TSC_L_subfr);

   /*-----------------------------------------------------------------*
    * - Quantization of gains.                                        *
    *-----------------------------------------------------------------*/

    *ana++ = Ener_Qua(st, Aq,&st->exc[i_subfr],code, TSC_L_subfr, &gain_pit, &gain_code);

   /*-------------------------------------------------------*
    * - Find the total excitation                           *
    * - Update filter memory mem_w0 for finding the target  *
    *   vector in the next subframe.                        *
    *   The filter state mem_w0[] is found by filtering by  *
    *   1/Ap4(z) the error between res[i] and exc[i]        *
    *-------------------------------------------------------*/

    for (i = 0; i < TSC_L_subfr;  i++)
    {
      /* exc[i] = gain_pit*exc[i] + gain_code*code[i]; */
      /* exc[i]  in Q0   gain_pit in Q12               */
      /* code[i] in Q12  gain_cod in Q0                */
      L_temp = L_mult0(st->exc[i+i_subfr], gain_pit);
      L_temp = L_mac0(st, L_temp, code[i], gain_code);
      st->exc[i+i_subfr] = (Word16)L_shr_r(st, L_temp, (Word16)12);
    }

    for(i=0; i<TSC_L_subfr; i++)
      res[i] = sub(st, res[i], st->exc[i_subfr+i]);

    Syn_Filt(st, Ap4, res, code, TSC_L_subfr, st->mem_w0, (Word16)1);

   /* Note: we use vector code[] as output only as temporary vector */

   /*-------------------------------------------------------*
    * - find synthesis speech corresponding to exc[]        *
    *   This filter is to help debug only.                  *
    *-------------------------------------------------------*/

    Syn_Filt(st, Aq, &st->exc[i_subfr], &synth[i_subfr], TSC_L_subfr, st->mem_syn,
			(Word16)1);

    Aq += TSC_pp1;
  }

 /*--------------------------------------------------*
  * Update signal for next frame.                    *
  * -> shift to the left by TSC_L_frame:                 *
  *     speech[], wsp[] and  exc[]                   *
  *--------------------------------------------------*/

  for(i=0; i< TSC_L_total-TSC_L_frame; i++)
    st->old_speech[i] = st->old_speech[i+TSC_L_frame];

  for(i=0; i<TSC_pit_max; i++)
    st->old_wsp[i] = st->old_wsp[i+TSC_L_frame];

  for(i=0; i<TSC_pit_max+TSC_L_inter; i++)
    st->old_exc[i] = st->old_exc[i+TSC_L_frame];

  return;
}





//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- SUB_SC_D.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------




/************************************************************************
*
*	FILENAME		:	sub_sc_d.c
*
*	DESCRIPTION		:	Source coder/decoder sub-routines used in the
*					TETRA speech codec
*
************************************************************************
*
*	SUB-ROUTINES	:	- Bits2prm_Tetra(st, )
*					- Cal_Rr2()
*					- Clsp_334()
*					- Dec_Ener()
*					- D4i60_16()
*					- D_D4i60()
*					- D_Lsp334(st, )
*					- Ener_Qua()
*					- G_Code()
*					- G_Pitch()
*					- Inter8_M1_3()
*					- Inter8_1_3()
*					- Inter32_M1_3()
*					- Inter32_1_3()
*					- Lag_Max()
*					- Norm_Corr()
*					- Pitch_Fr()
*					- Pitch_Ol_Dec()
*					- Post_Process(st, )
*					- Pred_Lt()
*					- Pre_Process() (including Init_Pre_Process)
*					- Prm2bits_Tetra(st, )
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*
*					clsp_334.tab
*					ener_qua.tab
*
************************************************************************/

//#include "source.h"

/************************************************************************
*
*	FILENAME		:	clsp_334.tab
*
*	DESCRIPTION		:	LSP codebooks:  split VQ with 3 codebooks.
*					The LSPs are in the cosine domain (-1, 1) in Q15
*
* 					codebook    	vector dimension        number of vectors
*					~~~~~~~~   	~~~~~~~~~~~~~         ~~~~~~~~~~~~~~
*					        1			3		256	(8 bits)
*					        2			3		512	(9 bits)
*					        3			4		512	(9 bits)
*
************************************************************************/

#define TSC_dim_dic1		3
#define TSC_dim_dic2		3
#define TSC_dim_dic3		4
#define TSC_taille_dic1	256
#define TSC_taille_dic2	512
#define TSC_taille_dic3	512

/*-----------------------------------------*
 * 1st codebook:   lsp0 to lsp2            *
 *-----------------------------------------*/

static const Word16 dico1_clsp[TSC_taille_dic1*TSC_dim_dic1] = {
      16231 * 2,     16091 * 2,     15682 * 2,
      16198 * 2,     16032 * 2,     15577 * 2,
      16142 * 2,     15986 * 2,     15536 * 2,
      16190 * 2,     15922 * 2,     15515 * 2,
      16188 * 2,     16002 * 2,     15314 * 2,
      16099 * 2,     15919 * 2,     15368 * 2,
      16149 * 2,     15941 * 2,     15254 * 2,
      16228 * 2,     16079 * 2,     15245 * 2,
      16127 * 2,     15802 * 2,     15316 * 2,
      16151 * 2,     15855 * 2,     15122 * 2,
      16056 * 2,     15878 * 2,     15137 * 2,
      16092 * 2,     15782 * 2,     15027 * 2,
      16047 * 2,     15719 * 2,     15263 * 2,
      16070 * 2,     15640 * 2,     15074 * 2,
      16104 * 2,     15877 * 2,     14778 * 2,
      16118 * 2,     15740 * 2,     14805 * 2,
      15978 * 2,     15756 * 2,     14910 * 2,
      16058 * 2,     15653 * 2,     14737 * 2,
      16023 * 2,     15562 * 2,     14896 * 2,
      15923 * 2,     15643 * 2,     15103 * 2,
      15963 * 2,     15463 * 2,     15048 * 2,
      15964 * 2,     15472 * 2,     14741 * 2,
      16027 * 2,     15539 * 2,     14592 * 2,
      15822 * 2,     15577 * 2,     14759 * 2,
      15894 * 2,     15657 * 2,     14525 * 2,
      15897 * 2,     15322 * 2,     14813 * 2,
      15768 * 2,     15447 * 2,     14793 * 2,
      15970 * 2,     15379 * 2,     14477 * 2,
      15827 * 2,     15527 * 2,     14303 * 2,
      15798 * 2,     15189 * 2,     14552 * 2,
      15716 * 2,     15424 * 2,     14318 * 2,
      15684 * 2,     15270 * 2,     14723 * 2,
      15862 * 2,     15324 * 2,     14172 * 2,
      15606 * 2,     15259 * 2,     14402 * 2,
      15878 * 2,     15108 * 2,     14312 * 2,
      15675 * 2,     15173 * 2,     14035 * 2,
      15568 * 2,     15011 * 2,     14381 * 2,
      15475 * 2,     15137 * 2,     14183 * 2,
      15702 * 2,     15347 * 2,     13870 * 2,
      15791 * 2,     15147 * 2,     13806 * 2,
      15568 * 2,     15002 * 2,     13834 * 2,
      15558 * 2,     15249 * 2,     13725 * 2,
      15802 * 2,     14842 * 2,     14041 * 2,
      15747 * 2,     14957 * 2,     13522 * 2,
      15368 * 2,     14888 * 2,     13942 * 2,
      15537 * 2,     14783 * 2,     13594 * 2,
      15476 * 2,     15115 * 2,     13393 * 2,
      15563 * 2,     14924 * 2,     13195 * 2,
      15383 * 2,     14829 * 2,     13377 * 2,
      15328 * 2,     15000 * 2,     13392 * 2,
      15505 * 2,     14658 * 2,     13177 * 2,
      15221 * 2,     14686 * 2,     13400 * 2,
      15501 * 2,     14450 * 2,     13469 * 2,
      15338 * 2,     14510 * 2,     13020 * 2,
      15232 * 2,     14720 * 2,     12775 * 2,
      15349 * 2,     14936 * 2,     12698 * 2,
      15085 * 2,     14466 * 2,     13087 * 2,
      15572 * 2,     14715 * 2,     12716 * 2,
      15259 * 2,     14332 * 2,     12698 * 2,
      15380 * 2,     14205 * 2,     13160 * 2,
      15383 * 2,     14518 * 2,     12244 * 2,
      15641 * 2,     14328 * 2,     12633 * 2,
      15501 * 2,     14832 * 2,     12172 * 2,
      15081 * 2,     14564 * 2,     12175 * 2,
      15063 * 2,     14095 * 2,     12747 * 2,
      15528 * 2,     14047 * 2,     12215 * 2,
      15315 * 2,     14226 * 2,     11850 * 2,
      14982 * 2,     14294 * 2,     12094 * 2,
      15123 * 2,     13908 * 2,     12244 * 2,
      15280 * 2,     14519 * 2,     11473 * 2,
      15225 * 2,     14844 * 2,     11802 * 2,
      15093 * 2,     13935 * 2,     11431 * 2,
      15482 * 2,     13924 * 2,     11356 * 2,
      15475 * 2,     14577 * 2,     11183 * 2,
      15218 * 2,     14226 * 2,     10791 * 2,
      14737 * 2,     14154 * 2,     11457 * 2,
      14916 * 2,     14278 * 2,     10548 * 2,
      15167 * 2,     14584 * 2,     10383 * 2,
      15441 * 2,     14281 * 2,     10338 * 2,
      15217 * 2,     13701 * 2,     10458 * 2,
      14864 * 2,     13740 * 2,     10632 * 2,
      15158 * 2,     13956 * 2,      9868 * 2,
      15160 * 2,     14387 * 2,      9550 * 2,
      14820 * 2,     13827 * 2,      9553 * 2,
      15307 * 2,     14034 * 2,      9137 * 2,
      15561 * 2,     13729 * 2,      9617 * 2,
      15121 * 2,     13366 * 2,      9420 * 2,
      14834 * 2,     13347 * 2,      9964 * 2,
      15075 * 2,     13592 * 2,      8793 * 2,
      15032 * 2,     14126 * 2,      8529 * 2,
      14716 * 2,     13670 * 2,      8579 * 2,
      15158 * 2,     13141 * 2,      8223 * 2,
      14584 * 2,     13155 * 2,      8860 * 2,
      14934 * 2,     12714 * 2,      8774 * 2,
      14725 * 2,     12847 * 2,      9653 * 2,
      14797 * 2,     12944 * 2,      7872 * 2,
      15022 * 2,     13701 * 2,      7746 * 2,
      14460 * 2,     12456 * 2,      8752 * 2,
      14271 * 2,     12719 * 2,      7966 * 2,
      14616 * 2,     12296 * 2,      7929 * 2,
      14837 * 2,     11931 * 2,      9091 * 2,
      15162 * 2,     12854 * 2,      7074 * 2,
      14781 * 2,     13408 * 2,      7008 * 2,
      14513 * 2,     12747 * 2,      6582 * 2,
      14328 * 2,     11984 * 2,      6857 * 2,
      14732 * 2,     11567 * 2,      7593 * 2,
      14754 * 2,     11677 * 2,      6455 * 2,
      14795 * 2,     12447 * 2,      5975 * 2,
      13850 * 2,     11888 * 2,      7609 * 2,
      14927 * 2,     13097 * 2,      5461 * 2,
      14269 * 2,     11440 * 2,      5252 * 2,
      14752 * 2,     11489 * 2,      4625 * 2,
      15290 * 2,     11792 * 2,      5357 * 2,
      14609 * 2,     12271 * 2,      4241 * 2,
      14387 * 2,     11193 * 2,      3749 * 2,
      14900 * 2,     12992 * 2,      3558 * 2,
      14865 * 2,     12377 * 2,      2905 * 2,
      14600 * 2,     11232 * 2,      2715 * 2,
      15425 * 2,     10659 * 2,      3564 * 2,
      14985 * 2,     13418 * 2,      4349 * 2,
      15424 * 2,     13710 * 2,      2694 * 2,
      15248 * 2,     13980 * 2,      4419 * 2,
      15139 * 2,     13555 * 2,      6077 * 2,
      15878 * 2,     13355 * 2,      4567 * 2,
      15740 * 2,     13445 * 2,      6142 * 2,
      15483 * 2,     14469 * 2,      5488 * 2,
      15404 * 2,     13842 * 2,      7138 * 2,
      15091 * 2,     14172 * 2,      6855 * 2,
      15514 * 2,     14415 * 2,      7464 * 2,
      15828 * 2,     14810 * 2,      6056 * 2,
      15905 * 2,     14299 * 2,      7695 * 2,
      15616 * 2,     13743 * 2,      8439 * 2,
      15349 * 2,     14546 * 2,      8611 * 2,
      15768 * 2,     15204 * 2,      7352 * 2,
      15718 * 2,     14835 * 2,      8847 * 2,
      15715 * 2,     14298 * 2,      9159 * 2,
      15922 * 2,     15121 * 2,      8926 * 2,
      15659 * 2,     15242 * 2,      9593 * 2,
      15864 * 2,     14260 * 2,      9965 * 2,
      15598 * 2,     14749 * 2,     10202 * 2,
      15347 * 2,     14862 * 2,      9903 * 2,
      15832 * 2,     14938 * 2,     10314 * 2,
      16020 * 2,     15210 * 2,     10351 * 2,
      15772 * 2,     14601 * 2,     11044 * 2,
      15781 * 2,     15501 * 2,     10547 * 2,
      15812 * 2,     15191 * 2,     11258 * 2,
      15637 * 2,     14955 * 2,     11529 * 2,
      15658 * 2,     15382 * 2,     11479 * 2,
      15534 * 2,     15204 * 2,     11639 * 2,
      15943 * 2,     14965 * 2,     11691 * 2,
      15955 * 2,     15373 * 2,     11636 * 2,
      15413 * 2,     15029 * 2,     11640 * 2,
      15761 * 2,     15509 * 2,     12201 * 2,
      15834 * 2,     15090 * 2,     12389 * 2,
      15737 * 2,     15311 * 2,     12584 * 2,
      15857 * 2,     15616 * 2,     11897 * 2,
      15889 * 2,     15436 * 2,     12566 * 2,
      16008 * 2,     15590 * 2,     12292 * 2,
      16017 * 2,     15262 * 2,     12780 * 2,
      15540 * 2,     15106 * 2,     12672 * 2,
      15929 * 2,     15639 * 2,     12941 * 2,
      15699 * 2,     15448 * 2,     13140 * 2,
      15955 * 2,     15459 * 2,     13281 * 2,
      15984 * 2,     15782 * 2,     12858 * 2,
      15826 * 2,     15311 * 2,     13396 * 2,
      15608 * 2,     15291 * 2,     13082 * 2,
      15725 * 2,     15089 * 2,     13155 * 2,
      15790 * 2,     15525 * 2,     13586 * 2,
      15861 * 2,     15622 * 2,     13596 * 2,
      15960 * 2,     15103 * 2,     13294 * 2,
      15950 * 2,     15421 * 2,     13857 * 2,
      15938 * 2,     15709 * 2,     13665 * 2,
      16045 * 2,     15641 * 2,     13517 * 2,
      16113 * 2,     15435 * 2,     13783 * 2,
      15950 * 2,     15579 * 2,     14127 * 2,
      16032 * 2,     15763 * 2,     13984 * 2,
      16056 * 2,     15597 * 2,     14271 * 2,
      16071 * 2,     15282 * 2,     14159 * 2,
      15973 * 2,     15723 * 2,     14379 * 2,
      16086 * 2,     15755 * 2,     14429 * 2,
      16119 * 2,     15881 * 2,     14090 * 2,
      16205 * 2,     15680 * 2,     14372 * 2,
      16028 * 2,     15842 * 2,     14528 * 2,
      16199 * 2,     15761 * 2,     13813 * 2,
      16141 * 2,     15959 * 2,     14662 * 2,
      16198 * 2,     16013 * 2,     14080 * 2,
      16233 * 2,     16079 * 2,     14465 * 2,
      16203 * 2,     15948 * 2,     14878 * 2,
      16271 * 2,     16051 * 2,     13866 * 2,
      16287 * 2,     16154 * 2,     14482 * 2,
      16167 * 2,     15960 * 2,     13393 * 2,
      16267 * 2,     15940 * 2,     13407 * 2,
      16076 * 2,     15830 * 2,     13540 * 2,
      16237 * 2,     15561 * 2,     13288 * 2,
      16200 * 2,     15808 * 2,     12800 * 2,
      16126 * 2,     15534 * 2,     13042 * 2,
      16082 * 2,     15783 * 2,     12455 * 2,
      16139 * 2,     15482 * 2,     12064 * 2,
      16155 * 2,     15840 * 2,     11644 * 2,
      16238 * 2,     16026 * 2,     11633 * 2,
      16071 * 2,     15639 * 2,     11117 * 2,
      16186 * 2,     15311 * 2,     11146 * 2,
      15938 * 2,     15700 * 2,     10741 * 2,
      15948 * 2,     15495 * 2,     10329 * 2,
      16108 * 2,     15736 * 2,      9800 * 2,
      16056 * 2,     14735 * 2,     11245 * 2,
      15830 * 2,     15467 * 2,      9241 * 2,
      16167 * 2,     15951 * 2,      8972 * 2,
      16051 * 2,     15570 * 2,      8007 * 2,
      15888 * 2,     15592 * 2,      7559 * 2,
      15978 * 2,     15368 * 2,      7032 * 2,
      15846 * 2,     15467 * 2,      5711 * 2,
      16014 * 2,     15082 * 2,      4645 * 2,
      15831 * 2,     14913 * 2,      3458 * 2,
      15755 * 2,     13822 * 2,     10633 * 2,
      15883 * 2,     14259 * 2,     11598 * 2,
      15463 * 2,     13031 * 2,     10169 * 2,
      15821 * 2,     12195 * 2,      8137 * 2,
      15582 * 2,     12516 * 2,      9433 * 2,
      14821 * 2,     12969 * 2,     10527 * 2,
      15732 * 2,     13268 * 2,     11498 * 2,
      15146 * 2,     13327 * 2,     11267 * 2,
      14842 * 2,     12289 * 2,     10422 * 2,
      14567 * 2,     13505 * 2,     10780 * 2,
      14390 * 2,     12759 * 2,     11170 * 2,
      14858 * 2,     13614 * 2,     11657 * 2,
      14237 * 2,     13218 * 2,     11461 * 2,
      14608 * 2,     13504 * 2,     12080 * 2,
      14674 * 2,     13974 * 2,     12184 * 2,
      14866 * 2,     13720 * 2,     12773 * 2,
      14896 * 2,     14332 * 2,     12782 * 2,
      14840 * 2,     14067 * 2,     13132 * 2,
      15435 * 2,     13786 * 2,     12823 * 2,
      15121 * 2,     14302 * 2,     13483 * 2,
      15333 * 2,     14603 * 2,     13872 * 2,
      15772 * 2,     14524 * 2,     13052 * 2,
      15708 * 2,     14584 * 2,     12009 * 2,
      15797 * 2,     14669 * 2,     13637 * 2,
      15814 * 2,     14858 * 2,     12904 * 2,
      15962 * 2,     14568 * 2,     12497 * 2,
      15969 * 2,     15014 * 2,     13792 * 2,
      14039 * 2,     12511 * 2,      9773 * 2,
      14251 * 2,     11878 * 2,      9621 * 2,
      14100 * 2,     11691 * 2,      8583 * 2,
      14814 * 2,     10471 * 2,      7911 * 2,
      13965 * 2,     10725 * 2,      7012 * 2,
      15689 * 2,      9981 * 2,      6000 * 2,
      14388 * 2,      9740 * 2,      5769 * 2,
      14384 * 2,      9736 * 2,      4341 * 2,
      13970 * 2,      9084 * 2,      5028 * 2,
      15298 * 2,      8857 * 2,      4205 * 2,
      14071 * 2,      7942 * 2,      4472 * 2,
      13980 * 2,      9604 * 2,      3273 * 2,
      13383 * 2,      8140 * 2,      3828 * 2,
      13005 * 2,      7620 * 2,      4742 * 2,
      16310 * 2,     16203 * 2,     15451 * 2 };


/*-----------------------------------------*
 * 2nd codebook:   lsp3 to lsp5            *
 *-----------------------------------------*/

static const Word16 dico2_clsp[TSC_taille_dic2*TSC_dim_dic2] = {

       7564 * 2,      -293 * 2,     -2543 * 2,
       7820 * 2,      -150 * 2,     -1813 * 2,
       6913 * 2,     -1006 * 2,     -2297 * 2,
       7099 * 2,      -498 * 2,     -3337 * 2,
       6526 * 2,       308 * 2,     -2075 * 2,
       6953 * 2,       958 * 2,     -2388 * 2,
       6811 * 2,       758 * 2,     -3159 * 2,
       6448 * 2,        93 * 2,     -1294 * 2,
       7572 * 2,       843 * 2,     -1485 * 2,
       5803 * 2,        76 * 2,     -2785 * 2,
       5791 * 2,      -566 * 2,     -1850 * 2,
       5655 * 2,       944 * 2,     -1852 * 2,
       6035 * 2,      1575 * 2,     -2728 * 2,
       5226 * 2,       672 * 2,     -2515 * 2,
       6131 * 2,      1800 * 2,     -1816 * 2,
       6104 * 2,      1092 * 2,      -955 * 2,
       4868 * 2,       150 * 2,     -1704 * 2,
       4540 * 2,       928 * 2,     -2198 * 2,
       4954 * 2,      1589 * 2,     -1403 * 2,
       5133 * 2,      1825 * 2,     -2481 * 2,
       5266 * 2,      2505 * 2,     -2018 * 2,
       4140 * 2,      1163 * 2,     -1467 * 2,
       5312 * 2,      2208 * 2,      -846 * 2,
       4412 * 2,      2344 * 2,     -1137 * 2,
       4016 * 2,      1694 * 2,     -1020 * 2,
       3937 * 2,      2428 * 2,     -2025 * 2,
       4706 * 2,      3059 * 2,     -1381 * 2,
       3469 * 2,      1623 * 2,     -1956 * 2,
       4611 * 2,      2818 * 2,      -273 * 2,
       3830 * 2,      2126 * 2,      -140 * 2,
       5232 * 2,      1681 * 2,      -440 * 2,
       4523 * 2,      1482 * 2,       106 * 2,
       4461 * 2,       950 * 2,      -458 * 2,
       3228 * 2,      1395 * 2,      -718 * 2,
       4792 * 2,      2283 * 2,       392 * 2,
       4198 * 2,       318 * 2,      -990 * 2,
       5397 * 2,      3081 * 2,      -622 * 2,
       6034 * 2,      1867 * 2,       379 * 2,
       6276 * 2,      1978 * 2,      -806 * 2,
       5881 * 2,      2678 * 2,       452 * 2,
       6308 * 2,      2449 * 2,      -240 * 2,
       5782 * 2,      2367 * 2,      1051 * 2,
       6319 * 2,      3243 * 2,       153 * 2,
       6434 * 2,      2971 * 2,      -955 * 2,
       7053 * 2,      1856 * 2,      -216 * 2,
       7385 * 2,      2650 * 2,        64 * 2,
       6774 * 2,      2930 * 2,      1008 * 2,
       7416 * 2,      2162 * 2,       703 * 2,
       7282 * 2,      2723 * 2,      -833 * 2,
       7452 * 2,      3624 * 2,      -298 * 2,
       8031 * 2,      3206 * 2,       530 * 2,
       6889 * 2,      4018 * 2,       470 * 2,
       7485 * 2,      3735 * 2,      1191 * 2,
       7008 * 2,      2808 * 2,      1663 * 2,
       8040 * 2,      3276 * 2,      1625 * 2,
       7841 * 2,      4506 * 2,       984 * 2,
       8460 * 2,      4070 * 2,       391 * 2,
       6611 * 2,      4137 * 2,      1514 * 2,
       7487 * 2,      4459 * 2,      2017 * 2,
       7549 * 2,      3763 * 2,      2384 * 2,
       6420 * 2,      3454 * 2,      1907 * 2,
       6694 * 2,      4869 * 2,      1313 * 2,
       7852 * 2,      5242 * 2,      1730 * 2,
       6991 * 2,      5489 * 2,      1418 * 2,
       6873 * 2,      4940 * 2,      2548 * 2,
       6389 * 2,      4265 * 2,      2477 * 2,
       7703 * 2,      4357 * 2,      2885 * 2,
       7318 * 2,      5438 * 2,      2788 * 2,
       7257 * 2,      4987 * 2,      3443 * 2,
       8222 * 2,      5222 * 2,      2614 * 2,
       7645 * 2,      6104 * 2,      2139 * 2,
       8148 * 2,      5663 * 2,      3485 * 2,
       8607 * 2,      4884 * 2,      3527 * 2,
       7512 * 2,      5913 * 2,      4110 * 2,
       8568 * 2,      6064 * 2,      2396 * 2,
       8113 * 2,      6670 * 2,      2267 * 2,
       8790 * 2,      6611 * 2,      3338 * 2,
       9197 * 2,      5692 * 2,      3316 * 2,
       8272 * 2,      6549 * 2,      4051 * 2,
       9059 * 2,      6038 * 2,      4412 * 2,
       9554 * 2,      6359 * 2,      3624 * 2,
       9264 * 2,      5327 * 2,      4105 * 2,
       9149 * 2,      6926 * 2,      4494 * 2,
       9916 * 2,      6778 * 2,      4197 * 2,
       9600 * 2,      7338 * 2,      3435 * 2,
      10366 * 2,      5712 * 2,      3980 * 2,
      10362 * 2,      6304 * 2,      4858 * 2,
      10383 * 2,      7083 * 2,      4912 * 2,
      10467 * 2,      7485 * 2,      4027 * 2,
      10573 * 2,      6057 * 2,      3302 * 2,
      10367 * 2,      6987 * 2,      3015 * 2,
      11091 * 2,      6561 * 2,      3914 * 2,
      11279 * 2,      7390 * 2,      4289 * 2,
      11378 * 2,      5806 * 2,      4560 * 2,
      11504 * 2,      6657 * 2,      5030 * 2,
      11919 * 2,      6407 * 2,      3700 * 2,
      11338 * 2,      5343 * 2,      3847 * 2,
      11333 * 2,      5308 * 2,      3155 * 2,
      11051 * 2,      6450 * 2,      2535 * 2,
      11518 * 2,      7251 * 2,      2774 * 2,
      12044 * 2,      6329 * 2,      2477 * 2,
      12482 * 2,      5821 * 2,      4097 * 2,
      12478 * 2,      7102 * 2,      3587 * 2,
      12676 * 2,      5586 * 2,      3333 * 2,
      12171 * 2,      4918 * 2,      2877 * 2,
      12751 * 2,      6844 * 2,      2736 * 2,
      13176 * 2,      5963 * 2,      3744 * 2,
      12802 * 2,      6393 * 2,      4533 * 2,
      12692 * 2,      4994 * 2,      4094 * 2,
      13017 * 2,      4667 * 2,      3408 * 2,
      13343 * 2,      5710 * 2,      2375 * 2,
      12868 * 2,      5296 * 2,      2008 * 2,
      13364 * 2,      4239 * 2,      2361 * 2,
      13843 * 2,      5743 * 2,      4255 * 2,
      14049 * 2,      5097 * 2,      3468 * 2,
      14033 * 2,      4939 * 2,      2835 * 2,
      13898 * 2,      6536 * 2,      3115 * 2,
      14169 * 2,      4961 * 2,      1705 * 2,
      13710 * 2,      4832 * 2,      1148 * 2,
      14596 * 2,      5508 * 2,      2662 * 2,
      13187 * 2,      6223 * 2,       960 * 2,
      13822 * 2,      6583 * 2,       583 * 2,
      14436 * 2,      7264 * 2,      1534 * 2,
      13770 * 2,      7598 * 2,      1965 * 2,
      13170 * 2,      7340 * 2,      1891 * 2,
      13941 * 2,      4561 * 2,       129 * 2,
      13431 * 2,      5471 * 2,      -449 * 2,
      14335 * 2,      5782 * 2,      -545 * 2,
      12866 * 2,      5691 * 2,      -477 * 2,
      12509 * 2,      6194 * 2,      1001 * 2,
      12692 * 2,      7437 * 2,       704 * 2,
      12722 * 2,      4956 * 2,       901 * 2,
      12151 * 2,      6291 * 2,       -75 * 2,
      12353 * 2,      4481 * 2,       -63 * 2,
      11761 * 2,      6052 * 2,       985 * 2,
      12284 * 2,      4617 * 2,      1797 * 2,
      11734 * 2,      4431 * 2,       388 * 2,
      11295 * 2,      5975 * 2,      -390 * 2,
      11032 * 2,      4944 * 2,       743 * 2,
      11399 * 2,      5089 * 2,      2116 * 2,
      10964 * 2,      6187 * 2,      1187 * 2,
      10863 * 2,      4166 * 2,      1606 * 2,
      10205 * 2,      5235 * 2,      1461 * 2,
      10772 * 2,      4436 * 2,      -223 * 2,
      10305 * 2,      5642 * 2,        70 * 2,
      10393 * 2,      6369 * 2,      1675 * 2,
      10597 * 2,      4319 * 2,      2393 * 2,
      10600 * 2,      5163 * 2,      2801 * 2,
       9741 * 2,      4539 * 2,       336 * 2,
       9347 * 2,      4692 * 2,      1139 * 2,
       9785 * 2,      3625 * 2,      1401 * 2,
       9006 * 2,      5480 * 2,      1621 * 2,
       9538 * 2,      5265 * 2,      2600 * 2,
       8992 * 2,      4694 * 2,      2420 * 2,
       8956 * 2,      4140 * 2,      1972 * 2,
       9472 * 2,      3917 * 2,      2709 * 2,
       9301 * 2,      4477 * 2,      3132 * 2,
       8874 * 2,      3348 * 2,      2170 * 2,
       9827 * 2,      5042 * 2,      3445 * 2,
       9890 * 2,      6059 * 2,      2463 * 2,
      10518 * 2,      3522 * 2,      2283 * 2,
      10797 * 2,      4366 * 2,      3190 * 2,
      10728 * 2,      4775 * 2,      3680 * 2,
       9809 * 2,      2897 * 2,      1746 * 2,
      11591 * 2,      3608 * 2,      1929 * 2,
      10834 * 2,      3491 * 2,       761 * 2,
      11979 * 2,      3948 * 2,      2794 * 2,
      11799 * 2,      2977 * 2,      1849 * 2,
      11804 * 2,      3224 * 2,      1025 * 2,
      10578 * 2,      2760 * 2,       924 * 2,
      12052 * 2,      2472 * 2,      1517 * 2,
      11763 * 2,      1897 * 2,       874 * 2,
      10935 * 2,      2556 * 2,       -10 * 2,
      12339 * 2,      2548 * 2,       178 * 2,
      11812 * 2,      2374 * 2,      -613 * 2,
      10869 * 2,      1311 * 2,       203 * 2,
      11242 * 2,      3714 * 2,      -618 * 2,
      10672 * 2,      1941 * 2,     -1043 * 2,
      11974 * 2,       768 * 2,      -371 * 2,
      10582 * 2,      1020 * 2,      -737 * 2,
      11222 * 2,      1190 * 2,     -1751 * 2,
      11652 * 2,      3167 * 2,     -1686 * 2,
      11771 * 2,      1803 * 2,     -2377 * 2,
      12157 * 2,       276 * 2,     -1261 * 2,
      11281 * 2,      1351 * 2,     -3052 * 2,
      11760 * 2,        77 * 2,     -1978 * 2,
      10461 * 2,       599 * 2,     -2519 * 2,
      12320 * 2,      -312 * 2,     -2603 * 2,
      12908 * 2,      1043 * 2,     -1899 * 2,
      12826 * 2,       931 * 2,     -3215 * 2,
      12987 * 2,      1269 * 2,     -1053 * 2,
      13424 * 2,       538 * 2,     -2342 * 2,
      13102 * 2,       -57 * 2,     -4076 * 2,
      12263 * 2,      2454 * 2,     -3882 * 2,
      13816 * 2,       796 * 2,     -3719 * 2,
      13973 * 2,      1713 * 2,     -2871 * 2,
      13869 * 2,      1446 * 2,     -1618 * 2,
      13292 * 2,      2675 * 2,     -1898 * 2,
      13868 * 2,      2841 * 2,     -1281 * 2,
      14411 * 2,      2364 * 2,     -3252 * 2,
      13521 * 2,      3517 * 2,     -3271 * 2,
      14425 * 2,      2348 * 2,     -1729 * 2,
      14732 * 2,      2891 * 2,     -2389 * 2,
      14526 * 2,      3821 * 2,     -1287 * 2,
      14039 * 2,      4544 * 2,     -2746 * 2,
      13926 * 2,      5111 * 2,     -1292 * 2,
      13529 * 2,      3616 * 2,      -239 * 2,
      14613 * 2,      5736 * 2,     -2230 * 2,
      14968 * 2,      5253 * 2,     -1021 * 2,
      14559 * 2,      3776 * 2,       465 * 2,
      14617 * 2,      4378 * 2,     -4224 * 2,
      13368 * 2,      6485 * 2,     -1784 * 2,
      12599 * 2,      4340 * 2,     -1435 * 2,
      12694 * 2,      3951 * 2,     -2784 * 2,
      13073 * 2,      3789 * 2,       588 * 2,
      12852 * 2,      2836 * 2,      -620 * 2,
      13700 * 2,      1906 * 2,       -33 * 2,
      13417 * 2,      2505 * 2,       809 * 2,
      12998 * 2,      3554 * 2,      1467 * 2,
      14015 * 2,      3141 * 2,      1469 * 2,
      12867 * 2,      3536 * 2,      2271 * 2,
      13756 * 2,       929 * 2,      -384 * 2,
      10188 * 2,      1745 * 2,       592 * 2,
      11419 * 2,      5535 * 2,     -1711 * 2,
       9637 * 2,      3260 * 2,       332 * 2,
      10133 * 2,      4611 * 2,     -1013 * 2,
       9707 * 2,      2996 * 2,      -871 * 2,
      10284 * 2,      3161 * 2,     -2270 * 2,
       9710 * 2,      1748 * 2,      -461 * 2,
       9184 * 2,      2352 * 2,       138 * 2,
       9166 * 2,      3710 * 2,     -1743 * 2,
       9163 * 2,      4407 * 2,      -475 * 2,
       8561 * 2,      3500 * 2,      -614 * 2,
       8676 * 2,      2539 * 2,     -1572 * 2,
       8520 * 2,      2158 * 2,      -676 * 2,
       7797 * 2,      3393 * 2,     -1237 * 2,
       7991 * 2,      1776 * 2,     -1668 * 2,
       8159 * 2,      3239 * 2,     -2484 * 2,
       8140 * 2,      4597 * 2,      -982 * 2,
       6990 * 2,      2510 * 2,     -1870 * 2,
       7026 * 2,      4004 * 2,     -1953 * 2,
       7136 * 2,      1609 * 2,     -1203 * 2,
       7217 * 2,      1904 * 2,     -2725 * 2,
       6773 * 2,      3042 * 2,     -2878 * 2,
       5944 * 2,      2926 * 2,     -1591 * 2,
       6396 * 2,      3845 * 2,      -761 * 2,
       5666 * 2,      3021 * 2,     -3053 * 2,
       5164 * 2,      3788 * 2,     -1911 * 2,
       5680 * 2,      4283 * 2,     -1046 * 2,
       6350 * 2,      4968 * 2,     -1518 * 2,
       5489 * 2,      3743 * 2,       -28 * 2,
       6183 * 2,      4396 * 2,       173 * 2,
       7036 * 2,      4742 * 2,      -310 * 2,
       5005 * 2,      3366 * 2,       513 * 2,
       5862 * 2,      3686 * 2,      1244 * 2,
       5779 * 2,      4527 * 2,      1341 * 2,
       5432 * 2,      3100 * 2,      1205 * 2,
       7416 * 2,      5187 * 2,       465 * 2,
       7287 * 2,      5658 * 2,      -643 * 2,
       8414 * 2,      5026 * 2,       586 * 2,
       7991 * 2,      5959 * 2,       759 * 2,
       8790 * 2,      2939 * 2,      1176 * 2,
       9149 * 2,      5825 * 2,       175 * 2,
       8623 * 2,      6672 * 2,      1016 * 2,
       8299 * 2,      6034 * 2,     -1339 * 2,
       9158 * 2,      6865 * 2,      -646 * 2,
       9716 * 2,      6566 * 2,       979 * 2,
       9119 * 2,      7508 * 2,      1103 * 2,
       9216 * 2,      6732 * 2,      2187 * 2,
       9767 * 2,      7791 * 2,       192 * 2,
      10117 * 2,      7427 * 2,      1579 * 2,
       9646 * 2,      8099 * 2,      2184 * 2,
      10789 * 2,      7148 * 2,      -115 * 2,
      10538 * 2,      8395 * 2,       577 * 2,
      10523 * 2,      7909 * 2,      2683 * 2,
      11293 * 2,      7494 * 2,      1440 * 2,
      11121 * 2,      8682 * 2,      1598 * 2,
      10387 * 2,      8722 * 2,      2843 * 2,
      11092 * 2,      7990 * 2,      3176 * 2,
      11790 * 2,      8511 * 2,      2266 * 2,
      10995 * 2,      9415 * 2,      2452 * 2,
      11666 * 2,      8212 * 2,      3760 * 2,
      11999 * 2,      9089 * 2,      3774 * 2,
      12268 * 2,      9210 * 2,      2028 * 2,
      11821 * 2,      9800 * 2,      3288 * 2,
      11424 * 2,      9183 * 2,      4356 * 2,
      12664 * 2,      9221 * 2,      3093 * 2,
      12343 * 2,      9929 * 2,      3844 * 2,
      11318 * 2,     10061 * 2,      4382 * 2,
      12362 * 2,      9114 * 2,      4908 * 2,
      12753 * 2,     10503 * 2,      3698 * 2,
      12985 * 2,      9731 * 2,      4733 * 2,
      11892 * 2,      9626 * 2,      5588 * 2,
      12458 * 2,     10370 * 2,      5601 * 2,
      11854 * 2,     10622 * 2,      5262 * 2,
      12878 * 2,      8886 * 2,      5806 * 2,
      12559 * 2,      9585 * 2,      6659 * 2,
      13068 * 2,     10386 * 2,      6203 * 2,
      12633 * 2,     11087 * 2,      6026 * 2,
      12186 * 2,     10894 * 2,      6938 * 2,
      11920 * 2,      9993 * 2,      7144 * 2,
      13104 * 2,     11211 * 2,      6967 * 2,
      13164 * 2,     10142 * 2,      7665 * 2,
      12367 * 2,      9323 * 2,      7804 * 2,
      12585 * 2,     10879 * 2,      8219 * 2,
      12317 * 2,     10251 * 2,      8449 * 2,
      12720 * 2,     11567 * 2,      8102 * 2,
      13234 * 2,      9369 * 2,      8463 * 2,
      11548 * 2,      9421 * 2,      8254 * 2,
      13249 * 2,      8973 * 2,      7039 * 2,
      12176 * 2,      8597 * 2,      6758 * 2,
      11461 * 2,      9342 * 2,      6748 * 2,
      10742 * 2,      8788 * 2,      7518 * 2,
      10806 * 2,      9124 * 2,      6235 * 2,
      10997 * 2,      8275 * 2,      6159 * 2,
      11560 * 2,      8613 * 2,      5437 * 2,
      12283 * 2,      8066 * 2,      5739 * 2,
      11068 * 2,      7556 * 2,      5368 * 2,
      10110 * 2,      8304 * 2,      6289 * 2,
      10311 * 2,      7972 * 2,      5163 * 2,
      10441 * 2,      7206 * 2,      5823 * 2,
      10880 * 2,      8468 * 2,      4478 * 2,
      11908 * 2,      7316 * 2,      5254 * 2,
      10635 * 2,      9226 * 2,      4414 * 2,
       9984 * 2,      8256 * 2,      4205 * 2,
       9275 * 2,      7592 * 2,      5010 * 2,
      12135 * 2,      8005 * 2,      4117 * 2,
      12553 * 2,      7346 * 2,      4798 * 2,
      12699 * 2,      8151 * 2,      2816 * 2,
      13136 * 2,      8345 * 2,      4614 * 2,
      12192 * 2,      7482 * 2,      1678 * 2,
      13413 * 2,      7358 * 2,      3703 * 2,
      13403 * 2,      6763 * 2,      5157 * 2,
      13300 * 2,      7614 * 2,      5971 * 2,
      13940 * 2,      8404 * 2,      3784 * 2,
      13777 * 2,      8074 * 2,      5441 * 2,
      13465 * 2,      9081 * 2,      2751 * 2,
      13567 * 2,      9647 * 2,      4283 * 2,
      14428 * 2,      7463 * 2,      4258 * 2,
      14332 * 2,      9351 * 2,      4182 * 2,
      14393 * 2,      8218 * 2,      2867 * 2,
      14276 * 2,      9080 * 2,      6060 * 2,
      14076 * 2,      9711 * 2,      2470 * 2,
      13630 * 2,      9696 * 2,      5808 * 2,
      14679 * 2,      9709 * 2,      2242 * 2,
      14160 * 2,     10997 * 2,      4409 * 2,
      13309 * 2,     10327 * 2,      2581 * 2,
      13224 * 2,     10865 * 2,      4126 * 2,
      14552 * 2,     11321 * 2,      4526 * 2,
      13808 * 2,     11509 * 2,      4970 * 2,
      13475 * 2,     10927 * 2,      5633 * 2,
      13018 * 2,     11618 * 2,      4262 * 2,
      13931 * 2,     12153 * 2,      3241 * 2,
      13376 * 2,     12029 * 2,      5980 * 2,
      12334 * 2,     10978 * 2,      3108 * 2,
      13728 * 2,     10999 * 2,      1834 * 2,
      14883 * 2,     11658 * 2,      3750 * 2,
      13940 * 2,     10644 * 2,      6944 * 2,
      14216 * 2,     12662 * 2,      6727 * 2,
      13652 * 2,     12582 * 2,      7202 * 2,
      13647 * 2,     11561 * 2,      8042 * 2,
      13073 * 2,     12210 * 2,      7825 * 2,
      13021 * 2,     11362 * 2,      9396 * 2,
      13850 * 2,     12600 * 2,      9870 * 2,
      13321 * 2,     12442 * 2,     10862 * 2,
      12641 * 2,     10336 * 2,      1486 * 2,
      12993 * 2,     11301 * 2,      1024 * 2,
      12026 * 2,     10336 * 2,      1353 * 2,
      14249 * 2,     10493 * 2,      1174 * 2,
      12910 * 2,      8766 * 2,      1401 * 2,
      13576 * 2,      9065 * 2,       929 * 2,
      13237 * 2,     10168 * 2,        45 * 2,
      12659 * 2,      9572 * 2,      -237 * 2,
      12200 * 2,      8811 * 2,       539 * 2,
      11676 * 2,      9627 * 2,       701 * 2,
      13175 * 2,      7994 * 2,      -647 * 2,
      14274 * 2,      8614 * 2,      1030 * 2,
      11770 * 2,      7827 * 2,       402 * 2,
      12397 * 2,      7494 * 2,      -904 * 2,
      11242 * 2,      8488 * 2,      -714 * 2,
      11883 * 2,      8572 * 2,     -2307 * 2,
      12942 * 2,      8137 * 2,     -2821 * 2,
      14047 * 2,      7675 * 2,     -1923 * 2,
      12011 * 2,      6420 * 2,     -3262 * 2,
      10184 * 2,      7669 * 2,     -1818 * 2,
      10155 * 2,      6096 * 2,     -1676 * 2,
      13048 * 2,      5362 * 2,     -4666 * 2,
      13653 * 2,      6231 * 2,     -5460 * 2,
      11925 * 2,      4015 * 2,     -5430 * 2,
      10735 * 2,      3346 * 2,     -3524 * 2,
       9349 * 2,      4598 * 2,     -3399 * 2,
      11031 * 2,      2014 * 2,     -5157 * 2,
       9330 * 2,      2646 * 2,     -5465 * 2,
       9374 * 2,      1306 * 2,     -3874 * 2,
      11556 * 2,       299 * 2,     -4747 * 2,
       9264 * 2,      1787 * 2,     -2708 * 2,
      10938 * 2,       -42 * 2,     -3854 * 2,
       9077 * 2,       249 * 2,     -2965 * 2,
       9951 * 2,      -643 * 2,     -4696 * 2,
       8695 * 2,       -70 * 2,     -3650 * 2,
       8193 * 2,       664 * 2,     -4738 * 2,
       8925 * 2,      -879 * 2,     -2873 * 2,
       8632 * 2,       644 * 2,     -2251 * 2,
      10049 * 2,      -345 * 2,     -1940 * 2,
       8812 * 2,     -1699 * 2,     -4282 * 2,
       7585 * 2,     -1492 * 2,     -3374 * 2,
       7321 * 2,      -662 * 2,     -4259 * 2,
       8513 * 2,      -862 * 2,     -5444 * 2,
       6425 * 2,       426 * 2,     -4076 * 2,
       7536 * 2,      1957 * 2,     -3948 * 2,
       6558 * 2,      1461 * 2,     -4873 * 2,
       5999 * 2,      -265 * 2,     -5358 * 2,
       5657 * 2,       348 * 2,     -4658 * 2,
       5488 * 2,      1123 * 2,     -3537 * 2,
       5150 * 2,        51 * 2,     -3476 * 2,
       5222 * 2,      -943 * 2,     -4507 * 2,
       4581 * 2,       362 * 2,     -4293 * 2,
       4868 * 2,      1676 * 2,     -5157 * 2,
       4431 * 2,      1543 * 2,     -4063 * 2,
       3850 * 2,       845 * 2,     -4686 * 2,
       4059 * 2,      -107 * 2,     -3395 * 2,
       4527 * 2,      1273 * 2,     -2949 * 2,
       3735 * 2,       744 * 2,     -3093 * 2,
       3082 * 2,      1200 * 2,     -4127 * 2,
       3405 * 2,      1879 * 2,     -3231 * 2,
       4594 * 2,      2463 * 2,     -3711 * 2,
       2999 * 2,       103 * 2,     -3506 * 2,
       2380 * 2,       690 * 2,     -3330 * 2,
       3339 * 2,       658 * 2,     -2070 * 2,
       2691 * 2,      1302 * 2,     -1971 * 2,
       3346 * 2,       -52 * 2,     -1960 * 2,
       2814 * 2,      -501 * 2,     -3150 * 2,
       2052 * 2,        79 * 2,     -2475 * 2,
       2320 * 2,      -560 * 2,     -2143 * 2,
       2676 * 2,       627 * 2,     -1273 * 2,
       2770 * 2,     -1345 * 2,     -2640 * 2,
       4487 * 2,       171 * 2,     -2513 * 2,
       4378 * 2,      -681 * 2,     -2521 * 2,
       3230 * 2,     -1229 * 2,     -3412 * 2,
       2693 * 2,     -2022 * 2,     -3359 * 2,
       4176 * 2,      -888 * 2,     -3993 * 2,
       3523 * 2,     -2141 * 2,     -3983 * 2,
       3004 * 2,      -293 * 2,     -4569 * 2,
       3375 * 2,     -1360 * 2,     -4767 * 2,
       4100 * 2,     -1735 * 2,     -5277 * 2,
       4204 * 2,      -302 * 2,     -5248 * 2,
       2085 * 2,      -979 * 2,     -4878 * 2,
       2646 * 2,     -1469 * 2,     -5896 * 2,
       3291 * 2,      -357 * 2,     -5984 * 2,
       2098 * 2,       152 * 2,     -5496 * 2,
       3851 * 2,     -1336 * 2,     -6741 * 2,
       2411 * 2,     -1002 * 2,     -6998 * 2,
       4221 * 2,     -2346 * 2,     -6013 * 2,
       1591 * 2,     -1869 * 2,     -6561 * 2,
       1203 * 2,      -796 * 2,     -6075 * 2,
       1744 * 2,     -2522 * 2,     -5654 * 2,
       2310 * 2,     -3164 * 2,     -6344 * 2,
       1385 * 2,     -1923 * 2,     -4877 * 2,
       1045 * 2,     -2988 * 2,     -6886 * 2,
       1018 * 2,     -3544 * 2,     -5868 * 2,
       2581 * 2,     -3527 * 2,     -5340 * 2,
        441 * 2,     -2748 * 2,     -5208 * 2,
       1814 * 2,     -3242 * 2,     -4644 * 2,
       1316 * 2,     -3962 * 2,     -6696 * 2,
       1858 * 2,     -2301 * 2,     -4189 * 2,
        -49 * 2,     -3715 * 2,     -6801 * 2,
        -16 * 2,     -4707 * 2,     -6171 * 2,
       1248 * 2,     -4592 * 2,     -7268 * 2,
       1512 * 2,     -3646 * 2,     -7795 * 2,
      -1013 * 2,     -2955 * 2,     -6031 * 2,
        -92 * 2,     -2204 * 2,     -7396 * 2,
        255 * 2,     -4976 * 2,     -7893 * 2,
       3177 * 2,     -4114 * 2,     -6632 * 2,
       1920 * 2,     -5053 * 2,     -9114 * 2,
       3872 * 2,     -3048 * 2,     -7641 * 2,
       4455 * 2,     -3457 * 2,     -6789 * 2,
       3806 * 2,     -2934 * 2,     -4924 * 2,
       5509 * 2,     -2502 * 2,     -6742 * 2,
       5692 * 2,     -1390 * 2,     -5591 * 2,
       6772 * 2,     -2446 * 2,     -4848 * 2,
       5667 * 2,     -1925 * 2,     -3815 * 2,
       4441 * 2,       625 * 2,     -6528 * 2,
       5532 * 2,     -1052 * 2,     -3049 * 2,
       7978 * 2,      -858 * 2,     -6705 * 2,
       3643 * 2,      1692 * 2,     -5218 * 2,
       5975 * 2,      2763 * 2,     -4375 * 2,
       2361 * 2,       891 * 2,     -5213 * 2,
       1718 * 2,      -991 * 2,     -4008 * 2,
       1651 * 2,       -61 * 2,     -3633 * 2,
       4489 * 2,      3255 * 2,     -2982 * 2,
       1208 * 2,      -990 * 2,     -2970 * 2,
       -338 * 2,     -2055 * 2,     -4560 * 2,
       6498 * 2,      1132 * 2,      -164 * 2,
       7336 * 2,       403 * 2,      -675 * 2,
       8818 * 2,      1050 * 2,      -967 * 2,
       8652 * 2,      1072 * 2,      -293 * 2,
       8291 * 2,      1571 * 2,       303 * 2,
       9525 * 2,       235 * 2,     -1134 * 2,
       9652 * 2,      1049 * 2,     -1668 * 2,
       9069 * 2,      2296 * 2,       828 * 2,
       9526 * 2,      2265 * 2,      1359 * 2,
      10787 * 2,     -1109 * 2,     -3223 * 2,
       7767 * 2,      5486 * 2,     -3086 * 2,
      11690 * 2,     -1931 * 2,     -4657 * 2,
       8997 * 2,     -3150 * 2,     -5482 * 2,
      12107 * 2,      -782 * 2,     -5986 * 2,
      13013 * 2,       511 * 2,     -5164 * 2,
      14172 * 2,       840 * 2,     -4841 * 2,
      14794 * 2,      7792 * 2,     -3028 * 2,
      14596 * 2,      7581 * 2,      -813 * 2,
       8808 * 2,      7531 * 2,      2912 * 2,
      13790 * 2,     10971 * 2,     -1294 * 2 };


/*-----------------------------------------*
 * 3nd codebook:   lsp6 to lsp9            *
 *-----------------------------------------*/

static const Word16 dico3_clsp[TSC_taille_dic3*TSC_dim_dic3] = {

       2258 * 2,     -3027 * 2,    -14154 * 2,    -15275 * 2,
       2617 * 2,     -2421 * 2,    -13423 * 2,    -14666 * 2,
       1491 * 2,     -4636 * 2,    -13353 * 2,    -15066 * 2,
       4111 * 2,     -2037 * 2,    -13313 * 2,    -14464 * 2,
       2450 * 2,     -2406 * 2,    -12501 * 2,    -14033 * 2,
       2373 * 2,     -4492 * 2,    -11954 * 2,    -15115 * 2,
       3913 * 2,     -3857 * 2,    -11915 * 2,    -15152 * 2,
       3423 * 2,     -1758 * 2,    -11741 * 2,    -15053 * 2,
       3765 * 2,     -3882 * 2,    -12788 * 2,    -13519 * 2,
       1562 * 2,     -1771 * 2,    -11644 * 2,    -14135 * 2,
       1075 * 2,     -4880 * 2,    -12649 * 2,    -14058 * 2,
       3638 * 2,     -3092 * 2,    -11893 * 2,    -13136 * 2,
       1839 * 2,     -3454 * 2,    -12190 * 2,    -13087 * 2,
        541 * 2,     -3451 * 2,    -11602 * 2,    -13865 * 2,
       1557 * 2,     -4501 * 2,    -11314 * 2,    -12939 * 2,
        872 * 2,     -5493 * 2,    -11547 * 2,    -14155 * 2,
       -297 * 2,     -5352 * 2,    -12385 * 2,    -13504 * 2,
        123 * 2,     -5666 * 2,    -11671 * 2,    -12907 * 2,
        618 * 2,     -5376 * 2,    -10600 * 2,    -13291 * 2,
      -1370 * 2,     -3979 * 2,    -11880 * 2,    -13663 * 2,
      -1024 * 2,     -5949 * 2,    -10829 * 2,    -13536 * 2,
      -1721 * 2,     -5448 * 2,    -11365 * 2,    -13776 * 2,
       -641 * 2,     -6511 * 2,    -11878 * 2,    -14002 * 2,
      -1514 * 2,     -6290 * 2,    -12143 * 2,    -13142 * 2,
      -1583 * 2,     -4767 * 2,    -11486 * 2,    -12786 * 2,
      -1183 * 2,     -5072 * 2,    -10358 * 2,    -13940 * 2,
      -2687 * 2,     -4541 * 2,    -10928 * 2,    -13507 * 2,
      -2919 * 2,     -5785 * 2,    -11747 * 2,    -13494 * 2,
      -2538 * 2,     -6619 * 2,    -11027 * 2,    -13192 * 2,
      -2446 * 2,     -5705 * 2,    -11156 * 2,    -12698 * 2,
      -2043 * 2,     -5217 * 2,    -10126 * 2,    -13152 * 2,
      -3314 * 2,     -5992 * 2,    -10494 * 2,    -13277 * 2,
      -3380 * 2,     -4989 * 2,    -10032 * 2,    -13357 * 2,
      -4071 * 2,     -5683 * 2,    -11128 * 2,    -13624 * 2,
      -3464 * 2,     -4918 * 2,    -11331 * 2,    -13973 * 2,
      -3331 * 2,     -6676 * 2,    -11238 * 2,    -13913 * 2,
      -3686 * 2,     -6545 * 2,    -11198 * 2,    -12740 * 2,
      -4633 * 2,     -5939 * 2,    -10395 * 2,    -13663 * 2,
      -4443 * 2,     -6846 * 2,    -11017 * 2,    -13196 * 2,
      -4066 * 2,     -6849 * 2,    -10301 * 2,    -13913 * 2,
      -4291 * 2,     -6407 * 2,    -11398 * 2,    -14119 * 2,
      -5002 * 2,     -6434 * 2,    -11583 * 2,    -13789 * 2,
      -5310 * 2,     -6544 * 2,    -10539 * 2,    -13746 * 2,
      -5112 * 2,     -7394 * 2,    -11164 * 2,    -13567 * 2,
      -4827 * 2,     -6971 * 2,    -10001 * 2,    -13369 * 2,
      -5657 * 2,     -7140 * 2,    -10365 * 2,    -13493 * 2,
      -5190 * 2,     -8027 * 2,    -10599 * 2,    -13415 * 2,
      -5907 * 2,     -7080 * 2,    -11187 * 2,    -14166 * 2,
      -6284 * 2,     -7633 * 2,    -10747 * 2,    -13860 * 2,
      -6051 * 2,     -7745 * 2,    -11537 * 2,    -13737 * 2,
      -5341 * 2,     -8328 * 2,    -11248 * 2,    -13067 * 2,
      -6310 * 2,     -8473 * 2,    -11362 * 2,    -13596 * 2,
      -6467 * 2,     -8073 * 2,    -10692 * 2,    -12970 * 2,
      -7180 * 2,     -8492 * 2,    -10931 * 2,    -13343 * 2,
      -6977 * 2,     -8046 * 2,    -11613 * 2,    -14006 * 2,
      -5705 * 2,     -7955 * 2,    -12189 * 2,    -13287 * 2,
      -6190 * 2,     -8455 * 2,    -11824 * 2,    -12820 * 2,
      -6802 * 2,     -8892 * 2,    -12057 * 2,    -13555 * 2,
      -5609 * 2,     -9078 * 2,    -11130 * 2,    -13067 * 2,
      -7282 * 2,     -9224 * 2,    -11146 * 2,    -13507 * 2,
      -5587 * 2,     -8876 * 2,    -12187 * 2,    -13602 * 2,
      -6157 * 2,     -9509 * 2,    -12117 * 2,    -13104 * 2,
      -5840 * 2,     -9545 * 2,    -11672 * 2,    -13785 * 2,
      -5068 * 2,     -8786 * 2,    -11683 * 2,    -13401 * 2,
      -6344 * 2,     -9137 * 2,    -12671 * 2,    -13555 * 2,
      -6655 * 2,    -10179 * 2,    -11980 * 2,    -13370 * 2,
      -6242 * 2,    -10092 * 2,    -12522 * 2,    -13432 * 2,
      -5821 * 2,     -9271 * 2,    -12552 * 2,    -14012 * 2,
      -7085 * 2,     -9616 * 2,    -12276 * 2,    -13968 * 2,
      -7061 * 2,     -9362 * 2,    -12863 * 2,    -13880 * 2,
      -6592 * 2,     -8328 * 2,    -12460 * 2,    -13845 * 2,
      -6286 * 2,     -8804 * 2,    -12945 * 2,    -14148 * 2,
      -6710 * 2,    -10050 * 2,    -13020 * 2,    -14176 * 2,
      -6713 * 2,    -10374 * 2,    -12927 * 2,    -13653 * 2,
      -6869 * 2,     -8974 * 2,    -13313 * 2,    -14071 * 2,
      -6771 * 2,     -9578 * 2,    -13010 * 2,    -14529 * 2,
      -7836 * 2,     -9434 * 2,    -13091 * 2,    -14286 * 2,
      -7582 * 2,    -10082 * 2,    -13311 * 2,    -14192 * 2,
      -7685 * 2,     -8768 * 2,    -12554 * 2,    -14190 * 2,
      -6660 * 2,     -9185 * 2,    -13614 * 2,    -14547 * 2,
      -7287 * 2,     -9994 * 2,    -13761 * 2,    -14424 * 2,
      -6119 * 2,    -10003 * 2,    -13519 * 2,    -14326 * 2,
      -7610 * 2,     -8870 * 2,    -13578 * 2,    -14780 * 2,
      -7944 * 2,     -9765 * 2,    -13469 * 2,    -14817 * 2,
      -7068 * 2,    -10465 * 2,    -13490 * 2,    -14913 * 2,
      -6437 * 2,     -9902 * 2,    -13996 * 2,    -14643 * 2,
      -6860 * 2,     -9320 * 2,    -14192 * 2,    -14997 * 2,
      -6033 * 2,     -9723 * 2,    -13600 * 2,    -15112 * 2,
      -7312 * 2,    -10482 * 2,    -14253 * 2,    -15107 * 2,
      -6026 * 2,    -10794 * 2,    -13901 * 2,    -14777 * 2,
      -7911 * 2,    -10718 * 2,    -14030 * 2,    -14709 * 2,
      -6908 * 2,    -11331 * 2,    -13939 * 2,    -15039 * 2,
      -5948 * 2,    -10946 * 2,    -13395 * 2,    -14592 * 2,
      -5903 * 2,     -9715 * 2,    -14491 * 2,    -15166 * 2,
      -7071 * 2,    -11312 * 2,    -14446 * 2,    -15360 * 2,
      -6732 * 2,    -10351 * 2,    -13546 * 2,    -15582 * 2,
      -7547 * 2,    -11517 * 2,    -13733 * 2,    -14499 * 2,
      -8135 * 2,    -11430 * 2,    -14394 * 2,    -14974 * 2,
      -8436 * 2,    -11182 * 2,    -13843 * 2,    -14806 * 2,
      -7833 * 2,    -12031 * 2,    -14002 * 2,    -15207 * 2,
      -8091 * 2,    -10713 * 2,    -13635 * 2,    -14465 * 2,
      -8509 * 2,    -11642 * 2,    -13383 * 2,    -15333 * 2,
      -9007 * 2,    -11657 * 2,    -14100 * 2,    -14737 * 2,
      -8815 * 2,    -11159 * 2,    -13312 * 2,    -14700 * 2,
      -9207 * 2,    -11913 * 2,    -13759 * 2,    -14560 * 2,
      -9285 * 2,    -10831 * 2,    -13848 * 2,    -14694 * 2,
      -8597 * 2,    -11303 * 2,    -13516 * 2,    -14163 * 2,
      -9264 * 2,    -11873 * 2,    -14171 * 2,    -15153 * 2,
     -10039 * 2,    -11402 * 2,    -13745 * 2,    -14729 * 2,
      -9459 * 2,    -11244 * 2,    -14384 * 2,    -14989 * 2,
      -9566 * 2,    -10760 * 2,    -13370 * 2,    -14409 * 2,
      -8917 * 2,    -10214 * 2,    -13810 * 2,    -14916 * 2,
      -9932 * 2,    -11413 * 2,    -13003 * 2,    -14482 * 2,
      -8704 * 2,    -10367 * 2,    -13309 * 2,    -14295 * 2,
     -10370 * 2,    -11827 * 2,    -13525 * 2,    -14295 * 2,
      -9682 * 2,    -11721 * 2,    -13258 * 2,    -13944 * 2,
      -8606 * 2,    -11043 * 2,    -13115 * 2,    -13931 * 2,
      -9323 * 2,    -10436 * 2,    -12611 * 2,    -14175 * 2,
      -8799 * 2,    -10894 * 2,    -12675 * 2,    -13938 * 2,
      -8745 * 2,    -11134 * 2,    -12565 * 2,    -14705 * 2,
      -8699 * 2,    -11343 * 2,    -12686 * 2,    -13532 * 2,
      -8615 * 2,    -10184 * 2,    -12313 * 2,    -13908 * 2,
      -8013 * 2,    -10163 * 2,    -12923 * 2,    -13918 * 2,
      -7827 * 2,    -10632 * 2,    -12342 * 2,    -13684 * 2,
      -8801 * 2,    -10562 * 2,    -12012 * 2,    -13314 * 2,
      -7782 * 2,     -9941 * 2,    -12520 * 2,    -13508 * 2,
      -7184 * 2,    -10374 * 2,    -12445 * 2,    -14056 * 2,
      -8329 * 2,     -9450 * 2,    -12339 * 2,    -14009 * 2,
      -8129 * 2,     -9827 * 2,    -11637 * 2,    -13458 * 2,
      -7572 * 2,     -9148 * 2,    -12059 * 2,    -13390 * 2,
      -8677 * 2,     -9708 * 2,    -12866 * 2,    -14477 * 2,
      -7594 * 2,    -10662 * 2,    -12969 * 2,    -14495 * 2,
      -6901 * 2,     -8991 * 2,    -12128 * 2,    -14424 * 2,
      -7217 * 2,     -8242 * 2,    -12871 * 2,    -14508 * 2,
      -7034 * 2,    -10885 * 2,    -13420 * 2,    -14158 * 2,
      -6348 * 2,     -9423 * 2,    -12826 * 2,    -14969 * 2,
      -6111 * 2,    -10844 * 2,    -12297 * 2,    -14546 * 2,
      -5493 * 2,    -10258 * 2,    -12720 * 2,    -14164 * 2,
      -5372 * 2,     -9350 * 2,    -12931 * 2,    -14455 * 2,
      -5654 * 2,     -9935 * 2,    -12744 * 2,    -15081 * 2,
      -5464 * 2,     -8605 * 2,    -13009 * 2,    -14723 * 2,
      -5012 * 2,    -10523 * 2,    -13054 * 2,    -14874 * 2,
      -4898 * 2,     -9333 * 2,    -13567 * 2,    -14618 * 2,
      -5347 * 2,     -9040 * 2,    -12207 * 2,    -15051 * 2,
      -4353 * 2,     -9886 * 2,    -12281 * 2,    -14405 * 2,
      -4458 * 2,     -8878 * 2,    -12945 * 2,    -14333 * 2,
      -4073 * 2,     -8797 * 2,    -12507 * 2,    -14647 * 2,
      -3601 * 2,     -8943 * 2,    -13112 * 2,    -14775 * 2,
      -4518 * 2,     -9080 * 2,    -13530 * 2,    -15080 * 2,
      -4159 * 2,     -9772 * 2,    -12836 * 2,    -14063 * 2,
      -4909 * 2,     -8973 * 2,    -12211 * 2,    -14094 * 2,
      -4961 * 2,     -8262 * 2,    -12126 * 2,    -14627 * 2,
      -4409 * 2,     -8219 * 2,    -12531 * 2,    -14028 * 2,
      -4435 * 2,     -7684 * 2,    -12834 * 2,    -14629 * 2,
      -3529 * 2,     -7962 * 2,    -12876 * 2,    -14255 * 2,
      -3751 * 2,     -8655 * 2,    -13494 * 2,    -14204 * 2,
      -5246 * 2,     -8294 * 2,    -13381 * 2,    -14163 * 2,
      -5690 * 2,     -8131 * 2,    -12666 * 2,    -14052 * 2,
      -4438 * 2,     -7338 * 2,    -13222 * 2,    -13975 * 2,
      -5584 * 2,     -7365 * 2,    -13205 * 2,    -14332 * 2,
      -4946 * 2,     -7365 * 2,    -12549 * 2,    -13708 * 2,
      -3934 * 2,     -7160 * 2,    -12882 * 2,    -13679 * 2,
      -5238 * 2,     -8579 * 2,    -12846 * 2,    -13617 * 2,
      -4910 * 2,     -6631 * 2,    -12587 * 2,    -14270 * 2,
      -5542 * 2,     -7117 * 2,    -12123 * 2,    -13862 * 2,
      -4012 * 2,     -7223 * 2,    -12212 * 2,    -13973 * 2,
      -5125 * 2,     -8049 * 2,    -11943 * 2,    -13886 * 2,
      -3879 * 2,     -7351 * 2,    -12375 * 2,    -13382 * 2,
      -4600 * 2,     -6946 * 2,    -12021 * 2,    -13168 * 2,
      -4605 * 2,     -7897 * 2,    -11739 * 2,    -13337 * 2,
      -4082 * 2,     -8445 * 2,    -12068 * 2,    -13597 * 2,
      -3665 * 2,     -7463 * 2,    -11557 * 2,    -13688 * 2,
      -3867 * 2,     -7581 * 2,    -11625 * 2,    -12775 * 2,
      -3722 * 2,     -7888 * 2,    -10925 * 2,    -13470 * 2,
      -3031 * 2,     -8207 * 2,    -12287 * 2,    -13101 * 2,
      -2639 * 2,     -7377 * 2,    -11750 * 2,    -13247 * 2,
      -3291 * 2,     -8811 * 2,    -11872 * 2,    -12968 * 2,
      -3759 * 2,     -8500 * 2,    -10902 * 2,    -12955 * 2,
      -4233 * 2,     -8614 * 2,    -11788 * 2,    -12483 * 2,
      -3832 * 2,     -8933 * 2,    -11332 * 2,    -12544 * 2,
      -2269 * 2,     -8265 * 2,    -11190 * 2,    -12700 * 2,
      -2151 * 2,     -8973 * 2,    -11367 * 2,    -13237 * 2,
      -3700 * 2,     -8547 * 2,    -10965 * 2,    -12117 * 2,
      -2425 * 2,     -8480 * 2,    -10385 * 2,    -13061 * 2,
      -4135 * 2,     -7675 * 2,    -10670 * 2,    -12717 * 2,
      -2463 * 2,     -7452 * 2,    -10692 * 2,    -12495 * 2,
      -4116 * 2,     -7373 * 2,    -11309 * 2,    -12239 * 2,
      -3509 * 2,     -7888 * 2,    -10270 * 2,    -11879 * 2,
      -3150 * 2,     -7234 * 2,    -10169 * 2,    -13063 * 2,
      -3005 * 2,     -7089 * 2,    -10930 * 2,    -11839 * 2,
      -3384 * 2,     -7758 * 2,     -9501 * 2,    -12527 * 2,
      -4448 * 2,     -7776 * 2,    -10006 * 2,    -12915 * 2,
      -4547 * 2,     -6600 * 2,    -10270 * 2,    -12426 * 2,
      -3027 * 2,     -6337 * 2,    -10287 * 2,    -12059 * 2,
      -3216 * 2,     -6586 * 2,     -9655 * 2,    -12890 * 2,
      -4121 * 2,     -6902 * 2,     -9225 * 2,    -12883 * 2,
      -2915 * 2,     -5839 * 2,     -9437 * 2,    -12387 * 2,
      -4060 * 2,     -5770 * 2,     -9461 * 2,    -13002 * 2,
      -2924 * 2,     -7160 * 2,     -8755 * 2,    -12284 * 2,
      -5042 * 2,     -6390 * 2,     -9057 * 2,    -13100 * 2,
      -4201 * 2,     -5887 * 2,     -8410 * 2,    -13116 * 2,
      -2889 * 2,     -6473 * 2,     -8294 * 2,    -13110 * 2,
      -3461 * 2,     -5187 * 2,     -8362 * 2,    -12949 * 2,
      -4651 * 2,     -6493 * 2,     -8009 * 2,    -12209 * 2,
      -5323 * 2,     -7101 * 2,     -9166 * 2,    -12660 * 2,
      -3232 * 2,     -5772 * 2,     -7693 * 2,    -11727 * 2,
      -3954 * 2,     -6435 * 2,     -8949 * 2,    -11244 * 2,
      -2691 * 2,     -6185 * 2,     -8260 * 2,    -11380 * 2,
      -2648 * 2,     -5363 * 2,     -8576 * 2,    -11475 * 2,
      -5102 * 2,     -6875 * 2,     -8511 * 2,    -11190 * 2,
      -3824 * 2,     -6717 * 2,     -9987 * 2,    -11387 * 2,
      -2425 * 2,     -6443 * 2,     -9186 * 2,    -10955 * 2,
      -2216 * 2,     -7118 * 2,     -9656 * 2,    -11662 * 2,
      -2027 * 2,     -5643 * 2,     -9985 * 2,    -11669 * 2,
      -1582 * 2,     -6236 * 2,     -9060 * 2,    -12331 * 2,
      -2068 * 2,     -7313 * 2,    -10351 * 2,    -11582 * 2,
       -755 * 2,     -6844 * 2,     -9710 * 2,    -11707 * 2,
      -1518 * 2,     -7717 * 2,    -10100 * 2,    -12380 * 2,
       -744 * 2,     -6089 * 2,    -10481 * 2,    -11823 * 2,
      -1635 * 2,     -5886 * 2,    -10856 * 2,    -12158 * 2,
       -404 * 2,     -7526 * 2,    -10553 * 2,    -12547 * 2,
      -1336 * 2,     -7077 * 2,    -11289 * 2,    -12512 * 2,
        500 * 2,     -6955 * 2,    -10065 * 2,    -12438 * 2,
        -17 * 2,     -5673 * 2,    -10843 * 2,    -12286 * 2,
       -815 * 2,     -6924 * 2,     -9235 * 2,    -12999 * 2,
         88 * 2,     -7143 * 2,    -11346 * 2,    -13241 * 2,
      -1661 * 2,     -7418 * 2,    -10376 * 2,    -13501 * 2,
      -1128 * 2,     -8022 * 2,    -11858 * 2,    -12875 * 2,
      -1097 * 2,     -8331 * 2,    -11141 * 2,    -13565 * 2,
      -1839 * 2,     -7119 * 2,    -11451 * 2,    -13819 * 2,
       -331 * 2,     -7501 * 2,     -9970 * 2,    -13978 * 2,
      -2238 * 2,     -6683 * 2,    -10486 * 2,    -13947 * 2,
       -583 * 2,     -7350 * 2,    -11045 * 2,    -14570 * 2,
      -2354 * 2,     -8320 * 2,    -11714 * 2,    -13950 * 2,
      -1573 * 2,     -8210 * 2,    -12170 * 2,    -13626 * 2,
      -2298 * 2,     -7840 * 2,    -11162 * 2,    -14578 * 2,
      -2478 * 2,     -6535 * 2,    -11432 * 2,    -14457 * 2,
      -2071 * 2,     -6924 * 2,    -12417 * 2,    -14293 * 2,
      -1423 * 2,     -6187 * 2,    -11723 * 2,    -14768 * 2,
      -2987 * 2,     -6901 * 2,    -12229 * 2,    -13837 * 2,
      -2254 * 2,     -8289 * 2,    -12343 * 2,    -14663 * 2,
      -3354 * 2,     -7942 * 2,    -12155 * 2,    -14284 * 2,
      -3437 * 2,     -6531 * 2,    -12137 * 2,    -14527 * 2,
      -2882 * 2,     -7340 * 2,    -11837 * 2,    -15079 * 2,
      -3021 * 2,     -7084 * 2,    -12860 * 2,    -14905 * 2,
      -2649 * 2,     -5888 * 2,    -12910 * 2,    -14579 * 2,
      -1922 * 2,     -6363 * 2,    -12567 * 2,    -15093 * 2,
      -3917 * 2,     -7895 * 2,    -12288 * 2,    -15009 * 2,
      -4268 * 2,     -7384 * 2,    -11802 * 2,    -14440 * 2,
      -3114 * 2,     -5648 * 2,    -12388 * 2,    -15124 * 2,
      -4192 * 2,     -6273 * 2,    -12999 * 2,    -14330 * 2,
      -3604 * 2,     -7229 * 2,    -13430 * 2,    -14509 * 2,
      -4506 * 2,     -5775 * 2,    -12538 * 2,    -14447 * 2,
      -4797 * 2,     -6306 * 2,    -13096 * 2,    -14957 * 2,
      -3669 * 2,     -5302 * 2,    -12349 * 2,    -14277 * 2,
      -3340 * 2,     -5796 * 2,    -13124 * 2,    -14068 * 2,
      -3591 * 2,     -5380 * 2,    -13440 * 2,    -14707 * 2,
      -2836 * 2,     -6184 * 2,    -13714 * 2,    -14435 * 2,
      -4821 * 2,     -6920 * 2,    -13606 * 2,    -14505 * 2,
      -3866 * 2,     -6360 * 2,    -14015 * 2,    -14772 * 2,
      -4129 * 2,     -6015 * 2,    -12428 * 2,    -13774 * 2,
      -2525 * 2,     -7214 * 2,    -13270 * 2,    -14143 * 2,
      -2688 * 2,     -5151 * 2,    -12312 * 2,    -14060 * 2,
      -2841 * 2,     -4505 * 2,    -12737 * 2,    -14704 * 2,
      -2167 * 2,     -4432 * 2,    -13236 * 2,    -14351 * 2,
      -1674 * 2,     -5285 * 2,    -13542 * 2,    -14729 * 2,
      -1682 * 2,     -5312 * 2,    -12819 * 2,    -13912 * 2,
      -2689 * 2,     -5766 * 2,    -12461 * 2,    -13395 * 2,
       -522 * 2,     -5970 * 2,    -12628 * 2,    -14436 * 2,
        -37 * 2,     -6079 * 2,    -13267 * 2,    -14165 * 2,
      -1251 * 2,     -3344 * 2,    -12940 * 2,    -14242 * 2,
      -2042 * 2,     -3653 * 2,    -12979 * 2,    -14880 * 2,
      -1563 * 2,     -4470 * 2,    -11756 * 2,    -14774 * 2,
       -252 * 2,     -4094 * 2,    -12046 * 2,    -14881 * 2,
      -2590 * 2,     -4088 * 2,    -11559 * 2,    -14093 * 2,
      -1817 * 2,     -3406 * 2,    -11146 * 2,    -14214 * 2,
      -2347 * 2,     -4840 * 2,    -10723 * 2,    -14719 * 2,
       -333 * 2,     -4997 * 2,    -10524 * 2,    -14779 * 2,
      -1048 * 2,     -3769 * 2,    -10195 * 2,    -14339 * 2,
      -2502 * 2,     -5724 * 2,    -10197 * 2,    -14570 * 2,
      -2069 * 2,     -4009 * 2,     -9460 * 2,    -14313 * 2,
      -3075 * 2,     -4500 * 2,     -9738 * 2,    -14703 * 2,
      -3665 * 2,     -5650 * 2,    -10968 * 2,    -14696 * 2,
      -3038 * 2,     -5592 * 2,     -9134 * 2,    -14096 * 2,
      -1683 * 2,     -5635 * 2,     -8743 * 2,    -14555 * 2,
      -3478 * 2,     -6306 * 2,     -9632 * 2,    -14834 * 2,
      -2310 * 2,     -6941 * 2,     -9412 * 2,    -14113 * 2,
      -1789 * 2,     -7180 * 2,    -10026 * 2,    -14917 * 2,
      -3649 * 2,     -7174 * 2,    -10798 * 2,    -14964 * 2,
      -4884 * 2,     -6782 * 2,     -9331 * 2,    -14539 * 2,
      -4130 * 2,     -7783 * 2,     -9704 * 2,    -14296 * 2,
      -3209 * 2,     -8267 * 2,    -10448 * 2,    -14139 * 2,
      -1523 * 2,     -8503 * 2,    -10050 * 2,    -14744 * 2,
      -4967 * 2,     -7233 * 2,    -10980 * 2,    -14622 * 2,
      -4405 * 2,     -8268 * 2,    -11101 * 2,    -14164 * 2,
      -4444 * 2,     -7988 * 2,    -11414 * 2,    -15169 * 2,
      -5731 * 2,     -7986 * 2,    -10506 * 2,    -14819 * 2,
      -3533 * 2,     -8866 * 2,    -11731 * 2,    -14582 * 2,
      -3904 * 2,     -9181 * 2,    -11170 * 2,    -13855 * 2,
      -2560 * 2,     -9499 * 2,    -11078 * 2,    -14844 * 2,
      -5102 * 2,     -9641 * 2,    -11597 * 2,    -14571 * 2,
      -5762 * 2,     -8671 * 2,    -11671 * 2,    -14272 * 2,
      -5555 * 2,     -7699 * 2,    -12183 * 2,    -14734 * 2,
      -3257 * 2,     -8774 * 2,    -12130 * 2,    -15150 * 2,
      -6385 * 2,     -8877 * 2,    -10879 * 2,    -14952 * 2,
      -2639 * 2,     -9484 * 2,    -12495 * 2,    -14502 * 2,
      -4573 * 2,    -10038 * 2,    -12044 * 2,    -13615 * 2,
      -2839 * 2,     -9657 * 2,    -11759 * 2,    -13773 * 2,
      -4563 * 2,     -8769 * 2,    -12475 * 2,    -13275 * 2,
      -3931 * 2,     -9965 * 2,    -12686 * 2,    -13508 * 2,
      -2446 * 2,     -9300 * 2,    -12342 * 2,    -13708 * 2,
      -2678 * 2,     -8421 * 2,    -12910 * 2,    -13811 * 2,
      -4683 * 2,     -9297 * 2,    -12269 * 2,    -12973 * 2,
      -2442 * 2,    -10261 * 2,    -13015 * 2,    -14122 * 2,
      -4892 * 2,     -9479 * 2,    -11754 * 2,    -12898 * 2,
      -5478 * 2,     -9556 * 2,    -13200 * 2,    -13844 * 2,
      -4707 * 2,    -10400 * 2,    -13392 * 2,    -14169 * 2,
      -3466 * 2,     -9494 * 2,    -13865 * 2,    -14327 * 2,
      -3353 * 2,     -8233 * 2,    -13869 * 2,    -14575 * 2,
      -5703 * 2,     -8553 * 2,    -13673 * 2,    -14414 * 2,
      -4650 * 2,     -9988 * 2,    -14209 * 2,    -14789 * 2,
      -5839 * 2,     -8786 * 2,    -13931 * 2,    -14949 * 2,
      -5016 * 2,     -7751 * 2,    -13859 * 2,    -14823 * 2,
      -3606 * 2,     -8092 * 2,    -14241 * 2,    -14949 * 2,
      -4620 * 2,     -7842 * 2,    -14347 * 2,    -15262 * 2,
      -3504 * 2,     -7468 * 2,    -13683 * 2,    -15120 * 2,
      -4662 * 2,     -8042 * 2,    -12990 * 2,    -15281 * 2,
      -5991 * 2,     -7721 * 2,    -13275 * 2,    -15019 * 2,
      -5360 * 2,     -7044 * 2,    -13474 * 2,    -15424 * 2,
      -4058 * 2,     -8833 * 2,    -12783 * 2,    -15391 * 2,
      -3911 * 2,     -8840 * 2,    -13626 * 2,    -15721 * 2,
      -4785 * 2,     -9441 * 2,    -13014 * 2,    -15382 * 2,
      -5789 * 2,     -8844 * 2,    -12818 * 2,    -15475 * 2,
      -3056 * 2,     -8107 * 2,    -12624 * 2,    -15369 * 2,
      -3757 * 2,     -6875 * 2,    -12637 * 2,    -15623 * 2,
      -2557 * 2,     -8635 * 2,    -13284 * 2,    -15087 * 2,
      -3702 * 2,     -9042 * 2,    -12656 * 2,    -15986 * 2,
      -4066 * 2,    -10051 * 2,    -12271 * 2,    -15442 * 2,
      -2533 * 2,    -10409 * 2,    -13255 * 2,    -15462 * 2,
      -5357 * 2,    -10545 * 2,    -13401 * 2,    -15642 * 2,
      -4009 * 2,    -10553 * 2,    -13868 * 2,    -15118 * 2,
      -3004 * 2,     -9401 * 2,    -14284 * 2,    -15227 * 2,
      -3917 * 2,    -10672 * 2,    -14563 * 2,    -15613 * 2,
      -1729 * 2,     -9887 * 2,    -13941 * 2,    -14722 * 2,
      -5558 * 2,    -10041 * 2,    -14600 * 2,    -15447 * 2,
      -4949 * 2,    -11441 * 2,    -14321 * 2,    -14933 * 2,
      -3688 * 2,    -11414 * 2,    -13850 * 2,    -14526 * 2,
      -5148 * 2,    -12634 * 2,    -14279 * 2,    -15609 * 2,
      -5913 * 2,    -11034 * 2,    -14790 * 2,    -15975 * 2,
      -6840 * 2,    -12206 * 2,    -14079 * 2,    -15646 * 2,
      -7445 * 2,    -11687 * 2,    -14823 * 2,    -15405 * 2,
      -7657 * 2,    -11124 * 2,    -14105 * 2,    -15748 * 2,
      -7460 * 2,    -11012 * 2,    -13241 * 2,    -15590 * 2,
      -7816 * 2,    -11903 * 2,    -15113 * 2,    -15635 * 2,
      -7283 * 2,    -12023 * 2,    -14984 * 2,    -16015 * 2,
      -8586 * 2,    -12702 * 2,    -14420 * 2,    -15665 * 2,
      -8768 * 2,    -11668 * 2,    -14871 * 2,    -16066 * 2,
      -8938 * 2,    -11004 * 2,    -14763 * 2,    -15378 * 2,
      -8517 * 2,    -10255 * 2,    -14626 * 2,    -15736 * 2,
      -9509 * 2,    -10924 * 2,    -14371 * 2,    -15669 * 2,
      -9560 * 2,    -12148 * 2,    -14681 * 2,    -15573 * 2,
      -9524 * 2,    -11637 * 2,    -13931 * 2,    -15687 * 2,
      -9952 * 2,    -11446 * 2,    -14749 * 2,    -16022 * 2,
     -10589 * 2,    -11733 * 2,    -14404 * 2,    -15410 * 2,
     -10599 * 2,    -12143 * 2,    -14897 * 2,    -15877 * 2,
      -9402 * 2,    -12839 * 2,    -14943 * 2,    -15516 * 2,
      -8972 * 2,    -12682 * 2,    -14666 * 2,    -15233 * 2,
     -10771 * 2,    -13173 * 2,    -14355 * 2,    -15621 * 2,
     -10297 * 2,    -13564 * 2,    -14733 * 2,    -15734 * 2,
     -10822 * 2,    -12937 * 2,    -14477 * 2,    -15042 * 2,
     -11142 * 2,    -12785 * 2,    -15139 * 2,    -15905 * 2,
     -10101 * 2,    -12649 * 2,    -14283 * 2,    -14883 * 2,
     -10035 * 2,    -12902 * 2,    -13731 * 2,    -15264 * 2,
     -11322 * 2,    -13805 * 2,    -15031 * 2,    -15382 * 2,
     -12119 * 2,    -13184 * 2,    -14533 * 2,    -15272 * 2,
     -10453 * 2,    -12308 * 2,    -14055 * 2,    -14623 * 2,
     -10909 * 2,    -12027 * 2,    -13531 * 2,    -15149 * 2,
     -10569 * 2,    -12502 * 2,    -13783 * 2,    -14390 * 2,
     -12295 * 2,    -13452 * 2,    -14519 * 2,    -16014 * 2,
      -9070 * 2,    -13916 * 2,    -15100 * 2,    -15778 * 2,
      -9173 * 2,    -13031 * 2,    -15055 * 2,    -16027 * 2,
     -10759 * 2,    -13738 * 2,    -15398 * 2,    -16195 * 2,
     -10454 * 2,    -11762 * 2,    -15233 * 2,    -16185 * 2,
      -8931 * 2,    -13657 * 2,    -15532 * 2,    -15997 * 2,
     -10846 * 2,    -14286 * 2,    -15519 * 2,    -16176 * 2,
      -7775 * 2,    -10326 * 2,    -15224 * 2,    -15794 * 2,
      -8874 * 2,    -10860 * 2,    -13493 * 2,    -15558 * 2,
      -7842 * 2,    -10236 * 2,    -14047 * 2,    -15982 * 2,
      -8241 * 2,     -9518 * 2,    -14118 * 2,    -15253 * 2,
      -7417 * 2,     -8839 * 2,    -14502 * 2,    -15988 * 2,
      -6511 * 2,     -8764 * 2,    -14745 * 2,    -15662 * 2,
      -6592 * 2,     -9218 * 2,    -13452 * 2,    -15816 * 2,
      -8666 * 2,    -10330 * 2,    -12824 * 2,    -15364 * 2,
      -6878 * 2,     -8214 * 2,    -14008 * 2,    -15242 * 2,
      -6472 * 2,     -9944 * 2,    -12580 * 2,    -15630 * 2,
      -6093 * 2,     -7786 * 2,    -14146 * 2,    -15705 * 2,
      -7140 * 2,     -9175 * 2,    -12134 * 2,    -15479 * 2,
      -7465 * 2,    -10291 * 2,    -12240 * 2,    -15235 * 2,
      -6556 * 2,     -8071 * 2,    -13251 * 2,    -14474 * 2,
      -6529 * 2,     -8031 * 2,    -12129 * 2,    -14681 * 2,
      -5876 * 2,    -10227 * 2,    -11466 * 2,    -15307 * 2,
      -5654 * 2,     -6943 * 2,    -12564 * 2,    -14736 * 2,
      -6353 * 2,     -7469 * 2,    -12201 * 2,    -14350 * 2,
      -5627 * 2,     -7728 * 2,    -11689 * 2,    -15788 * 2,
      -5197 * 2,     -6409 * 2,    -11655 * 2,    -14597 * 2,
      -4239 * 2,     -6118 * 2,    -11432 * 2,    -15123 * 2,
      -4436 * 2,     -9115 * 2,    -10839 * 2,    -15679 * 2,
      -3548 * 2,     -8408 * 2,    -10345 * 2,    -15621 * 2,
      -4227 * 2,     -7315 * 2,     -9749 * 2,    -15819 * 2,
      -3051 * 2,     -6359 * 2,    -10763 * 2,    -15801 * 2,
      -1744 * 2,     -8688 * 2,    -11352 * 2,    -15519 * 2,
      -1507 * 2,     -7789 * 2,    -12284 * 2,    -15213 * 2,
       -898 * 2,     -8081 * 2,    -12437 * 2,    -15816 * 2,
      -1027 * 2,    -10079 * 2,    -11965 * 2,    -15213 * 2,
       -734 * 2,     -7974 * 2,    -12549 * 2,    -14605 * 2,
        827 * 2,     -7472 * 2,    -12192 * 2,    -15324 * 2,
       -511 * 2,     -9088 * 2,    -11661 * 2,    -14363 * 2,
      -1296 * 2,     -8962 * 2,    -13040 * 2,    -14289 * 2,
      -1798 * 2,     -7790 * 2,    -13720 * 2,    -14889 * 2,
      -1027 * 2,     -7855 * 2,    -13842 * 2,    -15336 * 2,
       -438 * 2,     -6447 * 2,    -13350 * 2,    -15000 * 2,
      -1645 * 2,     -7871 * 2,    -13618 * 2,    -14276 * 2,
      -1796 * 2,     -7857 * 2,    -14187 * 2,    -14647 * 2,
       1419 * 2,     -7430 * 2,    -13326 * 2,    -14604 * 2,
        968 * 2,     -8714 * 2,    -14064 * 2,    -14746 * 2,
      -1389 * 2,     -8443 * 2,    -14591 * 2,    -15035 * 2,
      -2376 * 2,     -6547 * 2,    -14289 * 2,    -15131 * 2,
        681 * 2,     -6039 * 2,    -14174 * 2,    -15261 * 2,
      -1278 * 2,     -7426 * 2,    -14995 * 2,    -15460 * 2,
      -1650 * 2,     -5306 * 2,    -14072 * 2,    -15417 * 2,
       1179 * 2,     -5908 * 2,    -12475 * 2,    -14644 * 2,
        156 * 2,     -7873 * 2,    -12697 * 2,    -13694 * 2,
       -624 * 2,     -4770 * 2,    -12873 * 2,    -15838 * 2,
      -3129 * 2,     -4907 * 2,    -13467 * 2,    -15778 * 2,
      -4247 * 2,     -5883 * 2,    -14196 * 2,    -15500 * 2,
      -4086 * 2,     -7081 * 2,    -14616 * 2,    -15684 * 2,
      -2201 * 2,     -7651 * 2,    -12663 * 2,    -13344 * 2,
      -3453 * 2,     -6630 * 2,    -12004 * 2,    -12864 * 2,
        856 * 2,     -8323 * 2,    -11736 * 2,    -13622 * 2,
      -2646 * 2,     -7110 * 2,    -11587 * 2,    -12505 * 2,
      -5336 * 2,     -7532 * 2,    -11484 * 2,    -12670 * 2,
      -5033 * 2,     -8705 * 2,    -10678 * 2,    -12334 * 2,
      -4970 * 2,     -7667 * 2,    -10717 * 2,    -11807 * 2,
      -5795 * 2,     -8993 * 2,    -11290 * 2,    -12255 * 2,
      -4249 * 2,     -8531 * 2,     -9938 * 2,    -12187 * 2,
      -5502 * 2,     -7806 * 2,     -9865 * 2,    -12731 * 2,
      -5921 * 2,     -8823 * 2,    -10290 * 2,    -13093 * 2,
      -6168 * 2,     -8325 * 2,    -10247 * 2,    -11776 * 2,
      -6218 * 2,     -7581 * 2,     -9573 * 2,    -12938 * 2,
      -4793 * 2,     -7741 * 2,     -9545 * 2,    -11334 * 2,
      -6263 * 2,     -9720 * 2,    -11069 * 2,    -12886 * 2,
      -6931 * 2,     -9690 * 2,    -11618 * 2,    -12771 * 2,
      -4180 * 2,     -6051 * 2,     -8002 * 2,    -10238 * 2,
      -2494 * 2,     -5755 * 2,     -7791 * 2,    -13517 * 2,
      -2327 * 2,     -4509 * 2,     -8991 * 2,    -12941 * 2,
      -2882 * 2,     -5563 * 2,     -7443 * 2,    -10727 * 2,
      -1012 * 2,     -5319 * 2,     -9297 * 2,    -12505 * 2,
      -1512 * 2,     -5510 * 2,     -7537 * 2,    -11602 * 2,
      -1865 * 2,     -4284 * 2,     -8937 * 2,    -11243 * 2,
       -675 * 2,     -5828 * 2,     -9165 * 2,    -11035 * 2,
        121 * 2,     -5394 * 2,     -8187 * 2,    -12358 * 2,
      -1828 * 2,     -3798 * 2,     -7161 * 2,    -12607 * 2,
      -1250 * 2,     -3403 * 2,     -8971 * 2,    -12892 * 2,
         82 * 2,     -3228 * 2,     -8165 * 2,    -12383 * 2,
       -205 * 2,     -4263 * 2,     -9511 * 2,    -13198 * 2,
       -970 * 2,     -4326 * 2,    -10376 * 2,    -12376 * 2,
        771 * 2,     -5818 * 2,     -9306 * 2,    -12364 * 2,
       1200 * 2,     -4560 * 2,     -9731 * 2,    -11691 * 2,
        373 * 2,     -3855 * 2,    -10589 * 2,    -11908 * 2,
       -806 * 2,     -2391 * 2,     -9914 * 2,    -12909 * 2,
       -340 * 2,     -5531 * 2,    -10039 * 2,    -11121 * 2,
        -29 * 2,     -4047 * 2,     -9090 * 2,    -10202 * 2,
       -118 * 2,     -6292 * 2,     -9474 * 2,    -13728 * 2,
       2151 * 2,     -6644 * 2,    -10546 * 2,    -13000 * 2,
       3581 * 2,     -3940 * 2,    -10803 * 2,    -12157 * 2,
        281 * 2,     -3753 * 2,     -9785 * 2,    -14179 * 2,
       -970 * 2,     -4782 * 2,     -7680 * 2,    -13979 * 2,
        616 * 2,     -2056 * 2,     -9585 * 2,    -14312 * 2,
        301 * 2,     -1319 * 2,     -7377 * 2,    -12987 * 2,
      -2860 * 2,     -4579 * 2,     -7377 * 2,    -13190 * 2,
       1579 * 2,     -3254 * 2,     -5673 * 2,    -12155 * 2,
       -670 * 2,     -3002 * 2,     -6917 * 2,    -10909 * 2,
       1148 * 2,      -112 * 2,     -6676 * 2,    -13294 * 2,
      -1456 * 2,     -4559 * 2,     -6359 * 2,    -11251 * 2,
       2315 * 2,     -4495 * 2,     -6938 * 2,    -10622 * 2,
       -910 * 2,     -4965 * 2,     -7196 * 2,     -9729 * 2,
      -3313 * 2,     -5078 * 2,     -7349 * 2,     -9257 * 2,
       2995 * 2,        30 * 2,     -6398 * 2,    -14337 * 2,
       2728 * 2,      1352 * 2,     -6451 * 2,    -14372 * 2,
       4050 * 2,     -5770 * 2,    -10681 * 2,    -14211 * 2,
       3555 * 2,      2365 * 2,     -8304 * 2,    -14656 * 2,
       3378 * 2,       470 * 2,    -11318 * 2,    -14972 * 2,
       4295 * 2,       920 * 2,    -11156 * 2,    -14996 * 2,
       4690 * 2,      2206 * 2,    -10216 * 2,    -14946 * 2,
       5240 * 2,       661 * 2,    -11119 * 2,    -14998 * 2,
       6099 * 2,      1695 * 2,     -9999 * 2,    -14947 * 2,
       4890 * 2,     -1618 * 2,    -11589 * 2,    -15054 * 2,
       6603 * 2,     -2072 * 2,    -11463 * 2,    -15067 * 2,
       5215 * 2,      3325 * 2,     -9190 * 2,    -14925 * 2,
       4800 * 2,     -1018 * 2,    -12645 * 2,    -13591 * 2,
       5953 * 2,      3144 * 2,     -7673 * 2,    -14850 * 2,
       4467 * 2,      3391 * 2,     -6908 * 2,    -14827 * 2,
       6908 * 2,      4412 * 2,     -7205 * 2,    -14816 * 2,
       5090 * 2,      3105 * 2,     -5817 * 2,    -14749 * 2,
       5482 * 2,      3606 * 2,     -3761 * 2,    -14629 * 2,
       4139 * 2,      2850 * 2,     -1782 * 2,    -14239 * 2,
       5611 * 2,      4106 * 2,      -669 * 2,    -14497 * 2,
       7300 * 2,      3950 * 2,       736 * 2,    -13497 * 2,
       5376 * 2,     -5185 * 2,    -11826 * 2,    -15194 * 2,
       2924 * 2,     -6522 * 2,    -12052 * 2,    -15281 * 2,
      -6084 * 2,    -11147 * 2,    -12964 * 2,    -14025 * 2,
      -7054 * 2,    -11706 * 2,    -13153 * 2,    -14339 * 2 };

//#include "clsp_334.tab"

/************************************************************************
*
*	FILENAME		:	ener_qua.tab
*
*	DESCRIPTION		:	Energy codebook in Q8.
*
************************************************************************/

#define TSC_nb_qua_ener 64

static const Word16 t_qua_ener[64*2] = {
      12,       48,
      69,      287,
      98,      589,
     103,      865,
     406,      561,
      95,     1161,
     418,      913,
     465,     1288,
     705,     1056,
     822,      699,
    1036,     1016,
     842,     1353,
     555,     1646,
     681,     2119,
     891,     1758,
    1071,     1529,
    1236,     1281,
    1455,      962,
    1661,     1279,
    1416,     1514,
    1243,     1768,
    1110,     2088,
    1349,     2406,
    1480,     2035,
    1558,     1762,
    1770,     1554,
    2027,     1319,
    1952,     1030,
    2281,      914,
    2337,     1209,
    2451,     1441,
    2234,     1441,
    2033,     1592,
    2279,     1644,
    2109,     1808,
    2319,     1859,
    2469,     1656,
    2533,     1857,
    2414,     2084,
    2114,     2087,
    1834,     1887,
    1769,     2242,
    1977,     2560,
    2278,     2344,
    2343,     2749,
    2661,     2416,
    2857,     2037,
    3239,     2405,
    2848,     2836,
    3426,     2863,
    3475,     3416,
    2831,     3409,
    2410,     3322,
    2033,     3303,
    1609,     3514,
    1625,     2836,
     880,     2775,
    1122,     3500,
     605,     3517,
      88,     3532,
     111,     2956,
      89,     2436,
      88,     1938,
      74,     1526};

//#include "ener_qua.tab"


/**************************************************************************
*
*	ROUTINE				:	Bits2prm_Tetra
*
*	DESCRIPTION			:	Convert serial received bits to the encoder
*							parameter vector
*
**************************************************************************
*
*	USAGE				:	Bits2prm_Tetra(st, buffer_in,buffer_out)
*							(Routine_Name(input1,output1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Serial bits (137 + bfi)
*							- Format : Word16 - .. * 16 bit-samples
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Encoded parameters
*								         (23 parameters + bfi)
*							- Format : Word16 - .. * 16 bit-samples
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	The received parameters are :
*
*						- BFI	bad frame indicator		1 bit
*
*						- LSP	1st codebook 			8 bits
*							2nd codebook			9 bits
*							3nd codebook			9 bits
*
*						- for the 4 subframes :
*							pitch delay			8 bits (first)
*											5 bits (others)
*							codebook index		14 bits
*							pulse global sign		1 bit
*							pulse shift			1 bit
*							pitch and innovation gains	6 bits
*
**************************************************************************/

#define TSC_PRM_NO 23

void Bits2prm_Tetra(tetra_codec* st, Word16 *bits, Word16 prm[])
{
  Word16 i;
  static const Word16 bitno[TSC_PRM_NO] = {8, 9, 9,            /* split VQ LSP  */
                                 8, 14, 1, 1, 6,     /* subframe 1    */
                                 5, 14, 1, 1, 6,     /* subframe 2    */
                                 5, 14, 1, 1, 6,     /* subframe 3    */
                                 5, 14, 1, 1, 6};    /* subframe 4    */
  *prm++ = *bits++;     /* read BFI */

  for (i = 0; i < TSC_PRM_NO; i++)
  {
    prm[i] = bin2int(st, bitno[i], bits);
    bits  += bitno[i];
  }
}

/**************************************************************************
*
*	ROUTINE				:	Cal_Rr2
*
*	DESCRIPTION			:	Compute the autocorrelation matrix of
*						impulse response h
*							Only the even elements are stored in the matrix
*
**************************************************************************
*
*	USAGE				:	Cal_Rr2(buffer_1,buffer_2)
*							(Routine_Name(input1,input2))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Impulse response
*						- Format : Word16
*
*		INPUT2			:	- Description : Autocorrelation matrix
*								          (passed as a vector)
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	Matrice rr[TSC_dim_rr][TSC_dim_rr]
*
*							<------- TSC_dim_rr ------->
*							<------ TSC_L/2 ---->
*
*						^	x x x x x x ... x 0 ...  0	^
*						|	x x x x x x ... x 0 ...  0	|
*						|	x x x x x x ... x 0 ...  0	|
*						|	x x x x x x ... x 0 ...  0	|
*						TSC_L/2	x x x x x x ... x 0 ...  0	|
*						|	x x x x x x ... x 0 ...  0	TSC_dim_rr
*						|	| | | | | |       | |      |	|
*						v	x x x x x x ... x 0 ...  0	|
*							0 0 0 0 0 0 ... 0 0 ... 0	|
*							| | | | | |     | |     |		|
*							0 0 0 0 0 0 ... 0 0 ... 0	v
*
*						Firstly, rr[TSC_L/2-1][TSC_L/2-1] is computed, with
*						the main diagonal step up. Then, the same
*						is done for the other diagonals.
*
**************************************************************************/

/* Length of subframe = 60 */
#define TSC_L      60

/* Dimension of matrix rr[TSC_dim_rr][TSC_dim_rr]  = 32 */
//#define TSC_dim_rr 32

/* Increment on diagonal */
#define TSC_diag   33

/* Index of TSC_last element to compute = rr[TSC_L/2-1][TSC_L/2-1] = rr[29*32+29] */
#define TSC_last   957


void Cal_Rr2(tetra_codec* st, Word16 h[], Word16 *rr)
{
 Word16 i, k, dec;
 Word32 s;
 Word16 *ph, *phd;
 Word16 *rr_fin_c;
 Word16 *rr_fin_r;
 Word16 *rr_ij;
 Word16 *rr_ji;
 Word16 hs[TSC_L];

 /* Scaling for maximum precision */

 s = 0;
 for(i=0; i<TSC_L; i++)
   s = L_mac0(st, s, h[i], h[i]);

 k = norm_l(s);
 k = shr(st, sub(st, k, (Word16)1), (Word16)1);

 for(i=0; i<TSC_L; i++)
  hs[i] = shl(st, h[i], k);

 /* Compute rr[][] */

 rr_fin_r = &rr[TSC_last];                  /* pointer on rr[TSC_L/2-1][TSC_L/2-1]*/
 rr_fin_c = rr_fin_r;

 for (dec = 0 ;dec < TSC_L; dec+=2)
 {
   ph  = hs;                            /* pointer on hs[i]      */
   phd = &hs[dec];                      /* pointer on hs[dec+i]  */
   rr_ij = rr_fin_c;                    /* end of column j       */
   rr_ji = rr_fin_r;                    /* end of row j          */
   s = 0;
   for (k = 0; k < TSC_L-dec; k+=2)
   {
     s = L_mac(st, s, *ph++, *phd++);
     s = L_mac(st, s, *ph++, *phd++);
     *rr_ij = extract_h(s);
     *rr_ji = extract_h(s);
     rr_ij -= TSC_diag;                     /* decrement pointers */
     rr_ji -= TSC_diag;
   }
   rr_fin_c--;                          /* decrement column */
   rr_fin_r -= TSC_dim_rr;                  /* decrement row    */
 }
 return;
}


/**************************************************************************
*
*	ROUTINE				:	Clsp_334
*
*	DESCRIPTION			:	Split vector quantization of LSP parameters
*							Use a split-3 VQ in the cosine lsp domain (-1,1),
*							without weighting
*
**************************************************************************
*
*	USAGE				:	Clsp_334(buffer_in,buffer_out1,buffer_out2)
*							(Routine_Name(input1,output1,output2))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : LSPs in the cosine domain (-1,1)
*						- Format : Word16 - Q15
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Quantized LSPs in the cosine domain
*						- Format : Word16 - Q15
*
*		OUTPUT2			:	- Description : Indices of the 3 selected
*								          codebook entries
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
*
**************************************************************************/

void Clsp_334(tetra_codec* st, Word16 *lsp, Word16 *lsp_q, Word16 *indice)
{

 Word16 i, j, ind = 0, temp;
 const Word16 *p_dico;
 Word32 min, dist;

 //static Word16 lsp_old[10]={
 //              30000, 26000, 21000, 15000, 8000, 0,
//	 	  -8000,-15000,-21000,-26000};


/* Search dico1  lsp[0]-lsp[2] */

 p_dico = dico1_clsp;
 min    = MAX_32;

 for (i = 0; i< TSC_taille_dic1; i++)
 {
    temp  = sub(st, lsp[0], *p_dico++);
    dist = L_mult0(temp,temp);
    for(j=1; j<3; j++)
    {
      temp  = sub(st, lsp[j], *p_dico++);
      dist  = L_mac0(st, dist,temp,temp);
    }
    if (L_sub(st, dist, min) < 0) { min = dist; ind = i; }
 }

 indice[0] = ind;
 p_dico    = &dico1_clsp[ind * 3];
 lsp_q[0]  = *p_dico++ ;
 lsp_q[1]  = *p_dico++ ;
 lsp_q[2]  = *p_dico++ ;

/* Search dico2 lsp[3]-lsp[5] */

 p_dico = dico2_clsp;
 min    = MAX_32;

 for (i = 0; i< TSC_taille_dic2; i++)
 {
    temp  = sub(st, lsp[3], *p_dico++);
    dist  = L_mult0(temp,temp);
    for(j=4; j<6; j++)
    {
      temp  = sub(st, lsp[j], *p_dico++);
      dist  = L_mac0(st, dist,temp,temp);
    }
    if (L_sub(st, dist, min) < 0) { min = dist; ind = i; }
 }

 indice[1] = ind;
 p_dico    = &dico2_clsp[ind * 3];
 lsp_q[3]  = *p_dico++ ;
 lsp_q[4]  = *p_dico++ ;
 lsp_q[5]  = *p_dico++ ;

/* Search dico3 lsp[6]-lsp[9] */

 p_dico = dico3_clsp;
 min    = MAX_32;

 for (i = 0; i< TSC_taille_dic3; i++)
 {
    temp  = sub(st, lsp[6], *p_dico++);
    dist  = L_mult0(temp,temp);
    for(j=7; j<10; j++)
    {
      temp  = sub(st, lsp[j], *p_dico++);
      dist  = L_mac0(st, dist,temp,temp);
    }
    if (L_sub(st, dist, min) < 0) { min = dist; ind = i; }
 }

 indice[2] = ind;
 p_dico    = &dico3_clsp[ind * 4];
 lsp_q[6]  = *p_dico++ ;
 lsp_q[7]  = *p_dico++ ;
 lsp_q[8]  = *p_dico++ ;
 lsp_q[9]  = *p_dico++ ;

 /* Minimum distance between lsp_q[2] and lsp_q[3] */

 temp = 917;                    /* 917 = 0.028 in Q15 = 50hz around 1000hz */
 temp = sub(st, temp, lsp_q[2]);
 temp = add(st, temp, lsp_q[3]);    /* temp = 0.028 - (lsp_q[2]-lsp_q[3])      */
 if (temp > 0)
 {
   temp = shr(st, temp, (Word16)1);
   lsp_q[2] = add(st, lsp_q[2], temp);
   lsp_q[3] = sub(st, lsp_q[3], temp);
 }
 /* Minimum distance between lsp_q[5] and lsp_q[6] */

 temp = 1245;                   /* 1245= 0.038 in Q15 = 50hz around 1600hz */
 temp = sub(st, temp, lsp_q[5]);
 temp = add(st, temp, lsp_q[6]);    /* temp = 0.038 - (lsp_q[5]-lsp_q[6])      */
 if (temp > 0)
 {
   temp = shr(st, temp, (Word16)1);
   lsp_q[5] = add(st, lsp_q[5], temp);
   lsp_q[6] = sub(st, lsp_q[6], temp);
 }

 /* Verify if lsp_q[] are still in order */

 temp = 0;
 for(i=0; i<9; i++)
 {
    if(sub(st, lsp_q[i],lsp_q[i+1]) <= 0 )
    {
      temp = 1;
    }
 }

  /* If lsp_q[] are not in order keep old lsp_q[] */

 if(temp != 0)
 {
   for(i=0; i<10; i++)
   lsp_q[i] = st->lsp_old[i];
 }

 else
 {
   for(i=0; i<10; i++)
   st->lsp_old[i] = lsp_q[i];
 }
 return;
}

/**************************************************************************
*
*	ROUTINE				:	Dec_Ener
*
*	DESCRIPTION			:	Decode gains of pitch and innovative codebook
*
**************************************************************************
*
*	USAGE				:	Dec_Ener(index,bfi,buffer_in1,buffer_in2,buffer_in3,
*							L_subfr,gain_p,gain_c)
*							(Routine_Name(input1,input2,input3,input4,input5,
*							input6,output1,output2))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Index of energy quantizer
*						- Format : Word16
*
*	INPUT2			:	- Description : Bad frame indicator
*						- Format : Word16
*
*	INPUT3			:	- Description : LPC filter
*						- Format : Word16
*
*	INPUT4			:	- Description : Adaptive codebook
*						- Format : Word16
*
*		INPUT5			:	- Description : Innovation codebook
*							- Format : Word16
*
*		INPUT6			:	- Description : Subframe length
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Quantized pitch gain
*							- Format : Word16
*
*		OUTPUT2			:	- Description : Quantized code gain
*							- Format : Word16
*
*	RETURNED VALUE		:	Index of energy quantizer
*
**************************************************************************/

//  extern Word16 last_ener_pit;
//  extern Word16 last_ener_cod;

Word16 Dec_Ener(tetra_codec* st, Word16 index, Word16 bfi, Word16 A[], Word16 prd_lt[],
          Word16 code[], Word16 L_subfr, Word16 *gain_pit, Word16 *gain_cod)
{

  Word16  i, j;
  Word16  exp, frac;
  Word16  exp_lpc, ener_lpc;
  Word16  exp_plt, ener_plt, pred_pit;
  Word16  ener_c, pred_cod;
  Word32  L_tmp;


  /*------------------------------------------------------*
   * Energy of impulse response of 1/A(z) for length = 60 *
   *------------------------------------------------------*/

  L_tmp = Lpc_Gain(st, A);

  exp_lpc  = norm_l(L_tmp);
  ener_lpc = extract_h(L_shl(st, L_tmp, exp_lpc) );

  /*-------------------------------*
   * Energy of adaptive codebook   *
   *-------------------------------*/

  L_tmp = 1;				/* Avoid case of all-zeros */

  for(i=0; i<L_subfr; i++)
    L_tmp = L_mac0(st, L_tmp, prd_lt[i], prd_lt[i]);

  exp_plt  = norm_l(L_tmp);
  ener_plt = extract_h( L_shl(st, L_tmp, exp_plt) );

  /* ener_plt = Log2(ener_plt * ener_lpc) */

  L_tmp = L_mult0(ener_plt, ener_lpc);
  exp_plt = add(st, exp_plt, exp_lpc);

  Log2(st, L_tmp, &exp, &frac);

  L_tmp = Load_sh16(st, exp);
  L_tmp = add_sh(st, L_tmp, frac, (Word16)1);	/* Log2(ener_plt*ener_lpc) in Q16*/
  L_tmp = sub_sh16(st, L_tmp, exp_plt);		/* subtract exponant of ener_plt */

  /* Input on 15 bits */
  L_tmp = add_sh(st, L_tmp, (Word16)1710, (Word16)8); /* +6.68 in Q16
									for other scaling    */

  L_tmp = L_shr(st, L_tmp, (Word16)8);		/* Convert from Q16 to Q8        */
  ener_plt = extract_l(L_tmp);

  /*-------------------------------*
   * Energy coming from code       *
   *-------------------------------*/

  L_tmp = 0;
  for(i=0; i<L_subfr; i++)
    L_tmp = L_mac0(st, L_tmp, code[i], code[i]);

  ener_c = extract_h(L_tmp);

  /* ener_c = Log2(ener_c * ener_lpc) */

  L_tmp = L_mult0(ener_c, ener_lpc);

  Log2(st, L_tmp, &exp, &frac);

  L_tmp = Load_sh16(st, exp);
  L_tmp = add_sh(st, L_tmp, frac, (Word16)1);	/* Log2(ener_plt*ener_lpc) in Q16*/
  L_tmp = sub_sh16(st, L_tmp, exp_lpc);		/* subtract exponant of ener_lpc */

  /* Input on 15 bits */
  L_tmp = sub_sh(st, L_tmp, (Word16)4434, (Word16)8); /*-17.32 in Q16
									for other scaling    */

  L_tmp = L_shr(st, L_tmp, (Word16)8);		/* Convert from Q16 to Q8        */
  ener_c = extract_l(L_tmp);


  /*-----------------------------------------------*
   * Test for bfi.                                 *
   *                                               *
   *  if (bfi != 0)                                *
   *    ->last_ener_pit -= 0.5                     *
   *    ->last_ener_cod -= 0.5                     *
   *  else                                         *
   *  {                                            *
   *    decode new last_ener_pit et last_ener_cod  *
   *  }                                            *
   *-----------------------------------------------*/

   if(bfi != 0)
   {
     st->last_ener_pit = sub(st, st->last_ener_pit, (Word16)128);	/* -0.5 in Q8 */
     if(st->last_ener_pit < 0) st->last_ener_pit = 0;

     st->last_ener_cod = sub(st, st->last_ener_cod, (Word16)128);	/* -0.5 in Q8 */
     if(st->last_ener_cod < 0) st->last_ener_cod = 0;
   }
   else
   {

     /*-----------------------------------------------------------------*
      * Prediction on pitch energy.                                     *
      *  pred_pit = 0.50 * last_ener_pit + 0.25 * last_ener_cod - 3.0   *
      *  if(pred_pit < 0.0) pred_pit = 0.0;                             *
      *-----------------------------------------------------------------*/

     L_tmp = Load_sh(st, st->last_ener_pit, (Word16)8);	/* .5 last_ener_pit in Q9  */
     L_tmp = add_sh(st, L_tmp,st->last_ener_cod, (Word16)7); /*+.25 last_ener_code
											in Q9    */
     L_tmp = sub_sh(st, L_tmp, (Word16)768, (Word16)9); /* -3.0 in Q9          */
     if(L_tmp < 0) L_tmp = 0;
     pred_pit = store_hi(st, L_tmp, (Word16)7);	/* result in Q8 		   */

     /*-----------------------------------------------------------------*
      * Prediction on code energy.                                      *
      *  pred_cod = 0.50 * last_ener_cod + 0.25 * last_ener_pit - 3.0   *
      *  if(pred_cod < 0.0) pred_cod = 0.0;                             *
      *-----------------------------------------------------------------*/

     L_tmp = Load_sh(st, st->last_ener_cod, (Word16)8);	/* .5 last_ener_cod in Q9  */
     L_tmp = add_sh(st, L_tmp, st->last_ener_pit, (Word16)7); /*+.25 last_ener_pit
											in Q9    */
     L_tmp = sub_sh(st, L_tmp, (Word16)768, (Word16)9); /* -3.0 in Q9          */
     if(L_tmp < 0) L_tmp = 0;
     pred_cod = store_hi(st, L_tmp, (Word16)7);	/* result in Q8		   */

     /*-----------------------------------------------------------------*
      * Read quantized values.                                          *
      *-----------------------------------------------------------------*/

     j = shl(st, index, (Word16)1);
     st->last_ener_pit = add(st, t_qua_ener[j],   pred_pit);
     st->last_ener_cod = add(st, t_qua_ener[j+1], pred_cod);

     /* Limit energies ->for transmission errors */

     if(sub(st, st->last_ener_pit, (Word16)6912)>0)st->last_ener_pit = 6912;
									/* 6912 = 27 in Q8 */
     if(sub(st, st->last_ener_cod, (Word16)6400)>0)st->last_ener_cod = 6400;
									/* 6400 = 25 in Q8 */
  }


  /*---------------------------------------------------*
   *                                                   *
   *  Compute the quantized pitch gain.                *
   *                                                   *
   *     temp = 0.5 * ( last_ener_pit - ener_plt );    *
   *     temp = pow(2.0, temp);                        *
   *     if( temp > 1.2) temp = 1.2;                   *
   *     *gain_pit = temp;                             *
   *                                                   *
   *---------------------------------------------------*/

  L_tmp = Load_sh(st, st->last_ener_pit, (Word16)6);	 /* last_ener_pit/2 in Q15 */
  L_tmp = sub_sh(st, L_tmp, ener_plt, (Word16)6);	 /* - ener_plt/2    in Q15 */
  L_tmp = add_sh(st, L_tmp, (Word16)12, (Word16)15); /* to have gain in Q12    */
  L_extract(st, L_tmp, &exp, &frac);
  L_tmp = pow2(st, exp, frac);
  if( L_sub(st, L_tmp, (Word16)4915) > 0) L_tmp = 4915; /* 4915 = 1.2 in Q12   */
  *gain_pit = extract_l(L_tmp);

  /*---------------------------------------------------*
   *                                                   *
   *  Compute the innovative code gain.                *
   *                                                   *
   *     temp = 0.5 * ( last_ener_cod - ener_c );      *
   *     *gain_cod = pow(2.0, temp);                   *
   *                                                   *
   *---------------------------------------------------*/

  L_tmp = Load_sh(st, st->last_ener_cod, (Word16)6);	/* last_ener_cod/2 in Q15 */
  L_tmp = sub_sh(st, L_tmp, ener_c, (Word16)6);	/* - ener_c/2      in Q15 */
  L_extract(st, L_tmp, &exp, &frac);
  L_tmp = pow2(st, exp, frac);
  *gain_cod = extract_l(L_tmp);

  return index;
}


/**************************************************************************
*
*	ROUTINE				:	D4i60_16
*
*	DESCRIPTION			:	Innovation codebook search
*
**************************************************************************
*
*	USAGE				:	D4i60_16(buffer_in1,buffer_in2,buffer_in3,buffer_in4,
*							buffer_out1, buffer_out2,sign,shift)
*							(Routine_Name(input1,input2,input3,input4,
*							output1,output2,output3,output4))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Result of backward-filtering
*						- Format : Word16
*
*	INPUT2			:	- Description : Impulse response of noise filter
*							          (vector f[])
*						- Format : Word16 - (must be preceeded by 64 zeros)
*
*	INPUT3			:	- Description : Total impulse response (vector h[])
*						- Format : Word16 - (must be preceeded by 64 zeros)
*
*	INPUT4			:	- Description : Autocorrelations of h[]
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Innovation code found
*							  by convolving with f[]
*							- Format : Word16
*
*		OUTPUT2			:	- Description : Innovation code found
*							  by convolving with h[]
*							- Format : Word16
*
*		OUTPUT3			:	- Description : Innovation code gain sign
*							- Format : Word16
*
*		OUTPUT4			:	- Description : Shift of innovation code
*							- Format : Word16
*
*	RETURNED VALUE		:	Index of innovation code
*
**************************************************************************


**************************************************************************
*
*	COMMENTS			:	The code length is 60, containing 4 nonzero pulses
*
*						i0, i1, i2, i3, with fixed signs.
*						i0 can take 30 positions.
*						i1, i2 and i3 can have 8 possible positions :
*
*						i0 (+)	: 0,  2, 4.  .... 58
*						i1 (-)	: 2, 10, 18, 26, 34, 42, 50, 58
*						i2 (+)	: 4, 12, 20, 28, 36, 44, 52, (60)
*						i3 (-)	: 6, 14, 22, 30, 38, 46, 54, (62)
*
*						Positions 60 and 62 correspond, respectively,
*						to pulses i2 and i3 not used.
*						The code can be shifted 1 position to the right.
*
*	NOTES ON COMPLEXITY		:
*
*	1 :	In this routine all accesses to matrix rr[][] use indices. This was done for clarity and
*		readability of the program. In the actual real-time implementation all these
*		accesses use pointers.  Parallel autoincrement of pointers on dn[] and rr[][]
*		should be used, especially for the two most inner loops in the search which are very
*		critical.
*
*		As the number of addressing registers and the possible modes of autoincrement
*		can vary with each DSP, the addressing of matrix rr[][] should be optimized
*		for the DSP used.
*
*	2 :	Because of the thresholds it is difficult to estimate the average complexity.
*		Usually, we mesure it with the real-time implementation. The worst
*		case corresponds to "time <= 0".  The decrements of 3 and 4 on variable
*		"time" come from the relative number of operations when "seuil2" is exceeded and
*		the number of operations when "seuil1" is exceeded, excluding the operations when
*		"seuil2" is exceeded. We supose that a new maximum is found 1 out of 8 new tries.
*		The numbers 3 and 4 correspond to the real-time implementation and not to this
*		C code.
*
**************************************************************************/

#define TSC_lcode 60

/* Gain for impulse 0 (i0) = sqrt(2) = 1.4142                             */
/* We use different values to be compatible with the implementation       */
/* which uses 3 different Qx formats for this number to simplify scaling. */

#define TSC_Q11_gain_i0  2896
#define TSC_Q13_gain_i0 11585
#define TSC_Q14_gain_i0 23170

/*------------------------------------------------------------------*
 * The next 3 parameters control the time taken by the innovative   *
 * codebook search.                                                 *
 * The first two parameters control if a section of the innovative  *
 * codebook should be search or not. The TSC_last parameter controls    *
 * the absolute maximum time used by the innovative codebook search.*
 *------------------------------------------------------------------*/

 /* thresholds = 0.586 in Q15 */

#define TSC_threshold1 19200
#define TSC_threshold2 19200
#define TSC_max_time   350

Word16 D4i60_16(tetra_codec* st, Word16 dn[], Word16 f[], Word16 h[], Word16 rr[][32],
            Word16 cod[], Word16 y[], Word16 *sign, Word16 *shift_code)
{
 Word16 ip0, ip1, ip2, ip3;
 Word16 ii0, ii1, ii2, ii3;
 Word16 *p0, *p1, *p2, *p3;
 Word16 i,j;
 Word16 shif, shift;
 Word16 time, index;
 Word16 ps0,  ps1,  ps2,  ps3;
 Word16 ps0a, ps1a, ps2a;
 Word16 ps, ps3c, psc;
 Word16 alpha, alp1, alp2, alp3;
 Word16 seuil1,seuil2;
 Word16 max0, max1, max2, min0, min1, min2;
 Word32 L_tmp;
 Word32 alp0_32, alp1_32, alp2_32, ps2_32;

 /* set dn[60] - dn[63] to zero */
 dn[60]=dn[61]=dn[62]=dn[63] = 0;

 /* Find min. and max. of dn[] for positions of i0, i1, i2 */

 min0=max0=0;

 for(i = 0; i<60; i += 2)
 {
   j = add(st, dn[i], dn[i+1]);
   if      (sub(st, j,max0) > 0) max0 = j;
   else if (sub(st, j,min0) < 0) min0 = j;
 }
 max0 = shr(st, max0, (Word16)1);
 min0 = shr(st, min0, (Word16)1);


 /* Multiply max0 and min0 by gain_i0 in Q11 */

 max0 = store_hi(st, L_mult0(max0, (Word16)TSC_Q11_gain_i0), (Word16)5);
 min0 = store_hi(st, L_mult0(min0, (Word16)TSC_Q11_gain_i0), (Word16)5);

 min1=max1=0;
 for(i = 2; i<60; i += 8)
 {
   j = add(st, dn[i], dn[i+1]);
   if      (sub(st, j,max1) > 0) max1 = j;
   else if (sub(st, j,min1) < 0) min1 = j;
 }
 max1 = shr(st, max1, (Word16)1);
 min1 = shr(st, min1, (Word16)1);

 min2=max2=0;
 for(i = 4; i<60; i += 8)
 {
   j = add(st, dn[i], dn[i+1]);
   if      (sub(st, j,max2) > 0) max2 = j;
   else if (sub(st, j,min2) < 0) min2 = j;
 }
 max2 = shr(st, max2, (Word16)1);
 min2 = shr(st, min2, (Word16)1);


/*----------------------------------------------------------*
 * Find absolute maximum for combination dn[i0] and dn[i1]  *
 *      and for combination dn[i0], dn[i1] and dn[i2]       *
 *   max1 = max ( max0-min1 , max1-min0 )                   *
 *   max2 = max ( max0-min1+max2 , max1-min0-min2 )         *
 *----------------------------------------------------------*/

 max0 = sub(st, max0, min1);
 max2 = add(st, max0, max2);
 max1 = sub(st, max1, min0);
 j    = sub(st, max1, min2);

 if(sub(st, max0,max1) > 0) max1 = max0;
 if(sub(st, j,max2)    > 0) max2 = j;

 /* Set thresholds */
 seuil1 = mult(st, max1, (Word16)TSC_threshold1);
 seuil2 = mult(st, max2, (Word16)TSC_threshold2);

 /* Default values */

 ip0    = 0;
 ip1    = 2;
 ip2    = 4;
 ip3    = 6;
 shift  = 0;
 ps     = 0;
 psc    = 0;
 alpha  = 255;
 time   = TSC_max_time;

 /* Four loops to search innovation code. */

 for (ii0=0; ii0<30; ii0++)
 {
   ps0  = store_hi(st, L_mult0((Word16)TSC_Q11_gain_i0, dn[2*ii0]), (Word16)5);
   ps0a = store_hi(st, L_mult0((Word16)TSC_Q11_gain_i0, dn[2*ii0+1]), (Word16)5);
   alp0_32 = Load_sh(st, rr[ii0][ii0], (Word16)14);

   for (ii1=1; ii1<30; ii1+=4)
   {
     ps1  = sub(st, ps0, dn[2*ii1]);
     ps1a = sub(st, ps0a, dn[2*ii1+1]);

     /* alp1 = alp0*2 + rr[ii1][ii1] - 2.0*gain_i0*rr[ii0][ii1]  */
     /* The result is divided by 8 to avoid overflow.            */
     L_tmp = add_sh(st, alp0_32, rr[ii1][ii1], (Word16)13);

     L_tmp = L_msu0(st, L_tmp, (Word16)TSC_Q14_gain_i0, rr[ii0][ii1]);
     alp1  = extract_h(L_tmp);
     alp1_32 = Load_sh(st, alp1, (Word16)15);

    /* If( abs( (ps1+ps1a)/2 ) > seuil1 */

     L_tmp = Load_sh(st, ps1, (Word16)15);
     L_tmp = L_abs( add_sh(st, L_tmp, ps1a, (Word16)15) );
     L_tmp = sub_sh16(st, L_tmp, seuil1);
     if( L_tmp > 0)
     {
       for (ii2=2; ii2<31; ii2+=4)
       {
         ps2  = add(st, ps1, dn[2*ii2]);
         ps2a = add(st, ps1a, dn[2*ii2+1]);

         /* alp2=alp1+rr[ii2][ii2]+2.0*(gain_i0*rr[ii0][ii2]-rr[ii1][ii2]) */
         /* The result is divided by 16 to avoid overflow.                 */

         L_tmp = add_sh(st, alp1_32, rr[ii2][ii2], (Word16)12);

         L_tmp = L_mac0(st, L_tmp, (Word16)TSC_Q13_gain_i0, rr[ii0][ii2]);
         L_tmp = sub_sh(st, L_tmp, rr[ii1][ii2] , (Word16)13);
         alp2  = extract_h(L_tmp);
         alp2_32 = Load_sh16(st, alp2);

         /* If( abs( (ps2+ps2a)/2 ) > seuil2 */

         L_tmp = Load_sh(st, ps2, (Word16)15);
         L_tmp = L_abs( add_sh(st, L_tmp, ps2a, (Word16)15) );
         L_tmp = sub_sh16(st, L_tmp, seuil2);
         if( L_tmp > 0)
         {
           shif = 0;
           if( sub(st, abs_s(ps2a),abs_s(ps2)) > 0)
           {
             ps2  = ps2a;
             shif = 1;
           }

           ps2_32 = Load_sh(st, ps2, (Word16)15);

           for (ii3=3; ii3<32; ii3+=4)
           {
             /* ps3 = (ps2-dn[2*ii3+shift]) / 2                       */
             /* Need to work on 32 bits to avoid possible overflow */

             L_tmp = sub_sh(st, ps2_32, dn[2*ii3+shif], (Word16)15);

             ps3  = extract_h(L_tmp);

             /* alp3 = alp2 + rr[ii3][ii3] +                               */
             /*   2.0*(rr[ii1][ii3] - gain_i0*rr[ii0][ii3] - rr[ii2][ii3]) */
             /*                                                            */
             /* This part is the most critical part in all the coder       */

             L_tmp = add_sh(st, alp2_32, rr[ii3][ii3], (Word16)12);
             L_tmp = add_sh(st, L_tmp, rr[ii1][ii3], (Word16)13);
             L_tmp = L_msu0(st, L_tmp, (Word16)TSC_Q13_gain_i0, rr[ii0][ii3] );
             L_tmp = sub_sh(st, L_tmp, rr[ii2][ii3], (Word16)13);
             alp3  = extract_h(L_tmp);

             /* if( (ps3**2 * alpha) - (psc *alp3)  > 0 ) */
             /*    ->new maximum                          */

             ps3c = mult(st, ps3, ps3);
             L_tmp  = L_mult(st, ps3c, alpha);
             if( L_msu(st, L_tmp,psc,alp3) > 0)
             {
               ps    = ps3;
               psc   = ps3c;
               alpha = alp3;
               ip0 = 2*ii0;
               ip1 = 2*ii1;
               ip2 = 2*ii2;
               ip3 = 2*ii3;
               shift = shif;
             } /*  end of if(.. > 0) */
           } /* end of for ii3 = */

           time = sub(st, time, (Word16)3);	/* See note on complexity */
           if(time <= 0 ) goto fin;		/* Maximum time finish    */

         } /* end if >seuil2 */
       } /* end of for ii2 = */

       time = sub(st, time, (Word16)4);		/* See note on complexity */
       if(time <= 0 ) goto fin;		/* Maximum time finish    */

     } /* end if >seuil1 */
   } /* end of for ii1 = */
 } /* end of for ii0 = */

 fin:

 /* Convolve code with f[] */

 f -= shift;            /* Operations on pointers !!! */
 p0 = f - ip0;
 p1 = f - ip1;
 p2 = f - ip2;
 p3 = f - ip3;

 /* cod[i] = p0[i] * gain_i0  - p1[i] + p2[i] - p3[i]; */

 if(ps >= 0)
 {
   *sign = 0;
   for (i = 0; i < TSC_lcode; i++)
   {
     L_tmp  = L_mult0(p0[i], (Word16)TSC_Q11_gain_i0);
     L_tmp  = sub_sh(st, L_tmp, p1[i], (Word16)11);
     L_tmp  = add_sh(st, L_tmp, p2[i], (Word16)11);
     L_tmp  = sub_sh(st, L_tmp, p3[i], (Word16)11);
     cod[i] = store_hi(st, L_tmp, (Word16)5);
   }
 }
 else
 {
   *sign = 1;
   for (i = 0; i < TSC_lcode; i++)
   {
     L_tmp  = L_mult0(p0[i], (Word16)TSC_Q11_gain_i0);
     L_tmp  = sub_sh(st, L_tmp, p1[i], (Word16)11);
     L_tmp  = add_sh(st, L_tmp, p2[i], (Word16)11);
     L_tmp  = sub_sh(st, L_tmp, p3[i], (Word16)11);
     L_tmp  = L_negate(L_tmp);
     cod[i] = store_hi(st, L_tmp, (Word16)5);
   }
 }

 /* Convolve code with h[] */

 h -= shift;            /* Operations on pointers !!! */
 p0 = h - ip0;
 p1 = h - ip1;
 p2 = h - ip2;
 p3 = h - ip3;

 /* y[i] = p0[i] * gain_i0  - p1[i] + p2[i] - p3[i]; */

 if(ps >= 0)
 {
   for (i = 0; i < TSC_lcode; i++)
   {
     L_tmp = L_mult0(p0[i], (Word16)TSC_Q11_gain_i0);
     L_tmp = sub_sh(st, L_tmp, p1[i], (Word16)11);
     L_tmp = add_sh(st, L_tmp, p2[i], (Word16)11);
     L_tmp = sub_sh(st, L_tmp, p3[i], (Word16)11);
     y[i]  = store_hi(st, L_tmp, (Word16)5);
   }
 }
 else
 {
   for (i = 0; i < TSC_lcode; i++)
   {
     L_tmp = L_mult0(p0[i], (Word16)TSC_Q11_gain_i0);
     L_tmp = sub_sh(st, L_tmp, p1[i], (Word16)11);
     L_tmp = add_sh(st, L_tmp, p2[i], (Word16)11);
     L_tmp = sub_sh(st, L_tmp, p3[i], (Word16)11);
     L_tmp = L_negate(L_tmp);
     y[i]  = store_hi(st, L_tmp, (Word16)5);
   }
 }



 /* Return parameters */

 *shift_code = shift;

 /* index = (ip0>>1) + ( (ip1>>3)<<5 ) + ( (ip2>>3)<<8 ) + ( (ip3>>3)<<11 )*/

 index = shr(st, ip0, (Word16)1);
 index = add(st, index, shl(st,  shr(st,  ip1,(Word16)3 ),(Word16)5 ));
 index = add(st, index, shl(st,  shr(st,  ip2,(Word16)3 ),(Word16)8 ));
 index = add(st, index, shl(st,  shr(st,  ip3,(Word16)3 ),(Word16)11 ));

 return index;
}


/**************************************************************************
*
*	ROUTINE				:	D_D4i60
*
*	DESCRIPTION			:	Decode innovative codebook d4i60_16
*
**************************************************************************
*
*	USAGE				:	D_D4i60(index,sign, shift,buffer_in,buffer_out)
*							(Routine_Name(input1,input2,input3,input4,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Index of codebook
*						- Format : Word16
*
*	INPUT2			:	- Description : Sign of codebook
*						- Format : Word16
*
*	INPUT3			:	- Description : Shift of codebook
*						- Format : Word16
*
*	INPUT4			:	- Description : Noise filter
*						- Format : Word16 - (must be preceeded by 64 zeros)
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Innovative vector
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
*
**************************************************************************/

#define TSC_lcode 60

/* Gain for impulse 0 = sqrt(2) = 1.4142 = 2896 in Q11 */

#define TSC_Q11_gain_i0  2896

void D_D4i60(tetra_codec* st, Word16 index, Word16 sign, Word16 shift, Word16 F[], Word16 cod[])
{
  Word16 i, pos0, pos1, pos2, pos3;
  Word16 *p0, *p1, *p2, *p3;
  Word32 L_tmp;

  /* Position of the 4 impulses */

  pos0 = shl(st,  (Word16)(index & 31), (Word16)1);

  pos1 = shr(st,  (Word16)(index & 224), (Word16)2);
  pos1 = add(st, pos1, (Word16)2);

  pos2 = shr(st,  (Word16)(index & 1792), (Word16)5);
  pos2 = add(st, pos2, (Word16)4);

  pos3 = shr(st,  (Word16)(index & 14336), (Word16)8);
  pos3 = add(st, pos3, (Word16)6);


  /* Convolve code with F[]                             */
  /* cod[i] = p0[i] * gain_i0  - p1[i] + p2[i] - p3[i]; */

  F -= shift;            /* Operations on pointers !!!  */
  p0 = F - pos0;
  p1 = F - pos1;
  p2 = F - pos2;
  p3 = F - pos3;

  if(sign == 0)
  {
    for (i = 0; i < TSC_lcode; i++)
    {
      L_tmp  = L_mult0(p0[i], (Word16)TSC_Q11_gain_i0);
      L_tmp  = sub_sh(st, L_tmp, p1[i], (Word16)11);
      L_tmp  = add_sh(st, L_tmp, p2[i], (Word16)11);
      L_tmp  = sub_sh(st, L_tmp, p3[i], (Word16)11);
      cod[i] = store_hi(st, L_tmp, (Word16)5);
    }
  }
  else
  {
    for (i = 0; i < TSC_lcode; i++)
    {
      L_tmp  = L_mult0(p0[i], (Word16)TSC_Q11_gain_i0);
      L_tmp  = sub_sh(st, L_tmp, p1[i], (Word16)11);
      L_tmp  = add_sh(st, L_tmp, p2[i], (Word16)11);
      L_tmp  = sub_sh(st, L_tmp, p3[i], (Word16)11);
      L_tmp  = L_negate(L_tmp);
      cod[i] = store_hi(st, L_tmp, (Word16)5);
    }
  }

  return;
}

/**************************************************************************
*
*	ROUTINE				:	D_Lsp334
*
*	DESCRIPTION			:	Decoding: Split vector quantization of
*							LSP parameters
*
**************************************************************************
*
*	USAGE				:	D_Lsp334(st, buffer_in1,buffer_out1,buffer_in2)
*							(Routine_Name(input1,output1,input2))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : indices of the three selected
*							          codebook entries
*						- Format : Word16
*
*		INPUT2			:	- Description : Previous LSP values
*								          (in the cosine domain)
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Quantized LSPs (in the cosine domain)
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
*
**************************************************************************/

void D_Lsp334(tetra_codec* st, Word16 indice[], Word16 lsp[], Word16 lsp_old[])
{
   Word16 i/*, *p*/, temp;
   const Word16* p;

   p      = &dico1_clsp[ indice[0] * 3];
   lsp[0] = *p++ ;
   lsp[1] = *p++ ;
   lsp[2] = *p++ ;

   p      = &dico2_clsp[ indice[1] * 3];
   lsp[3] = *p++ ;
   lsp[4] = *p++ ;
   lsp[5] = *p++ ;

   p      = &dico3_clsp[ indice[2] * 4];
   lsp[6] = *p++ ;
   lsp[7] = *p++ ;
   lsp[8] = *p++ ;
   lsp[9] = *p++ ;

   /* Minimum distance between lsp[2] and lsp[3] */

   temp = 917;			/* 917 = 0.028 in Q15 = 50hz around 1000hz */
   temp = sub(st, temp, lsp[2]);
   temp = add(st, temp, lsp[3]);	/* temp = 0.028 - (lsp[2]-lsp[3])      */
   if (temp > 0)
   {
     temp = shr(st, temp, (Word16)1);
     lsp[2] = add(st, lsp[2], temp);
     lsp[3] = sub(st, lsp[3], temp);
   }


   /* Minimum distance between lsp[5] and lsp[6] */

   temp = 1245;			/* 1245= 0.038 in Q15 = 50hz around 1600hz */
   temp = sub(st, temp, lsp[5]);
   temp = add(st, temp, lsp[6]);	/* temp = 0.038 - (lsp[5]-lsp[6])      */
   if (temp > 0)
   {
     temp = shr(st,  temp,(Word16)1 );
     lsp[5] = add(st, lsp[5], temp);
     lsp[6] = sub(st, lsp[6], temp);
   }

   /* Verify if lsp[] are still in order */

   temp = 0;
   for(i=0; i<9; i++)
   {
      if(sub(st, lsp[i],lsp[i+1]) <= 0 )
      {
        temp = 1;
      }
   }

  /* If lsp[] are not in order keep old lsp[] */

   if(temp != 0)
   {
     for(i=0; i<10; i++)
       lsp[i] = lsp_old[i];
   }
   return;
}


/**************************************************************************
*
*	ROUTINE				:	Ener_Qua
*
*	DESCRIPTION			:	Compute quantized gains of pitch and innovative
*						 	codebook resulting from the quantization of energies
*							of adaptive and innovation codebook
*
**************************************************************************
*
*	USAGE				:	Ener_Qua(buffer_in1,buffer_in2,buffer_in3,L_subfr,
*							gain_p,gain_c)
*							(Routine_Name(input1,input2,input3,input4,
*							arg5,arg6))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : LPC filter
*						- Format : Word16
*
*	INPUT2			:	- Description : Adaptive codebook
*						- Format : Word16
*
*	INPUT3			:	- Description : Innovation codebook
*						- Format : Word16
*
*	INPUT4			:	- Description : Subframe length
*						- Format : Word16
*
*		ARG5				:	- Description : Pitch gain
*							- Format : Word16
*
*		ARG6				:	- Description : Code gain
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		ARG5				:	- Description : Quantized pitch gain
*							- Format : Word16
*
*		ARG6				:	- Description : Quantized code gain
*							- Format : Word16
*
*	RETURNED VALUE		:	Index of quantization
*
**************************************************************************/

//  extern Word16 last_ener_pit;
//  extern Word16 last_ener_cod;

Word16 Ener_Qua(tetra_codec* st, Word16 A[], Word16 prd_lt[], Word16 code[], Word16 L_subfr,
                Word16 *gain_pit, Word16 *gain_cod)
{

  Word16  i, j, index = 0, tmp/*, *p*/;
  const Word16* p;
  Word16  exp, frac;
  Word16  exp_lpc, ener_lpc;
  Word16  exp_plt, ener_plt, ener_pit, err_pit, pred_pit;
  Word16  ener_c, ener_cod, err_cod, pred_cod;
  Word32  L_tmp, dist, dist_min;


  /*------------------------------------------------------*
   * Energy of impulse response of 1/A(z) for length = 60 *
   *------------------------------------------------------*/

  L_tmp = Lpc_Gain(st, A);

  exp_lpc  = norm_l(L_tmp);
  ener_lpc = extract_h( L_shl(st, L_tmp, exp_lpc) );

  /*-------------------------------*
   * Energy coming from pitch      *
   *-------------------------------*/

  L_tmp = 1;            /* Avoid case of all-zeros */

  for(i=0; i<L_subfr; i++)
    L_tmp = L_mac0(st, L_tmp, prd_lt[i], prd_lt[i]);

  exp_plt  = norm_l(L_tmp);
  ener_plt = extract_h( L_shl(st, L_tmp, exp_plt) );

  /* ener_plt = Log2(ener_plt * ener_lpc) */

  L_tmp = L_mult0(ener_plt, ener_lpc);
  exp_plt = add(st, exp_plt, exp_lpc);

  Log2(st, L_tmp, &exp, &frac);

  L_tmp = Load_sh16(st, exp);
  L_tmp = add_sh(st, L_tmp, frac, (Word16)1);	/*Log2(ener_plt*ener_lpc) in Q16 */
  L_tmp = sub_sh16(st, L_tmp, exp_plt);		/* subtract exponant of ener_plt */

  /* Input on 15 bits */
  L_tmp = add_sh(st, L_tmp, (Word16)1710, (Word16)8);	/* +6.68 in Q16
									for other scaling    */

  /* The scaling factor +6.68 includes many scalings, due to the fact that
     impulse response of the LP filter is in Q10 (see routine Lpc_Gain), and
     that the original quantization gain table was designed using -26 dB data
     reduced to 13 bits while the final coder works with -22dB signal reduced
     to 15 bits */

  L_tmp = L_shr(st, L_tmp, (Word16)8);		/* Convert from Q16 to Q8        */
  ener_plt = extract_l(L_tmp);

  /* ener_pit = Log2(gain_pit**2) + ener_plt */

  L_tmp = 1;					/* Avoid Log2(0) 			   */
  L_tmp = L_mac0(st, L_tmp, *gain_pit, *gain_pit);

  Log2(st, L_tmp, &exp, &frac);

  L_tmp = Load_sh16(st, exp);
  L_tmp = add_sh(st, L_tmp, frac, (Word16)1);	/* Log2(gain_pit**2)    in Q16   */
  L_tmp = sub_sh16(st, L_tmp, (Word16)24);	/* -24.0 -> gain_pit**2 is in Q24*/

  /* Pitch gain (variable "gain_pit")is in Q12, so its square is in Q24,
     which means that the real value is multplied by 2**24. Thus, in applying
     the function Log2(), a value of 24 is being added. A scaling factor of
     -24 is then needed to compensate  */

  L_tmp = L_shr(st, L_tmp, (Word16)8);		/* Convert from Q16 to Q8        */
  ener_pit = extract_l(L_tmp);
  ener_pit = add(st, ener_pit, ener_plt);


  /*-------------------------------*
   * Energy coming from code       *
   *-------------------------------*/

  L_tmp = 0;
  for(i=0; i<L_subfr; i++)
    L_tmp = L_mac0(st, L_tmp, code[i], code[i]);

  ener_c = extract_h(L_tmp);

  /* ener_c = Log2(ener_c * ener_lpc) */

  L_tmp = L_mult0(ener_c, ener_lpc);

  Log2(st, L_tmp, &exp, &frac);

  L_tmp = Load_sh16(st, exp);
  L_tmp = add_sh(st, L_tmp, frac, (Word16)1);	/* Log2(ener_plt*ener_lpc) in Q16*/
  L_tmp = sub_sh16(st, L_tmp, exp_lpc);		/* subtract exponant of ener_lpc */

  /* Input on 15 bits */
  L_tmp = sub_sh(st, L_tmp, (Word16)4434, (Word16)8);	/*-17.32 in Q16
									for other scaling    */

  /* The vector "code[i]" is in Q12, so its square (variable "ener_c") is in
     Q24, which means that the real value is multplied by 2**24. Thus, in
     applying the function Log2(), a value of 24 is being added. A scaling
     factor of -24 is then needed to compensate, but as the other scaling
     factor of +6.68 (as explained above) has also to be applied, the final
     resulting factor is -17.32 */

  L_tmp = L_shr(st, L_tmp, (Word16)8);		/* Convert from Q16 to Q8        */
  ener_c = extract_l(L_tmp);

  /* ener_cod = Log2(gain_code**2) + ener_c */

  L_tmp = 1;					/* Avoid Log2(0) 			   */
  L_tmp = L_mac0(st, L_tmp, *gain_cod, *gain_cod);

  Log2(st, L_tmp, &exp, &frac);

  L_tmp = Load_sh16(st, exp);
  L_tmp = add_sh(st, L_tmp, frac, (Word16)1);	/* Log2(gain_code**2) in Q16	   */
  L_tmp = L_shr(st, L_tmp, (Word16)8);		/* Convert from Q16 to Q8        */
  ener_cod = extract_l(L_tmp);
  ener_cod = add(st, ener_cod, ener_c);

  /*-----------------------------------------------------------------*
   * Prediction on pitch energy and prediction error                 *
   *  pred_pit = 0.50 * last_ener_pit + 0.25 * last_ener_cod - 3.0   *
   *  if(pred_pit < 0.0) pred_pit = 0.0;                             *
   *  err_pit = ener_pit - pred_pit;                                 *
   *-----------------------------------------------------------------*/

  L_tmp = Load_sh(st, st->last_ener_pit, (Word16)8); /* .5 last_ener_pit in Q9     */
  L_tmp = add_sh(st, L_tmp, st->last_ener_cod, (Word16)7); /*+.25 last_ener_code
											in Q9    */
  L_tmp = sub_sh(st, L_tmp, (Word16)768, (Word16)9); /* -3.0 in Q9             */
  if(L_tmp < 0) L_tmp = 0;
  pred_pit = store_hi(st, L_tmp, (Word16)7);	/* result in Q8 			   */
  err_pit  = sub(st, ener_pit, pred_pit);


  /*-----------------------------------------------------------------*
   * Prediction on code energy and prediction error                  *
   *  pred_cod = 0.50 * last_ener_cod + 0.25 * last_ener_pit - 3.0   *
   *  if(pred_cod < 0.0) pred_cod = 0.0;                             *
   *  err_cod = ener_cod - pred_cod;                                 *
   *-----------------------------------------------------------------*/

  L_tmp = Load_sh(st, st->last_ener_cod, (Word16)8); /* .5 last_ener_cod  in Q9    */
  L_tmp = add_sh(st, L_tmp, st->last_ener_pit, (Word16)7); /*+.25 last_ener_pit
											in Q9    */
  L_tmp = sub_sh(st, L_tmp, (Word16)768, (Word16)9); /* -3.0 in Q9             */
  if(L_tmp < 0) L_tmp = 0;
  pred_cod = store_hi(st, L_tmp, (Word16)7);	/* result in Q8 			   */
  err_cod  = sub(st, ener_cod, pred_cod);


  /*-----------------------------------------------------------------*
   * Codebook search.                                                *
   * Find vector which minimizes:                                    *
   *                                                                 *
   *   Min k :  (err_pit - dico[0,k])**2 +  (err_cod - dico[1,k])**2 *
   *                                                                 *
   *-----------------------------------------------------------------*/

  dist_min = MAX_32;
  p = t_qua_ener;

  for (i = 0; i< TSC_nb_qua_ener; i++)
  {
     tmp  = sub(st, *p++, err_pit);
     dist = L_mult0(tmp, tmp);
     tmp  = sub(st, *p++, err_cod);
     dist = L_mac0(st, dist, tmp, tmp);



     if (L_sub(st, dist, dist_min) < 0) { dist_min = dist; index = i; }
  }

  j = shl(st, index, (Word16)1);
  st->last_ener_pit = add(st, t_qua_ener[j],   pred_pit);
  st->last_ener_cod = add(st, t_qua_ener[j+1], pred_cod);

  /* Limit energies ->for transmission errors */

  if(sub(st, st->last_ener_pit, (Word16)6912) > 0) st->last_ener_pit = 6912;
									/* 6912 = 27 in Q8 */
  if(sub(st, st->last_ener_cod, (Word16)6400) > 0) st->last_ener_cod = 6400;
									/* 6400 = 25 in Q8 */

  /*---------------------------------------------------*
   *                                                   *
   *  Compute the quantized pitch gain.                *
   *                                                   *
   *     temp = 0.5 * ( last_ener_pit - ener_plt );    *
   *     temp = pow(2.0, temp);                        *
   *     if( temp > 1.2) temp = 1.2;                   *
   *     *gain_pit = temp;                             *
   *                                                   *
   *---------------------------------------------------*/

  L_tmp = Load_sh(st, st->last_ener_pit, (Word16)6);	/* last_ener_pit/2 in Q15 */
  L_tmp = sub_sh(st, L_tmp, ener_plt, (Word16)6);	/* - ener_plt/2    in Q15 */
  L_tmp = add_sh(st, L_tmp, (Word16)12, (Word16)15); /* to have gain in Q12   */
  L_extract(st, L_tmp, &exp, &frac);
  L_tmp = pow2(st, exp, frac);
  if( L_sub(st, L_tmp, (Word32)4915) > 0) L_tmp = 4915; /* 4915 = 1.2 in Q12  */
  *gain_pit = extract_l(L_tmp);


  /*---------------------------------------------------*
   *                                                   *
   *  Compute the innovative code gain.                *
   *                                                   *
   *     temp = 0.5 * ( last_ener_cod - ener_c );      *
   *     *gain_cod = pow(2.0, temp);                   *
   *                                                   *
   *---------------------------------------------------*/

  L_tmp = Load_sh(st, st->last_ener_cod, (Word16)6);	/* last_ener_cod/2 in Q15 */
  L_tmp = sub_sh(st, L_tmp, ener_c, (Word16)6);	/* - ener_c/2      in Q15 */
  L_extract(st, L_tmp, &exp, &frac);
  L_tmp = pow2(st, exp, frac);
  *gain_cod = extract_l(L_tmp);

  return index;
}


/**************************************************************************
*
*	ROUTINE				:	G_Code
*
*	DESCRIPTION			:	Compute the gain of innovative code
*
**************************************************************************
*
*	USAGE				:	G_Code(buffer_in1,buffer_in2,L_subfr)
*							(Routine_Name(input1,input2,input3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Code target
*						- Format : Word16
*
*		INPUT2			:	- Description : Filtered innovation code (with sign)
*							- Format : Word16
*
*	INPUT3			:	- Description : Length of subframe
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	Gain of innovation code
*
**************************************************************************/

Word16 G_Code(tetra_codec* st, Word16 xn2[], Word16 y2[], Word16 L_subfr)
{
   Word16 i;
   Word16 xy, yy, exp_xy, exp_yy, gain;
   Word32 s, L_tmp;

   st->Overflow = 0;

/* Compute scalar product <xn2[],y2[]> */

   s = 1;					/* Avoid case of all zeros */
   for(i=0; i<L_subfr; i++)
     s = L_mac0(st, s, xn2[i], y2[i]);

   exp_xy = norm_l(s);
   xy     = extract_h( L_shl(st, s, exp_xy) );

/* Compute scalar product <y2[],y2[]> */

   s = 1;					/* Avoid case of all zeros */
   for(i=0; i<L_subfr; i++)
     s = L_mac0(st, s, y2[i], y2[i]);

   exp_yy = norm_l(s);
   yy     = extract_h( L_shl(st, s, exp_yy) );


/* Test if Overflow */

   if(st->Overflow == 1)
   {
      /* Compute scalar product <xn2[],y2[]> */

      s = 1;				/* Avoid case of all zeros */
      for(i=0; i<L_subfr; i++)
      {
        L_tmp = L_mult0(xn2[i], y2[i]);
	L_tmp = L_shr(st, L_tmp, (Word16)6);
        s = L_add(st, s, L_tmp);
      }

      exp_xy = norm_l(s);
      xy     = extract_h( L_shl(st, s, exp_xy) );

      /* Compute scalar product <y2[],y2[]> */

      s = 1;				/* Avoid case of all zeros */
      for(i=0; i<L_subfr; i++)
      {
        L_tmp = L_mult0(y2[i], y2[i]);
	L_tmp = L_shr(st, L_tmp, (Word16)6);
        s = L_add(st, s, L_tmp);
      }

      exp_yy = norm_l(s);
      yy     = extract_h( L_shl(st, s, exp_yy) );

   }

/* If (xy < 0) gain = 0  */

   if( xy <= 0) return( (Word16) 0);

/* compute gain = xy/yy */

   xy = shr(st, xy, (Word16)1);			/* Be sure xy < yy */
   gain = div_s(st,  xy, yy);

   i = add(st, exp_xy, (Word16)2);		/* Denormalization of division */
   i = sub(st, i, exp_yy);

   gain = shr(st, gain, i);

   return(gain);
}


/**************************************************************************
*
*	ROUTINE				:	G_Pitch
*
*	DESCRIPTION			:	Compute the gain of pitch. Result in Q12
*							if (gain < 0)  gain =0
*								if (gain >1.2) gain =1.2
*
**************************************************************************
*
*	USAGE				:	G_Pitch(buffer_in1,buffer_in2,L_subfr)
*							(Routine_Name(input1,input2,input3))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Pitch target
*						- Format : Word16
*
*		INPUT2			:	- Description : Filtered adaptive codebook
*							- Format : Word16
*
*	INPUT3			:	- Description : Length of subframe
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	Gain of pitch lag in Q12 saturated to 1.2
*
**************************************************************************/

Word16 G_Pitch(tetra_codec* st, Word16 xn[], Word16 y1[], Word16 L_subfr)
{
   Word16 i;
   Word16 xy, yy, exp_xy, exp_yy, gain;
   Word32 s, L_tmp;

   st->Overflow = 0;

/* Compute scalar product <xn[],y1[]> */

   s = 1;					/* Avoid case of all zeros */
   for(i=0; i<L_subfr; i++)
     s = L_mac0(st, s, xn[i], y1[i]);

   exp_xy = norm_l(s);
   xy     = extract_h( L_shl(st, s, exp_xy) );

/* Compute scalar product <y1[],y1[]> */

   s = 1;					/* Avoid case of all zeros */
   for(i=0; i<L_subfr; i++)
     s = L_mac0(st, s, y1[i], y1[i]);

   exp_yy = norm_l(s);
   yy     = extract_h( L_shl(st, s, exp_yy) );


/* Test if Overflow */

   if(st->Overflow != 0)
   {
     /* Compute scalar product <xn[],y1[]> */

      s = 1;				/* Avoid case of all zeros */
      for(i=0; i<L_subfr; i++)
      {
        L_tmp = L_mult0(xn[i], y1[i]);
        L_tmp = L_shr(st, L_tmp, (Word16)6);
        s = L_add(st, s, L_tmp);
      }

      exp_xy = norm_l(s);
      xy     = extract_h( L_shl(st, s, exp_xy) );

      /* Compute scalar product <y1[],y1[]> */

      s = 1;				/* Avoid case of all zeros */
      for(i=0; i<L_subfr; i++)
      {
        L_tmp = L_mult0(y1[i], y1[i]);
        L_tmp = L_shr(st, L_tmp, (Word16)6);
        s = L_add(st, s, L_tmp);
      }

      exp_yy = norm_l(s);
      yy     = extract_h( L_shl(st, s, exp_yy) );

   }

/* If (xy < 4) gain = 0 */

   if( sub(st, xy, (Word16)4) < 0) return( (Word16) 0);

/* compute gain = xy/yy */

   xy = shr(st, xy, (Word16)1);			/* Be sure xy < yy */
   gain = div_s(st,  xy, yy);

   i = add(st, exp_xy, (Word16)2);		/* Denormalization of division */
   i = sub(st, i, exp_yy);

   gain = shr(st, gain, i);

/* if(gain >1.2) gain = 1.2  in Q12 */

   if( sub(st, gain, (Word16)4915) > 0) gain = 4915;

   return(gain);
}


/**************************************************************************
*
*	ROUTINE				:	Inter8_M1_3
*
*	DESCRIPTION			:	Fractional interpolation -1/3 with 8 coefficients
*
**************************************************************************
*
*	USAGE				:	Inter8_M1_3(buffer_in)
*							(Routine_Name(input1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Buffer containing the function
*								          to be interpolated
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Interpolated value
*						- Format : Word32
*
**************************************************************************/

static Word32 Inter8_M1_3(tetra_codec* st, Word16 x[])
{
  Word16 i;
  Word32 s;
  static const Word16 coef[8] = {-236,1050,-3572,12714,26674,-5217,1630,-384};

  s = 0;
  for(i=0; i<8; i++)
    s = L_mac0(st, s, x[i-4], coef[i]);

  return(s);
}


/**************************************************************************
*
*	ROUTINE				:	Inter8_1_3
*
*	DESCRIPTION			:	Fractional interpolation 1/3 with 8 coefficients
*
**************************************************************************
*
*	USAGE				:	Inter8_1_3(buffer_in)
*							(Routine_Name(input1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Buffer containing the function
*								          to be interpolated
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Interpolated value
*						- Format : Word32
*
**************************************************************************/

static Word32 Inter8_1_3(tetra_codec* st, Word16 x[])
{
  Word16 i;
  Word32 s;
  static const Word16 coef[8] = {-384,1630,-5217,26674,12714,-3572,1050,-236};

  s = 0;
  for(i=0; i<8; i++)
    s = L_mac0(st, s, x[i-3], coef[i]);

  return(s);
}


/**************************************************************************
*
*	ROUTINE				:	Inter32_M1_3
*
*	DESCRIPTION			:	Fractional interpolation -1/3 with 32 coefficients
*
**************************************************************************
*
*	USAGE				:	Inter32_M1_3(buffer_in)
*							(Routine_Name(input1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Buffer containing the function
*								          to be interpolated
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Interpolated value
*						- Format : Word16
*
*	COMMENTS			:	For long term prediction, it must be noted that
*							exc[-(T0-1/3)] corresponds to exc[-T0+1/3]
*
**************************************************************************/

static Word16 Inter32_M1_3(tetra_codec* st, Word16 x[])
{
  Word16 i;
  Word32 s;

  static const Word16 coef[32] = {
              -49,    66,   -96,   142,  -207,   294,  -407,   553,  -739,
       981, -1304,  1758, -2452,  3688, -6669, 27072, 13496, -5287,  3179,
     -2182,  1587, -1185,   893,  -672,   500,  -366,   263,  -183,   125,
       -84,    59,   -47 };

  s = 0;
  for(i=0; i<32; i++)
    s = L_mac0(st, s, x[i-15], coef[i]);

  s = L_add(st, s, s);
  i = tsc_round(st, s);
  return(i);
}


/**************************************************************************
*
*	ROUTINE				:	Inter32_1_3
*
*	DESCRIPTION			:	Fractional interpolation 1/3 with 32 coefficients
*
**************************************************************************
*
*	USAGE				:	Inter32_1_3(buffer_in)
*							(Routine_Name(input1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Buffer containing the function
*								          to be interpolated
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Interpolated value
*						- Format : Word16
*
*	COMMENTS			:	For long term prediction, it must be noted that
*							exc[-(T0+1/3)] corresponds to exc[-T0-1/3]
*
**************************************************************************/

static Word16 Inter32_1_3(tetra_codec* st, Word16 x[])
{
  Word16 i;
  Word32 s;

  static const Word16 coef[32] = {
      -47,    59,   -84,   125,  -183,   263,  -366,   500,  -672,   893,
    -1185,  1587, -2182,  3179, -5287, 13496, 27072, -6669,  3688, -2452,
     1758, -1304,   981,  -739,   553,  -407,   294,  -207,   142,   -96,
       66,   -49};

  s = 0;
  for(i=0; i<32; i++)
    s = L_mac0(st, s, x[i-16], coef[i]);

  s = L_add(st, s, s);
  i = tsc_round(st, s);
  return(i);
}


/**************************************************************************
*
*	ROUTINE				:	Lag_Max
*
*	DESCRIPTION			:	Find the lag that has maximum correlation with
*							the input buffer sig_dec[]
*
**************************************************************************
*
*	USAGE				:	Lag_Max(buffer_in1,buffer_in2,
*							L_frame, lag_max,lag_min,cor_max)
*							(Routine_Name(input1,input2,
*							input3,input4,input5,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description :	- Signal used to compute
*								the open loop pitch
*								- Buffer_in1[-142] to buffer_in1[-1]
*								should be known
*						- Format : Word16
*
*	INPUT2			:	- Description : Decimated signal (buffer sig_dec[])
*						- Format : Word16
*
*	INPUT3			:	- Description : Length of frame to compute pitch
*						- Format : Word16
*
*	INPUT4			:	- Description : Maximum lag
*						- Format : Word16
*
*		INPUT5			:	- Description : Minimum lag
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Maximum of normalized correlation
*								         of lag found
*							- Format : Word16
*
*	RETURNED VALUE		:	- Description : Lag found
*							- Format : Word16
*
**************************************************************************/

Word16 Lag_Max(tetra_codec* st, Word16 signal[], Word16 sig_dec[], Word16 L_frame,
               Word16 lag_max, Word16 lag_min, Word16 *cor_max)
{
  Word16  i, j;
  Word16  *p, *p1;
  Word32  max, t0;
  Word16  max_h, max_l, ener_h, ener_l;
  Word16  p_max = 0;

  max = MIN_32;


  for (i = lag_max; i >= lag_min; i--)
  {
    p  = sig_dec;
    p1 = &signal[-i];
    t0 = 0;

    for (j=0; j<L_frame; j+=2, p++, p1+=2)
      t0 = L_mac0(st, t0, *p, *p1);

    if (L_sub(st, t0,max) >= 0)
    {
      max    = t0;
      p_max = i;
    }
  }

  max = L_shr(st, max, (Word16)1);	/* for special double precision format */
  L_extract(st, max, &max_h, &max_l);

  /* compute energie */

  t0 = 0;
  p = &signal[-p_max];
  for(i=0; i<L_frame; i+=2, p+=2)
    t0 = L_mac0(st, t0, *p, *p);

  /* 1/sqrt(energie),    result in Q30 */

  t0 = inv_sqrt(st, t0);
  L_extract(st, t0, &ener_h, &ener_l);

  /* max = max/sqrt(energie)                  */
  /* This result will always be on 16 bits !! */

  t0 = mpy_32(st, max_h, max_l, ener_h, ener_l);
  *cor_max = extract_l(t0);

  return(p_max);
}


/**************************************************************************
*
*	ROUTINE				:	Norm_Corr
*
*	DESCRIPTION			:	Find the normalized correlation
*						(correlation between the target vector and
*						the filtered past excitation divided by the square root
*							of  energy of filtered excitation)
*
**************************************************************************
*
*	USAGE				:	Norm_Corr(buffer_in1,buffer_in2,buffer_in3,
*							L_subfr,t_min,t_max,corr_norm)
*							(Routine_Name(input1,input2,input3,
*							input4,input5,input6,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Excitation buffer
*						- Format : Word16
*
*	INPUT2			:	- Description : Target vector
*						- Format : Word16
*
*	INPUT3			:	- Description : Impulse response of synthesis
*							          and weighting filters
*						- Format : Word16 - Q12
*
*	INPUT4			:	- Description : Length of subframe
*						- Format : Word16
*
*		INPUT5			:	- Description : Minimum value of pitch lag
*							- Format : Word16
*
*		INPUT6			:	- Description : Maximum value of pitch lag
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Normalized correlation
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

static void Norm_Corr(tetra_codec* st, Word16 exc[], Word16 xn[], Word16 h[], Word16 L_subfr,
               Word16 t_min, Word16 t_max, Word16 corr_norm[])
{
 Word16 i,j,k;
 Word16 corr_h, corr_l, norm_hi, norm_lo;
 Word32 s;
 Word16 excf[80];			/* Usally dynamic allocation of (L_subfr) */



 k = - t_min;

 /* compute the filtered excitation for the first delay t_min */

 Convolve(st, &exc[k], h, excf, L_subfr);

 /* loop for every possible period */

 for (i = t_min; i <= t_max; ++i)
 {

   /* Compute correlation between xn[] and excf[] */

   s = 0;
   for (j = 0; j < L_subfr; j++)
     s = L_mac0(st, s, xn[j], excf[j]);

   s = L_shr(st, s, (Word16)1);		/* Special double precision format */
   L_extract(st, s, &corr_h, &corr_l);


   /* Compute 1/sqrt(energie of excf[]) */

   s = 0;
   for (j = 0; j < L_subfr; j++)
     s = L_mac0(st, s, excf[j], excf[j]);

   s = inv_sqrt(st, s);			/* Result in Q30 */
   L_extract(st, s, &norm_hi, &norm_lo);


   /* Normalize correlation = correlation * (1/sqrt(energie)) */

   s = mpy_32(st, corr_h, corr_l, norm_hi, norm_lo);

   corr_norm[i] = extract_l(s);	/* Result is on 16 bits */


   /* modify the filtered excitation excf[] for the next iteration */

   if( sub(st, i, t_max) != 0)
   {
     k--;
     for (j = L_subfr-1; j > 0; j--)
     {
        s = L_mult0(exc[k], h[j]);
        s = add_sh(st, s, excf[j-1], (Word16)12);
        excf[j] = store_hi(st, s, (Word16)4);
     }
      excf[0] = exc[k];
   }
 }
 return;
}


/**************************************************************************
*
*	ROUTINE				:	Pitch_Fr
*
*	DESCRIPTION			:	Find the pitch period with 1/3 subsample resolution
*
**************************************************************************
*
*	USAGE				:	Pitch_Fr(buffer_in1,buffer_in2,buffer_in3,
*							L_subfr,t0_min,t0_max,i_subfr,pit_frac)
*							(Routine_Name(input1,input2,input3,
*							input4,input5,input6,input7,output1))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description : Excitation buffer
*						- Format : Word16
*
*	INPUT2			:	- Description : Target vector
*						- Format : Word16
*
*	INPUT3			:	- Description : Impulse response of synthesis
*							          and weighting filters
*						- Format : Word16 - Q12
*
*	INPUT4			:	- Description : Length of subframe
*						- Format : Word16
*
*		INPUT5			:	- Description : Minimum value in the searched range
*							- Format : Word16
*
*		INPUT6			:	- Description : Maximum value in the searched range
*							- Format : Word16
*
*		INPUT7			:	- Description : Indicator for first subframe
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Chosen fraction : (-1, 0, 1) / 3
*							- Format : Word16
*
*	RETURNED VALUE		:	- Description : Pitch period
*							- Format : Word16
*
**************************************************************************/

       /* TSC_Lg_inter = Length for fractionnal interpolation = nb.coeff/2 */
#define TSC_Lg_inter 4

Word16 Pitch_Fr(tetra_codec* st, Word16 exc[], Word16 xn[], Word16 h[], Word16 L_subfr,
                Word16 t0_min, Word16 t0_max, Word16 i_subfr,
                Word16 *pit_frac)
{
  Word16 i;
  Word16 t_min, t_max;
  Word16 max, lag, frac;
  Word16 *corr;
  Word32 corr_int, L_max;
  Word16 corr_v[40];	/* Total length = t0_max-t0_min+1+2*TSC_Lg_inter */

  /* Find interval to compute normalized correlation */

  t_min = sub(st, t0_min, (Word16)TSC_Lg_inter);
  t_max = add(st, t0_max, (Word16)TSC_Lg_inter);

  corr = &corr_v[-t_min];

  /* Compute normalized correlation between target and filtered excitation */

  Norm_Corr(st, exc, xn, h, L_subfr, t_min, t_max, corr);

  /* Find integer pitch */

  max = corr[t0_min];
  lag = t0_min;

  for(i= t0_min+1; i<=t0_max; i++)
  {
    if( sub(st, corr[i],max)>= 0)
    {
      max = corr[i];
      lag = i;
    }
  }
  /* If first subframe and lag > 84 do not search fractionnal pitch */

  if( (i_subfr == 0) && (sub(st, lag, (Word16)84) > 0) )
  {
    *pit_frac = 0;
    return(lag);
  }

  /* Test the fractions around T0 and choose the one which maximizes   */
  /* the interpolated normalized correlation.                          */

  frac  = 0;
  L_max = Load_sh(st, max, (Word16)15);

  /* test +1/3 */

  corr_int = Inter8_1_3(st, &corr[lag]);
  if(L_sub(st, corr_int, L_max) >= 0)
    { L_max=corr_int; frac= 1;}

  /* test +2/3  = 1-1/3 */

  corr_int = Inter8_M1_3(st, &corr[lag+1]);
  if(L_sub(st, corr_int, L_max) >= 0)
    { L_max=corr_int; frac= 2;}

  /* test -2/3 = 1/3 -1 */

  corr_int = Inter8_1_3(st, &corr[lag-1]);
  if(L_sub(st, corr_int, L_max) >= 0)
    { L_max=corr_int; frac= -2;}

  /* test -1/3 */

  corr_int = Inter8_M1_3(st, &corr[lag]);
  if(L_sub(st, corr_int, L_max) >= 0)
    { L_max=corr_int; frac= -1;}

  if(sub(st, frac, (Word16)2) == 0)
    { frac = -1; lag = add(st, lag, (Word16)1);};
  if(sub(st, frac, (Word16)-2) == 0)
    { frac = 1; lag = sub(st, lag, (Word16)1);};
  *pit_frac = frac;
  return(lag);
}

/**************************************************************************
*
*	ROUTINE				:	Pitch_Ol_Dec
*
*	DESCRIPTION			:	Compute the open loop pitch lag
*							(includes decimation of signal)
*
**************************************************************************
*
*	USAGE				:	Pitch_Ol_Dec(buffer_in,L_frame)
*							(Routine_Name(input1,input2))
*
*	INPUT ARGUMENT(S)		:
*
*	INPUT1			:	- Description :	- Signal used to compute
*								the open loop pitch
*								- Buffer_in[-142] to buffer_in[-1]
*								should be known
*						- Format : Word16
*
*	INPUT2			:	- Description : Length of frame to compute pitch
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	- Description : Open pitch lag as found
*							- Format : Word16
*
**************************************************************************/

/* Minimum and maximum pitch lag */

//#define TSC_pit_min 20
#define TSC_pit_max_2 142

/* Threshold to favor small pitch = 0.85 in Q15 */

#define TSC_seuil 27856

Word16 Pitch_Ol_Dec(tetra_codec* st, Word16 signal[], Word16 L_frame)
{
  Word16  i, j;
  Word16  max1, max2, max3;
  Word16  p_max1, p_max2, p_max3;
  Word32  t0;

  /* Decimated signal                                       */
  /* Can be allocated with memory allocation of(L_frame/2)  */

  Word16  sig_dec[120];


  /*--------------------------------------------------------*
   *  Verification for risk of overflow.                    *
   *--------------------------------------------------------*/

   st->Overflow = 0;
   t0 = 0;
   for(i= -TSC_pit_max_2; i< L_frame; i+=2)
     t0 = L_mac0(st, t0, signal[i], signal[i]);



  /*--------------------------------------------------------*
   * Decimation of signal and scaling.                      *
   *                                                        *
   *   if Overflow        -> sig_dec[i] = signal[i*2]>>6    *
   *   else if t0 < 1<<22 -> sig_dec[i] = signal[i*2]<<4    *
   *   else               -> sig_dec[i] = signal[i*2]       *
   *--------------------------------------------------------*/

   if(st->Overflow == 1)
   {
     for(i=0, j=0; i<L_frame; i+=2, j++)
       sig_dec[j] = shr(st, signal[i], (Word16)6);
   }
   else if (L_sub(st, t0, L_shl(st, (Word32)1, (Word16)22)) < 0 )
   {
     for(i=0, j=0; i<L_frame; i+=2, j++)
       sig_dec[j] = shl(st, signal[i], (Word16)4);
   }
   else
   {
     for(i=0, j=0; i<L_frame; i+=2, j++)
       sig_dec[j] = signal[i];
   }

  /*--------------------------------------------------------------------*
   *  The pitch lag search is divided in three sections.                *
   *  Each section cannot have a pitch multiple.                        *
   *  A maximum is found for each section.                              *
   *  The maximum of each section is compared to favor small lag.       *
   *                                                                    *
   *  First section:  lag delay = TSC_pit_max_2 to 80                         *
   *  Second section: lag delay = 79 to 40                              *
   *  Third section:  lag delay = 39 to 20                              *
   *--------------------------------------------------------------------*/

   p_max1 = Lag_Max(st, signal, sig_dec, L_frame, (Word16)TSC_pit_max_2,
									(Word16)80 , &max1);
   p_max2 = Lag_Max(st, signal, sig_dec, L_frame, (Word16)79     ,
									(Word16)40 , &max2);
   p_max3 = Lag_Max(st, signal, sig_dec, L_frame, (Word16)39     ,
									(Word16)20 , &max3);

  /*--------------------------------------------------------------------*
   * Compare the 3 sections maximum, and favor small lag.               *
   *--------------------------------------------------------------------*/

  if( sub(st, mult(st, (Word16)max1, (Word16)TSC_seuil), (Word16)max2)  < 0)
    { max1 = max2;  p_max1 = p_max2;}

  if( sub(st, mult(st, (Word16)max1, (Word16)TSC_seuil), (Word16)max3)  < 0)
    p_max1 = p_max3;


  return (p_max1);
}


/**************************************************************************
*
*	ROUTINE				:	Post_Process
*
*	DESCRIPTION			:	Post-processing of output speech
*							Multiplication by two of output speech
*							with saturation
*
**************************************************************************
*
*	USAGE				:	Post_Process(st, buffer,Length)
*							(Routine_Name(arg1,input2))
*
*	INPUT ARGUMENT(S)		:
*
*	ARG1				:	- Description : Input speech signal buffer
*						- Format : Word16
*
*		INPUT2			:	- Description : Length of signal
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		ARG1				:	- Description : Output speech signal buffer
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Post_Process(tetra_codec* st, Word16 signal[], Word16 lg)
{
  Word16 i;

  for(i=0; i<lg; i++)
    signal[i] = add(st, signal[i], signal[i]);
}


/**************************************************************************
*
*	ROUTINE				:	Pred_Lt
*
*	DESCRIPTION			:	Compute the result of long term prediction with
*							fractional interpolation
*
**************************************************************************
*
*	USAGE				:	Pred_Lt(buffer,T0,frac,L_subfr)
*							(Routine_Name(arg1,input2,input3,input4))
*
*	INPUT ARGUMENT(S)		:
*
*	ARG1				:	- Description : Excitation vector
*						- Format : Word16
*
*	INPUT2			:	- Description : Pitch lag
*						- Format : Word16
*
*	INPUT3			:	- Description : Fraction of pitch lag : (-1, 0, 1)  / 3
*						- Format : Word16 - Q12
*
*	INPUT4			:	- Description : Length of subframe
*						- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		ARG1				:	- Description : Interpolated signal contained in
*								         buffer[0..L_subfr-1]
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Pred_Lt(tetra_codec* st, Word16 exc[], Word16 T0, Word16 frac, Word16 L_subfr)
{
   Word16 i;

   if(frac == 0)
   {
      for (i = 0; i < L_subfr; i++)
         exc[i] = exc[i-T0];
   }
   else if( sub(st, frac, (Word16)1) == 0)
   {
      for (i = 0; i < L_subfr; i++)
         exc[i] = Inter32_1_3(st, &exc[i-T0]);
   }
   else if( sub(st, frac, (Word16)-1) == 0)
   {
      for (i = 0; i < L_subfr; i++)
         exc[i] = Inter32_M1_3(st, &exc[i-T0]);
   }

   return;
}


/**************************************************************************
*
*	ROUTINE				:	Pre_Process
*							(include Init_Pre_Process)
*
*	DESCRIPTION			:	Preprocessing of input speech
*							- Offset compensation
*								- Divide input by two
*
**************************************************************************
*
*	USAGE				:	Pre_Process(buffer,Length)
*							(Routine_Name(arg1,input2))
*
*	INPUT ARGUMENT(S)		:
*
*	ARG1				:	- Description : Input speech signal buffer
*						- Format : Word16
*
*		INPUT2			:	- Description : Length of signal
*							- Format : Word16
*
*	OUTPUT ARGUMENT(S)		:
*
*		ARG1				:	- Description : Output speech signal buffer
*							- Format : Word16
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	Algorithm :
*
*						y[i] = x[i]/2 - x[i-1]/2 + alpha * y[i-1]
*
*  						alpha = 32735 / 2**15
*  						y[i-1] is keep in double precision
*
*  						It is the same algorithm than in GSM except that :
*  						- Input is divided by two in the filtering process
*  						- No rounding
*
**************************************************************************/

/* Static values to be preserved between calls */

//static Word16 y_hi, y_lo, x0;

/* Initialization of static values */

void Init_Pre_Process(tetra_codec* st)
{
  st->y_hi = 0;
  st->y_lo = 0;
  st->x0   = 0;
}


/* Offset compensation and divide by 2 */

void Pre_Process(tetra_codec* st, Word16 signal[], Word16 lg)
{
  Word16 i, x1;
  Word32 L_tmp;

  for(i=0; i<lg; i++)
  {
     x1 = st->x0;
     st->x0 = signal[i];

     L_tmp     = Load_sh(st, st->x0, (Word16)15);
     L_tmp     = sub_sh(st, L_tmp, x1, (Word16)15);
     L_tmp     = L_mac(st, L_tmp, st->y_hi, (Word16)32735);
     L_tmp     = add_sh(st, L_tmp, mult(st, st->y_lo, (Word16)32735), (Word16)1);
     signal[i] = extract_h(L_tmp);
     st->y_hi      = extract_h(L_tmp);
     st->y_lo      = extract_l(sub_sh(st, L_shr(st, L_tmp,(Word16)1 ), st->y_hi, (Word16)15));
  }
  return;
}

/**************************************************************************
*
*	ROUTINE				:	Prm2bits_Tetra
*
*	DESCRIPTION			:	Convert the encoder parameter vector into
*							a vector of serial bits
*
**************************************************************************
*
*	USAGE				:	Prm2bits_Tetra(st, buffer_in,buffer_out)
*							(Routine_Name(input1,output1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Encoded parameters
*								         (23 parameters)
*							- Format : Word16 - .. * 16 bit-samples
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Serial bits (137 + bfi)
*							- Format : Word16 - .. * 16 bit-samples
*
*	RETURNED VALUE		:	None
*
*	COMMENTS			:	The transmitted parameters are :
*
*						- LSP	1st codebook 			8 bits
*							2nd codebook			9 bits
*							3nd codebook			9 bits
*
*						- for the 4 subframes :
*							pitch delay			8 bits (first)
*											5 bits (others)
*							codebook index		14 bits
*							pulse global sign		1 bit
*							pulse shift			1 bit
*							pitch and innovation gains	6 bits
*
**************************************************************************/

#define TSC_PRM_NO    23

void Prm2bits_Tetra(tetra_codec* st, Word16 prm[], Word16 bits[])
{
  Word16 i;
  static Word16 bitno[TSC_PRM_NO] = {8, 9, 9,            /* split VQ LSP  */
                                 8, 14, 1, 1, 6,     /* subframe 1    */
                                 5, 14, 1, 1, 6,     /* subframe 2    */
                                 5, 14, 1, 1, 6,     /* subframe 3    */
                                 5, 14, 1, 1, 6};    /* subframe 4    */

  *bits++ = 0;	/* bit[0] = 0, at receiver this bits indicate BFI */

  for (i = 0; i < TSC_PRM_NO; i++)
  {
    int2bin(st, prm[i], bitno[i], bits);
    bits += bitno[i];
  }
  return;
}





//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- SDEC_TET.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------





/************************************************************************
*
*	FILENAME		:	sdec_tet.c
*
*	DESCRIPTION		:	Main routines for speech source decoding
*
************************************************************************
*
*	SUB-ROUTINES	:	- Init_Decod_Tetra()
*					- Decod_Tetra()
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*
************************************************************************/

//#include "source.h"

/*--------------------------------------------------------*
 *       Decoder constants parameters.                    *
 *                                                        *
 *   TSC_L_frame     : Frame size.                            *
 *   TSC_L_subfr     : Sub-frame size.                        *
 *   TSC_p           : LPC order.                             *
 *   TSC_pp1         : LPC order+1                            *
 *   TSC_pit_min     : Minimum pitch lag.                     *
 *   TSC_pit_max     : Maximum pitch lag.                     *
 *   TSC_L_inter     : Length of filter for interpolation     *
 *   TSC_parm_size   : Lenght of vector parm[]                *
 *--------------------------------------------------------*/

//#define  TSC_L_frame  (Word16)240
//#define  TSC_L_subfr  (Word16)60
//#define  TSC_p        (Word16)10
//#define  TSC_pp1      (Word16)11
#define  TSC_pit_min  (Word16)20
//#define  TSC_pit_max  (Word16)143
//#define  TSC_L_inter  (Word16)15
//#define  TSC_parm_size (Word16)23


/*--------------------------------------------------------*
 *   LPC bandwidth expansion factors for noise filter.    *
 *      In Q15 = 0.75, 0.85                               *
 *--------------------------------------------------------*/

#define TSC_gamma3  (Word16)24576
#define TSC_gamma4  (Word16)27853


/*--------------------------------------------------------*
 *         Static memory allocation.                      *
 *--------------------------------------------------------*/

        /* Excitation vector */

//static Word16 old_exc[TSC_L_frame+TSC_pit_max+TSC_L_inter];
//static Word16 *exc;

        /* Spectral expansion factors */

//static Word16 F_gamma3[TSC_p];
//static Word16 F_gamma4[TSC_p];

        /* Lsp (Line spectral pairs in the cosine domain) */

//static Word16 lspold[TSC_p]={
//              30000, 26000, 21000, 15000, 8000, 0,
//		  -8000,-15000,-21000,
//			-26000};
//static Word16 lspnew[TSC_p];

	  /* Initial lsp values used after each time */
        /* a reset is executed */

//static Word16 lspold_init[TSC_p]={
//              30000, 26000, 21000, 15000, 8000, 0,
//		  -8000,-15000,-21000,-26000};

        /* Filter's memory */

//static Word16 mem_syn[TSC_p];

        /* Default parameters */

//static Word16 old_parm[TSC_parm_size], old_T0;

       /* Global definition */

//Word16 last_ener_cod;
//Word16 last_ener_pit;


/**************************************************************************
*
*	ROUTINE				:	Init_Decod_Tetra
*
*	DESCRIPTION			:	Initialization of variables for the speech decoder
*
**************************************************************************
*
*	USAGE				:	Init_Decod_Tetra()
*
*	INPUT ARGUMENT(S)		:	None
*
*	OUTPUT ARGUMENT(S)		:	None
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Init_Decod_Tetra(tetra_codec* st)
{
  Word16 i;

  st->old_T0 = 60;
  for(i=0; i<23; i++)
     st->old_parm[i] = 0;

  /* Initialize static pointer */

  st->exc    = st->old_exc + TSC_pit_max + TSC_L_inter;

  /* Initialize global variables */

  st->last_ener_cod = 0;
  st->last_ener_pit = 0;

  /* Static vectors to zero */

  for(i=0; i<TSC_pit_max + TSC_L_inter; i++)
    st->old_exc[i] = 0;

  for(i=0; i<TSC_p; i++)
    st->mem_syn[i] = 0;


  /* Initialisation of lsp values for first */
  /* frame lsp interpolation */

  for(i=0; i<TSC_p; i++)
    st->lspold[i] = lspold_init[i];


  /* Compute LPC spectral expansion factors */

  Fac_Pond(st, TSC_gamma3, st->F_gamma3);
  Fac_Pond(st, TSC_gamma4, st->F_gamma4);

 return;
}


/**************************************************************************
*
*	ROUTINE				:	Decod_Tetra
*
*	DESCRIPTION			:	Main speech decoder function
*
**************************************************************************
*
*	USAGE				:	Decod_Tetra(parm,synth)
*							(Routine_Name(input1,output1))
*
*	INPUT ARGUMENT(S)		:
*
*		INPUT1			:	- Description : Synthesis parameters
*							- Format : 24 * 16 bit-samples
*
*	OUTPUT ARGUMENT(S)		:
*
*		OUTPUT1			:	- Description : Synthesis
*							- Format : 240 * 16 bit-samples
*
*	RETURNED VALUE		:	None
*
**************************************************************************/

void Decod_Tetra(tetra_codec* st, Word16 parm[], Word16 synth[])
{
  /* LPC coefficients */

  Word16 A_t[(TSC_pp1)*4];		/* A(z) unquantized for the 4 subframes */
  Word16 Ap3[TSC_pp1];		/* A(z) with spectral expansion         */
  Word16 Ap4[TSC_pp1];		/* A(z) with spectral expansion         */
  Word16 *A;			/* Pointer on A_t                       */

  /* Other vectors */

  Word16 zero_F[TSC_L_subfr+64],  *F;
  Word16 code[TSC_L_subfr+4];

  /* Scalars */

  Word16 i, i_subfr;
  Word16 T0 = 0, T0_min = 0, T0_max, T0_frac = 0;
  Word16 gain_pit, gain_code, index;
  Word16 sign_code, shift_code;
  Word16 bfi, temp;
  Word32 L_temp;

  /* Initialization of F */

  F  = &zero_F[64];
  for(i=0; i<64; i++)
   zero_F[i] = 0;

  /* Test bfi */

  bfi = *parm++;

  if(bfi == 0)
  {
    D_Lsp334(st, &parm[0], st->lspnew, st->lspold);	/* lsp decoding   */

    for(i=0; i< TSC_parm_size; i++)		/* keep parm[] as old_parm */
      st->old_parm[i] = parm[i];
  }
  else
  {
    for(i=1; i<TSC_p; i++)
      st->lspnew[i] = st->lspold[i];

    for(i=0; i< TSC_parm_size; i++)		/* use old parm[] */
      parm[i] = st->old_parm[i];
  }

  parm += 3;			/* Advance synthesis parameters pointer */

  /* Interpolation of LPC for the 4 subframes */

  Int_Lpc4(st, st->lspold,   st->lspnew,   A_t);

  /* update the LSPs for the next frame */

  for(i=0; i<TSC_p; i++)
    st->lspold[i]   = st->lspnew[i];

/*------------------------------------------------------------------------*
 *          Loop for every subframe in the analysis frame                 *
 *------------------------------------------------------------------------*
 * The subframe size is TSC_L_subfr and the loop is repeated TSC_L_frame/TSC_L_subfr
 *
 *  times                                                                 *
 *     - decode the pitch delay                                           *
 *     - decode algebraic code                                            *
 *     - decode pitch and codebook gains                                  *
 *     - find the excitation and compute synthesis speech                 *
 *------------------------------------------------------------------------*/

  A = A_t;				/* pointer to interpolated LPC parameters */

  for (i_subfr = 0; i_subfr < TSC_L_frame; i_subfr += TSC_L_subfr)
  {

    index = *parm++;				/* pitch index */

    if (i_subfr == 0)				/* if first subframe */
    {
      if (bfi == 0)
      {						/* if bfi == 0 decode pitch */
         if (index < 197)
         {
           /* T0 = (index+2)/3 + 19; T0_frac = index - T0*3 + 58; */

           i = add(st, index, (Word16)2);
           i = mult(st, i, (Word16)10923);	/* 10923 = 1/3 in Q15 */
           T0 = add(st, i, (Word16)19);

           i = add(st, T0, add(st, T0, T0) );	/* T0*3 */
           i = sub(st, (Word16)58, (Word16)i);
           T0_frac = add(st, index, (Word16)i);
         }
         else
         {
           T0 = sub(st, index, (Word16)112);
           T0_frac = 0;
         }
      }
      else   /* bfi == 1 */
      {
        T0 = st->old_T0;
        T0_frac = 0;
      }


      /* find T0_min and T0_max for other subframes */

      T0_min = sub(st, T0, (Word16)5);
      if (T0_min < TSC_pit_min) T0_min = TSC_pit_min;
      T0_max = add(st, T0_min, (Word16)9);
      if (T0_max > TSC_pit_max)
      {
        T0_max = TSC_pit_max;
        T0_min = sub(st, T0_max, (Word16)9);
      }
    }

    else  /* other subframes */

    {
      if (bfi == 0)				/* if bfi == 0 decode pitch */
      {
         /* T0 = (index+2)/3 - 1 + T0_min; */

         i = add(st, index, (Word16)2);
         i = mult(st, i, (Word16)10923);	/* 10923 = 1/3 in Q15 */
         i = sub(st, i, (Word16)1);
         T0 = add(st, T0_min, i);

         /* T0_frac = index - 2 - i*3; */

         i = add(st, i, add(st, i,i) );		/* i*3 */
         T0_frac = sub(st,  index , add(st, i, (Word16)2) );
      }
    }

   /*-------------------------------------------------*
    * - Find the adaptive codebook vector.            *
    *-------------------------------------------------*/

    Pred_Lt(st, &st->exc[i_subfr], T0, T0_frac, TSC_L_subfr);

   /*-----------------------------------------------------*
    * - Compute noise filter F[].                         *
    * - Decode codebook sign and index.                   *
    * - Find the algebraic codeword.                      *
    *-----------------------------------------------------*/

    Pond_Ai(st, A, st->F_gamma3, Ap3);
    Pond_Ai(st, A, st->F_gamma4, Ap4);

    for (i = 0;   i <= TSC_p;      i++) F[i] = Ap3[i];
    for (i = TSC_pp1; i < TSC_L_subfr; i++) F[i] = 0;

    Syn_Filt(st, Ap4, F, F, TSC_L_subfr, &F[TSC_pp1], (Word16)0);

    /* Introduce pitch contribution with fixed gain of 0.8 to F[] */

    for (i = T0; i < TSC_L_subfr; i++)
    {
      temp = mult(st, F[i-T0], (Word16)26216);
      F[i] = add(st, F[i], temp);
    }

    index = *parm++;
    sign_code  = *parm++;
    shift_code = *parm++;

    D_D4i60(st, index, sign_code, shift_code, F, code);


   /*-------------------------------------------------*
    * - Decode pitch and codebook gains.              *
    *-------------------------------------------------*/

    index = *parm++;        /* index of energy VQ */

    Dec_Ener(st, index,bfi,A,&st->exc[i_subfr],code, TSC_L_subfr, &gain_pit, &gain_code);

   /*-------------------------------------------------------*
    * - Find the total excitation.                          *
    * - Find synthesis speech corresponding to exc[].       *
    *-------------------------------------------------------*/

    for (i = 0; i < TSC_L_subfr;  i++)
    {
      /* exc[i] = gain_pit*exc[i] + gain_code*code[i]; */
      /* exc[i]  in Q0   gain_pit in Q12               */
      /* code[i] in Q12  gain_cod in Q0                */

      L_temp = L_mult0(st->exc[i+i_subfr], gain_pit);
      L_temp = L_mac0(st, L_temp, code[i], gain_code);
      st->exc[i+i_subfr] = (Word16)L_shr_r(st, L_temp, (Word16)12);
    }

    Syn_Filt(st, A, &st->exc[i_subfr], &synth[i_subfr], TSC_L_subfr, st->mem_syn, (Word16)1);

    A  += TSC_pp1;    /* interpolated LPC parameters for next subframe */
  }

 /*--------------------------------------------------*
  * Update signal for next frame.                    *
  * -> shift to the left by TSC_L_frame  exc[]           *
  *--------------------------------------------------*/

  for(i=0; i<TSC_pit_max+TSC_L_inter; i++)
    st->old_exc[i] = st->old_exc[i+TSC_L_frame];

  st->old_T0 = T0;

  return;
}




//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-- TETRA_OP.C ----------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



/***********************************************************************
*
*	FILENAME		:	tetra_op.c
*
*	DESCRIPTION		:	TETRA  basic operators
*
************************************************************************
*
*	FUNCTIONS		:	- abs_s()
*					- add(st, )
*					- div_s(st, )
*					- extract_h()
*					- extract_l()
*					- L_abs()
*					- L_add(st, )
*					- L_deposit_h()
*					- L_deposit_l()
*					- L_mac()
*					- L_mac0()
*					- L_msu()
*					- L_msu0(st, )
*					- L_mult(st, )
*					- L_mult0()
*					- L_negate()
*					- L_shl(st, )
*					- L_shr(st, )
*					- L_shr_r(st, )
*					- L_sub(st, )
*					- mult(st, )
*					- mult_r(st, )
*					- negate()
*					- norm_l()
*					- norm_s()
*					- round()
*					- sature(st, )
*					- shl(st, )
*					- shr(st, )
*					- sub(st, )
*
*	COMMENTS		:	Only the operators used in the actual version
*					of the TETRA codec are included in this file
*
************************************************************************
*
*	INCLUDED FILES	:	source.h
*					stdio.h
*					stdlib.h
*
************************************************************************/


//#include <stdio.h>
//#include <stdlib.h>
//#include "source.h"

/*-----------------*
 * Local Functions *
 *-----------------*/

static Word16 sature(tetra_codec* st, Word32 L_var1);

/*-----------------------*
 * Constants and Globals *
 *-----------------------*/

//Flag Overflow =0;
//Flag Carry =0;


/************************************************************************
*
*	Function Name : abs_s
*
*	Purpose :
*
*		Absolute value of var1; abs_s(-32768) = 32767.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 abs_s(Word16 var1)
  {
   Word16 var_out;

   if (var1 == (Word16)0X8000 )
     {
      var_out = MAX_16;
     }
   else
     {
      if (var1 < 0)
 {
  var_out = -var1;
 }
      else
 {
  var_out = var1;
 }
     }

   return(var_out);
  }

/************************************************************************
*
*	Function Name : add
*
*	Purpose :
*
*		Performs the addition (var1+var2) with overflow control and saturation;|
*		the 16 bit result is set at +32767 when overflow occurs or at -32768
*		when underflow occurs.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 add(tetra_codec* st, Word16 var1,Word16 var2)
  {
   Word16 var_out;
   Word32 L_somme;

   L_somme = (Word32) var1 + var2;
   var_out = sature(st, L_somme);
   return(var_out);
  }


/************************************************************************
*
*	Function Name : div_s
*
*	Purpose :
*
*		Produces a result which is the fractionnal integer division of var1 by
*		var2; var1 and var2 must be positive and var2 must be greater or equal
*		to var1; the result is positive (leading bit equal to 0) and truncated
*		to 16 bits.
*		If var1 = var2 then div(var1,var2) = 32767.
*
*	Complexity Weight : 18
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= var1 <= var2 and var2 != 0.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : var1 <= var2 <= 0x0000 7fff and var2 != 0.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= var_out <= 0x0000 7fff.
*			It's a Q15 value (point between b15 and b14).
*
************************************************************************/

Word16 div_s(tetra_codec* st, Word16 var1, Word16 var2)
  {
   Word16 var_out = 0;
   Word16 iteration;
   Word32 L_num;
   Word32 L_denom;

   if ((var1 > var2) || (var1 < 0) || (var2 < 0))
     {
      //printf("Division Error\n");
      //exit(0);
	  return 0;
     }

   if (var2 == 0)
     {
      //printf("Division by 0, Fatal error \n");
      //exit(0);
	  return 0;
     }

   if (var1 == 0)
     {
      var_out = 0;
     }
   else
     {
      if (var1 == var2)
 {
  var_out = MAX_16;
 }
      else
 {
  L_num = L_deposit_l(var1);
  L_denom = L_deposit_l(var2);

  for(iteration=0;iteration<15;iteration++)
    {
     var_out <<=1;
     L_num <<= 1;

     if (L_num >= L_denom)
       {
        L_num = L_sub(st, L_num,L_denom);
        var_out = add(st,  var_out,(Word16)1 );
       }
    }
 }
     }

   return(var_out);
  }

/************************************************************************
*
*	Function Name : extract_h
*
*	Purpose :
*
*		Return the 16 MSB of L_var1.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32 ) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 extract_h(Word32 L_var1)
  {
   Word16 var_out;

   var_out = (Word16) (L_var1 >> 16);
   return(var_out);
  }

/************************************************************************
*
*	Function Name : extract_l
*
*	Purpose :
*
*		Return the 16 LSB of L_var1.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32 ) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 extract_l(Word32 L_var1)
  {
   Word16 var_out;

   var_out = (Word16) L_var1;
   return(var_out);
  }

/************************************************************************
*
*	Function Name : L_abs
*
*	Purpose :
*
*		Absolute value of L_var1; Saturate in case where the input is -214783648
*
*	Complexity Weight : 3
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x0000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_abs(Word32 L_var1)
  {
   Word32 L_var_out;

   if (L_var1 == MIN_32)
     {
      L_var_out = MAX_32;
     }
   else
     {
      if (L_var1 < 0)
 {
  L_var_out = -L_var1;
 }
      else
 {
  L_var_out = L_var1;
 }
     }

   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_add
*
*	Purpose :
*
*		32 bits addition of the two 32 bits variables (L_var1+L_var2) with
*		overflow control and saturation; the result is set at +214783647 when
*		overflow occurs or at -214783648 when underflow occurs.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*		L_var2
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var2 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_add(tetra_codec* st, Word32 L_var1, Word32 L_var2)
  {
   Word32 L_var_out;

   L_var_out = L_var1 + L_var2;

   if (((L_var1 ^ L_var2) & MIN_32) == 0)
     {
      if ((L_var_out ^ L_var1) & MIN_32)
        {
  L_var_out = (L_var1 < 0) ? MIN_32 : MAX_32;
         st->Overflow = 1;
        }
     }
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_deposit_h
*
*	Purpose :
*
*		Deposit the 16 bit var1 into the 16 MS bits of the 32 bit output. The
*		16 LS bits of the output are zeroed.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= var_out <= 0x7fff 0000.
*
************************************************************************/

Word32 L_deposit_h(Word16 var1)
  {
   Word32 L_var_out;

   L_var_out = (Word32) var1 << 16;
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_deposit_l
*
*	Purpose :
*
*		Deposit the 16 bit var1 into the 16 LS bits of the 32 bit output. The
*		16 MS bits of the output are sign extended.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0xFFFF 8000 <= L_var_out <= 0x0000 7fff.
*
************************************************************************/

Word32 L_deposit_l(Word16 var1)
  {
   Word32 L_var_out;

   L_var_out = (Word32) var1;
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_mac
*
*	Purpose :
*
*		Multiply var1 by var2 and shift the result left by 1. Add the 32 bit
*		result to L_var3 with saturation, return a 32 bit result:
*			L_mac(L_var3,var1,var2) = L_add(st, L_var3,(L_mult(st, var1,var2)).
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var3
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var3 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_mac(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2)
  {
   Word32 L_var_out;
   Word32 L_produit;

   L_produit = L_mult(st, var1,var2);
   L_var_out = L_add(st, L_var3,L_produit);
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_mac0
*
*	Purpose :
*
*		Multiply var1 by var2. Add the 32 bit result to L_var3.
*		Control saturation and set overflow_flag.
*		Same as L_mac without left shift.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var3
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var3 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_mac0(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2)
  {
   Word32 L_var_out;
   Word32 L_produit;

   L_produit = L_mult0(var1,var2);
   L_var_out = L_add(st, L_var3,L_produit);
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_msu
*
*	Purpose :
*
*		Multiply var1 by var2 and shift the result left by 1. Subtract the 32
*		bit result to L_var3 with saturation, return a 32 bit result:
*			L_msu(L_var3,var1,var2) = L_sub(st, L_var3,(L_mult(st, var1,var2)).
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var3
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var3 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_msu(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2)
  {
   Word32 L_var_out;
   Word32 L_produit;

   L_produit = L_mult(st, var1,var2);
   L_var_out = L_sub(st, L_var3,L_produit);
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_msu0
*
*	Purpose :
*
*		Multiply var1 by var2. Subtract the 32 bit result to L_var3.
*		Control saturation and set overflow_flag.
*		Same as L_msu without left shift.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var3
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var3 <= 0x7fff ffff.
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_msu0(tetra_codec* st, Word32 L_var3, Word16 var1, Word16 var2)
  {
   Word32 L_var_out;
   Word32 L_produit;

   L_produit = L_mult0(var1,var2);
   L_var_out = L_sub(st, L_var3,L_produit);
   return(L_var_out);
  }


/************************************************************************
*
*	Function Name : L_mult
*
*	Purpose :
*
*		L_mult is the 32 bit result of the multiplication of var1 times var2
*		with one shift left i.e.:
*			L_mult(st, var1,var2) = shl(st, (var1 times var2),1) and
*			L_mult(st, -32768,-32768) = 2147483647.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_mult(tetra_codec* st, Word16 var1,Word16 var2)
  {
   Word32 L_var_out;

   L_var_out = (Word32)var1 * (Word32)var2;
   if (L_var_out != (Word32)0x40000000)
     {
      L_var_out *= 2;
     }
   else
     {
      st->Overflow = 1;
      L_var_out = MAX_32;
     }

   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_mult0
*
*	Purpose :
*
*		L_mult is the 32 bit result of the multiplication of var1 times var2.
*		Same as L_mult without left shift.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_mult0(Word16 var1, Word16 var2)
{
  Word32 L_var_out;

  L_var_out = (Word32)var1 * (Word32)var2;

  return(L_var_out);
}


/************************************************************************
*
*	Function Name : L_negate
*
*	Purpose :
*
*		Negate the 32 bit variable L_var1 with saturation; saturate in the case
*		where input is -2147483648 (0x8000 0000).
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_negate(Word32 L_var1)
  {
   Word32 L_var_out;

   L_var_out = (L_var1 == MIN_32) ? MAX_32 : -L_var1;
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_shl
*
*	Purpose :
*
*		Arithmetically shift the 32 bit input L_var1 left var2 positions. Zero fill the var2 LSB
*		of the result. If var2 is negative, L_var1 right by -var2 arithmetically shift with sign
*		extension. Saturate the result in case of underflows or overflows.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/
Word32 L_shl(tetra_codec* st, Word32 L_var1, Word16 var2)
  {
   Word32 L_var_out = 0;
if (var2 <= 0)
     {
      L_var_out = L_shr(st, L_var1, (Word16)(-var2) );
     }
   else
     {
      for(;var2>0;var2--)
        {
         if (L_var1 > (Word32) 0X3fffffff)
           {
            st->Overflow = 1;
            L_var_out = MAX_32;
            break;
           }
         else
           {
            if (L_var1 < (Word32) 0xc0000000)
              {
               st->Overflow = 1;
               L_var_out = MIN_32;
               break;
              }
           }
         L_var1 *= 2;
         L_var_out = L_var1;
        }
     }
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_shr
*
*	Purpose :
*
*		Arithmetically shift the 32 bit input L_var1 right var2 positions with
*		sign extension. If var2 is negative, arithmetically shift L_var1 left
*		by -var2 and zero fill the var2 LSB of the result. Saturate the result
*		in case of underflows or overflows.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_shr(tetra_codec* st, Word32 L_var1, Word16 var2)
  {
   Word32 L_var_out;

   if (var2 < 0)
     {
      L_var_out = L_shl(st, L_var1,(Word16)(-var2) );
     }
   else
     {
      if (var2 >= 31)
 {
  L_var_out = (L_var1 < 0L) ? -1 : 0;
 }
      else
 {
  if (L_var1<0)
    {
     L_var_out = ~((~L_var1) >> var2);
    }
  else
    {
     L_var_out = L_var1 >> var2;
    }
 }
     }
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_shr_r
*
*	Purpose :
*
*		Same as L_shr(st, L_var1,var2)but with rounding. Saturate the result in case|
*		of underflows or overflows :
*		If var2 is greater than zero :
*			L_shr_r(st, var1,var2) = L_shr(st, L_add(st, L_var1,2**(var2-1)),var2)
*		If var2 is less than zero :
*			L_shr_r(st, var1,var2) = L_shr(st, L_var1,var2).
*
*	Complexity Weight : 3
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_shr_r(tetra_codec* st, Word32 L_var1,Word16 var2)
  {
   Word32 L_var_out;

   if (var2 > 31)
     {
      L_var_out = 0;
     }
   else
     {
      L_var_out = L_shr(st, L_var1,var2);
      if (var2 > 0)
 {
  if ( (L_var1 & ( (Word32)1 << (var2-1) )) != 0)
    {
     L_var_out++;
    }
 }
     }
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : L_sub
*
*	Purpose :
*
*		32 bits subtraction of the two 32 bits variables (L_var1-L_var2) with
*		overflow control and saturation; the result is set at +214783647 when
*		overflow occurs or at -214783648 when underflow occurs.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*		L_var2
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var2 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		L_var_out
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var_out <= 0x7fff ffff.
*
************************************************************************/

Word32 L_sub(tetra_codec* st, Word32 L_var1, Word32 L_var2)
  {
   Word32 L_var_out;

   L_var_out = L_var1 - L_var2;

   if (((L_var1 ^ L_var2) & MIN_32) != 0)
     {
      if ((L_var_out ^ L_var1) & MIN_32)
        {
  L_var_out = (L_var1 < 0L) ? MIN_32 : MAX_32;
         st->Overflow = 1;
        }
     }
   return(L_var_out);
  }

/************************************************************************
*
*	Function Name : mult
*
*	Purpose :
*
*		Performs the multiplication of var1 by var2 and gives a 16 bit result
*		which is scaled i.e.:
*			mult(st, var1,var2) = shr(st, (var1 times var2),15) and
*			mult(st, -32768,-32768) = 32767.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 mult(tetra_codec* st, Word16 var1, Word16 var2)
  {
   Word16 var_out;
   Word32 L_produit;

   L_produit = (Word32)var1 * (Word32)var2;

   L_produit = (L_produit & (Word32) 0xffff8000) >> 15;

   if (L_produit & (Word32) 0x00010000)
     L_produit = L_produit | (Word32) 0xffff0000;

   var_out = sature(st, L_produit);
   return(var_out);
  }

/************************************************************************
*
*	Function Name : mult_r
*
*	Purpose :
*
*		Same as mult with rounding, i.e.:
*			mult_r(st, var1,var2) = shr(st, ((var1*var2) + 16384),15) and
*			mult_r(st, -32768,-32768) = 32767.
*
*	Complexity Weight : 2
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 mult_r(tetra_codec* st, Word16 var1, Word16 var2)
  {
   Word16 var_out;
   Word32 L_produit_arr;

   L_produit_arr = (Word32)var1 * (Word32)var2; /* product */
   L_produit_arr += (Word32) 0x00004000;        /* round */
   L_produit_arr &= (Word32) 0xffff8000;
   L_produit_arr >>= 15;                        /* shift */

   if (L_produit_arr & (Word32) 0x00010000)   /*sign extend when necessary*/
     {
      L_produit_arr |= (Word32) 0xffff0000;
     }

   var_out = sature(st, L_produit_arr);
   return(var_out);
  }

/************************************************************************
*
*	Function Name : negate
*
*	Purpose :
*
*		Negate var1 with saturation, saturate in the case where input is -32768:
*			negate(var1) = sub(st, 0,var1).
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 negate(Word16 var1)
  {
   Word16 var_out;

   var_out = (var1 == MIN_16) ? MAX_16 : -var1;
   return(var_out);
  }

/************************************************************************
*
*	Function Name : norm_l
*
*	Purpose :
*
*		Produces the number of left shift needed to normalize the 32 bit variable
*		L_var1 for positive values on the interval with minimum of
*		1073741824 and maximum of 2147483647, and for negative values on the
*		interval with minimum of -2147483648 and maximum of -1073741824; in order
*		to normalize the result, the following operation must be done :
*			norm_L_var1 = L_shl(st, L_var1,norm_l(L_var1)).
*
*	Complexity Weight : 30
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= var_out <= 0x0000 001f.
*
************************************************************************/

Word16 norm_l(Word32 L_var1)
  {
   Word16 var_out;

   if (L_var1 == 0)
     {
      var_out = 0;
     }
   else
     {
      if (L_var1 == (Word32)0xffffffff)
 {
  var_out = 31;
 }
      else
 {
  if (L_var1 < 0)
    {
     L_var1 = ~L_var1;
    }

  for(var_out = 0;L_var1 < (Word32)0x40000000;var_out++)
    {
     L_var1 <<= 1;
    }
 }
     }

   return(var_out);
  }

/************************************************************************
*
*	Function Name : norm_s
*
*	Purpose :
*
*		Produces the number of left shift needed to normalize the 16 bit variable
*		var1 for positive values on the interval with minimum of 16384 and
*		maximum of 32767, and for negative values on the interval with minimum
*		of -32768 and maximum of -16384; in order to normalize the result, the
*		following operation must be done :
*			norm_var1 = shl(st, var1,norm_s(var1)).
*
*	Complexity Weight : 15
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0x0000 0000 <= var_out <= 0x0000 000f.
*
************************************************************************/

Word16 norm_s(Word16 var1)
  {
   Word16 var_out;

   if (var1 == 0)
     {
      var_out = 0;
     }
   else
     {
      if (var1 == (Word16) 0xffff)
 {
  var_out = 15;
 }
      else
 {
  if (var1 < 0)
    {
     var1 = ~var1;
    }

  for(var_out = 0; var1 < 0x4000; var_out++)
    {
     var1 <<= 1;
    }
 }
     }

   return(var_out);
  }

/************************************************************************
*
*	Function Name : round
*
*	Purpose :
*
*		Round the lower 16 bits of the 32 bit input number into its MS 16 bits
*		with saturation. Shift the resulting bits right by 16 and return the 16
*		bit number:
*			round(L_var1) = extract_h(L_add(st, L_var1,32768))
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32 ) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 tsc_round(tetra_codec* st, Word32 L_var1)
  {
   Word16 var_out;
   Word32 L_arrondi;

   L_arrondi = L_add(st, L_var1, (Word32)0x00008000);
   var_out = extract_h(L_arrondi);
   return(var_out);
  }

/************************************************************************
*
*	Function Name : sature
*
*	Purpose :
*
*		Limit the 32 bit input to the range of a 16 bit word.
*
*	Inputs :
*
*		L_var1
*			32 bit long signed integer (Word32) whose value falls in the
*			range : 0x8000 0000 <= L_var1 <= 0x7fff ffff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 sature(tetra_codec* st, Word32 L_var1)
  {
   Word16 var_out;

   if (L_var1 > 0X00007fffL)
     {
      st->Overflow = 1;
      var_out = MAX_16;
     }
   else if (L_var1 < (Word32)0xffff8000L)
     {
      st->Overflow = 1;
      var_out = MIN_16;
     }
   else
     {
      st->Overflow = 0;
      var_out = extract_l(L_var1);
     }

   return(var_out);
  }


/************************************************************************
*
*	Function Name : shl
*
*	Purpose :
*
*		Arithmetically shift the 16 bit input var1 left var2 positions.Zero fill|
*		the var2 LSB of the result. If var2 is negative, arithmetically shift
*		var1 right by -var2 with sign extension. Saturate the result in case of
*		underflows or overflows.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 shl(tetra_codec* st, Word16 var1,Word16 var2)
  {
   Word16 var_out;
   Word32 resultat;

   if (var2 < 0)
     {
      var_out = shr(st,  var1,(Word16)(-var2) );
     }
   else
     {
      resultat = (Word32) var1 * ((Word32) 1 << var2);
     if ((var2 > 15 && var1 != 0) || (resultat != (Word32)((Word16) resultat)))
        {
         st->Overflow = 1;
         var_out = (var1 > 0) ? MAX_16 : MIN_16;
        }
      else
        {
         var_out = extract_l(resultat);
        }
     }

   return(var_out);
  }

/************************************************************************
*
*	Function Name : shr
*
*	Purpose :
*
*		Arithmetically shift the 16 bit input var1 right var2 positions with
*		sign extension. If var2 is negative, arithmetically shift var1 left by
*		-var2 with sign extension. Saturate the result in case of underflows or
*		overflows.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 shr(tetra_codec* st, Word16 var1,Word16 var2)
  {
   Word16 var_out;

   if (var2 < 0)
     {
      var_out = shl(st,  var1,(Word16)(-var2) );
     }
   else
     {
      if (var2 >= 15)
        {
         var_out = (var1 < 0) ? -1 : 0;
        }
      else
        {
         if (var1 < 0)
           {
     var_out = ~(( ~var1) >> var2 );
           }
         else
           {
            var_out = var1 >> var2;
           }
        }
     }

   return(var_out);
  }

/************************************************************************
*
*	Function Name : sub
*
*	Purpose :
*
*			Performs the subtraction (var1+var2) with overflow control and satu-
*			ration; the 16 bit result is set at +32767 when overflow occurs or at
*			-32768 when underflow occurs.
*
*	Complexity Weight : 1
*
*	Inputs :
*
*		var1
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var1 <= 0x0000 7fff.
*
*		var2
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var2 <= 0x0000 7fff.
*
*	Outputs :
*
*		none
*
*	Returned Value :
*
*		var_out
*			16 bit short signed integer (Word16) whose value falls in the
*			range : 0xffff 8000 <= var_out <= 0x0000 7fff.
*
************************************************************************/

Word16 sub(tetra_codec* st, Word16 var1,Word16 var2)
  {
   Word16 var_out;
   Word32 L_diff;

   L_diff = (Word32) var1 - var2;
   var_out = sature(st, L_diff);
   return(var_out);
  }
