
#include "tetra-codec.h"

#include <stdlib.h>

#include "tetra-codec-impl.h"

//------------------------------------------------------------------------------

#define CODED_FRAME_BIT_COUNT 137

//------------------------------------------------------------------------------

tetra_codec* tetra_encoder_create(void)
{
	tetra_codec* st = malloc(sizeof(tetra_codec));
	if (st) {
		init_tetra_codec(st);
		Init_Pre_Process(st);
		Init_Coder_Tetra(st);
	}
	return st;
}

//------------------------------------------------------------------------------

tetra_codec* tetra_decoder_create(void)
{
	tetra_codec* st = malloc(sizeof(tetra_codec));
	if (st) {
		init_tetra_codec(st);
		Init_Decod_Tetra(st);
	}
	return st;
}

//------------------------------------------------------------------------------

void tetra_codec_destroy(tetra_codec* st)
{
	free(st);
}

//------------------------------------------------------------------------------

static void pack_coded_bits(const int16_t* serial, uint8_t* coded)
{
	int bit = 8;
	for (int i = 0; i != CODED_FRAME_BIT_COUNT; ++i) {
		*coded = (uint8_t)((*coded << 1) | serial[i]);
		if (--bit == 0) {
			++coded;
			bit = 8;
		}
	}
	*coded <<= 7;

}

void tetra_encode(tetra_codec* st, const int16_t* pcm, uint8_t* coded)
{
	for (int i = 0; i != TSC_L_frame; ++i) {
		st->new_speech[i] = pcm[i];
	}

	Pre_Process(st, st->new_speech, TSC_L_frame);

	int16_t ana[23];
	int16_t syn[TSC_L_frame];
	Coder_Tetra(st, ana, syn);

	int16_t serial[1 + CODED_FRAME_BIT_COUNT];
	Prm2bits_Tetra(st, ana, serial);

	pack_coded_bits(serial + 1, coded);
}

//------------------------------------------------------------------------------

static void unpack_coded_bits(const uint8_t* coded, int16_t* serial)
{
	int bit = 7;
	for (int i = 0; i != CODED_FRAME_BIT_COUNT; ++i) {
		serial[i] = (*coded >> bit) & 1;
		if (--bit < 0) {
			++coded;
			bit = 7;
		}
	}
}

void tetra_decode(tetra_codec* st, const uint8_t* coded, int16_t* pcm, int bfi)
{
	int16_t serial[1 + CODED_FRAME_BIT_COUNT];
	serial[0] = bfi ? 1 : 0;
	unpack_coded_bits(coded, serial + 1);

	int16_t parm[24];
	Bits2prm_Tetra(st, serial, parm);

	Decod_Tetra(st, parm, pcm);

	Post_Process(st, pcm, TSC_L_frame);
}

//------------------------------------------------------------------------------

