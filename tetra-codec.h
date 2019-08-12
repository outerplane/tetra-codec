#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _MSC_VER
#define TETRA_CODEC_EXPORT __declspec(dllexport)
#else
#define TETRA_CODEC_EXPORT
#endif

typedef struct tetra_codec tetra_codec;

TETRA_CODEC_EXPORT tetra_codec* tetra_encoder_create(void);
TETRA_CODEC_EXPORT tetra_codec* tetra_decoder_create(void);
TETRA_CODEC_EXPORT void tetra_codec_destroy(tetra_codec* st);

TETRA_CODEC_EXPORT void tetra_encode(tetra_codec* st, const int16_t* pcm, uint8_t* coded);
TETRA_CODEC_EXPORT void tetra_decode(tetra_codec* st, const uint8_t* coded, int16_t* pcm, int bfi);

#ifdef __cplusplus
}
#endif
