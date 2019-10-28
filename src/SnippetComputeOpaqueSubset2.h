#pragma once

#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"

static INLINED int ComputeOpaqueSubsetError2(const Area& area, __m128i mc, const __m128i mwater) noexcept
{
	__m128i merrorBlock = _mm_setzero_si128();

	const __m128i mhalf = _mm_set1_epi16(32);
	const __m128i msign = _mm_set1_epi16(-0x8000);
	const __m128i mweights = gWeightsGRB;
	const __m128i mfix = gFixWeightsGRB;

	mc = _mm_packus_epi16(mc, mc);

	__m128i mtx = gTableInterpolate2_U8[0];
	__m128i mty = gTableInterpolate2_U8[1];

	mtx = _mm_maddubs_epi16(mc, mtx);
	mty = _mm_maddubs_epi16(mc, mty);

	mtx = _mm_add_epi16(mtx, mhalf);
	mty = _mm_add_epi16(mty, mhalf);

	mtx = _mm_srli_epi16(mtx, 6);
	mty = _mm_srli_epi16(mty, 6);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_unpacklo_epi64(mpacked, mpacked);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m128i mx = _mm_sub_epi16(mpixel, mtx);
		__m128i my = _mm_sub_epi16(mpixel, mty);

		mx = _mm_mullo_epi16(mx, mx);
		my = _mm_mullo_epi16(my, my);

		mx = _mm_xor_si128(mx, msign);
		my = _mm_xor_si128(my, msign);

		mx = _mm_madd_epi16(mx, mweights);
		my = _mm_madd_epi16(my, mweights);

		mx = _mm_add_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));
		my = _mm_add_epi32(my, _mm_shuffle_epi32(my, _MM_SHUFFLE(2, 3, 0, 1)));

		mx = _mm_min_epi32(mx, my);
		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, mx);

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}

	return _mm_cvtsi128_si32(merrorBlock);
}

template<int shuffle>
static INLINED int ComputeOpaqueSubsetError2Pair(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
{
	__m128i merrorBlock = _mm_setzero_si128();

	const __m128i mhalf = _mm_set1_epi16(32);
	const __m128i msign = _mm_set1_epi16(-0x8000);

	mc = _mm_shuffle_epi32(mc, shuffle);
	mc = _mm_packus_epi16(mc, mc);

	__m128i mt = gTableInterpolate2GR_U8[0];

	mt = _mm_maddubs_epi16(mc, mt);

	mt = _mm_add_epi16(mt, mhalf);

	mt = _mm_srli_epi16(mt, 6);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_shufflelo_epi16(mpacked, shuffle);
		mpixel = _mm_unpacklo_epi64(mpixel, mpixel);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m128i mx = _mm_sub_epi16(mpixel, mt);

		mx = _mm_mullo_epi16(mx, mx);

		mx = _mm_xor_si128(mx, msign);

		mx = _mm_madd_epi16(mx, mweights);

		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));
		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, mx);

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}

	return _mm_cvtsi128_si32(merrorBlock);
}

static INLINED int ComputeOpaqueSubsetTable2(const Area& area, __m128i mc, uint64_t& indices) noexcept
{
	const __m128i mhalf = _mm_set1_epi16(32);

	mc = _mm_packus_epi16(mc, mc);

	__m128i mtx = gTableInterpolate2_U8[0];
	__m128i mty = gTableInterpolate2_U8[1];

	mtx = _mm_maddubs_epi16(mc, mtx);
	mty = _mm_maddubs_epi16(mc, mty);

	mtx = _mm_add_epi16(mtx, mhalf);
	mty = _mm_add_epi16(mty, mhalf);

	mtx = _mm_srli_epi16(mtx, 6);
	mty = _mm_srli_epi16(mty, 6);

	const __m128i mopaque = _mm_set_epi16(0, 0, 0, 255, 0, 0, 0, 255);
	mtx = _mm_or_si128(mtx, mopaque);
	mty = _mm_or_si128(mty, mopaque);

	Modulations state;
	_mm_store_si128((__m128i*)&state.Values_I16[0], mtx);
	_mm_store_si128((__m128i*)&state.Values_I16[2], mty);

	int error = ComputeSubsetTable(area, gWeightsAGRB, gFixWeightsAGRB, state, 4);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		uint64_t index = static_cast<uint32_t>(state.Best[i]);
		indices |= index << (area.Indices[i] << 1);
	}

	return error;
}
