#pragma once

#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"

#include "SnippetHorizontalSum4.h"

static INLINED int ComputeOpaqueSubsetError3(const Area& area, __m128i mc, const __m128i mwater) noexcept
{
	__m128i merrorBlock = _mm_setzero_si128();

#if defined(OPTION_AVX2)
	const __m256i vweights = _mm256_broadcastsi128_si256(gWeightsGRB);
	const __m128i mfix = gFixWeightsGRB;

	const __m256i vhalf = _mm256_set1_epi16(32);
	const __m256i vsign = _mm256_set1_epi16(-0x8000);

	mc = _mm_packus_epi16(mc, mc);
	__m256i vc = _mm256_broadcastsi128_si256(mc);

	__m256i vtx = *(const __m256i*)&gTableInterpolate3_U8[0];
	__m256i vty = *(const __m256i*)&gTableInterpolate3_U8[2];

	vtx = _mm256_maddubs_epi16(vc, vtx);
	vty = _mm256_maddubs_epi16(vc, vty);

	vtx = _mm256_add_epi16(vtx, vhalf);
	vty = _mm256_add_epi16(vty, vhalf);

	vtx = _mm256_srli_epi16(vtx, 6);
	vty = _mm256_srli_epi16(vty, 6);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m256i vpixel = _mm256_broadcastq_epi64(mpacked);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m256i vx = _mm256_sub_epi16(vpixel, vtx);
		__m256i vy = _mm256_sub_epi16(vpixel, vty);

		vx = _mm256_mullo_epi16(vx, vx);
		vy = _mm256_mullo_epi16(vy, vy);

		vx = _mm256_xor_si256(vx, vsign);
		vy = _mm256_xor_si256(vy, vsign);

		vx = _mm256_madd_epi16(vx, vweights);
		vy = _mm256_madd_epi16(vy, vweights);

		vx = _mm256_add_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));
		vy = _mm256_add_epi32(vy, _mm256_shuffle_epi32(vy, _MM_SHUFFLE(2, 3, 0, 1)));

		vx = _mm256_min_epi32(vx, vy);
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, _mm_min_epi32(_mm256_extracti128_si256(vx, 1), _mm256_castsi256_si128(vx)));

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}
#else
	const __m128i mhalf = _mm_set1_epi16(32);
	const __m128i msign = _mm_set1_epi16(-0x8000);
	const __m128i mweights = gWeightsGRB;
	const __m128i mfix = gFixWeightsGRB;

	mc = _mm_packus_epi16(mc, mc);

	__m128i mtx = gTableInterpolate3_U8[0];
	__m128i mty = gTableInterpolate3_U8[1];
	__m128i mtz = gTableInterpolate3_U8[2];
	__m128i mtw = gTableInterpolate3_U8[3];

	mtx = _mm_maddubs_epi16(mc, mtx);
	mty = _mm_maddubs_epi16(mc, mty);
	mtz = _mm_maddubs_epi16(mc, mtz);
	mtw = _mm_maddubs_epi16(mc, mtw);

	mtx = _mm_add_epi16(mtx, mhalf);
	mty = _mm_add_epi16(mty, mhalf);
	mtz = _mm_add_epi16(mtz, mhalf);
	mtw = _mm_add_epi16(mtw, mhalf);

	mtx = _mm_srli_epi16(mtx, 6);
	mty = _mm_srli_epi16(mty, 6);
	mtz = _mm_srli_epi16(mtz, 6);
	mtw = _mm_srli_epi16(mtw, 6);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_unpacklo_epi64(mpacked, mpacked);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m128i mx = _mm_sub_epi16(mpixel, mtx);
		__m128i my = _mm_sub_epi16(mpixel, mty);
		__m128i mz = _mm_sub_epi16(mpixel, mtz);
		__m128i mw = _mm_sub_epi16(mpixel, mtw);

		mx = _mm_mullo_epi16(mx, mx);
		my = _mm_mullo_epi16(my, my);
		mz = _mm_mullo_epi16(mz, mz);
		mw = _mm_mullo_epi16(mw, mw);

		mx = _mm_xor_si128(mx, msign);
		my = _mm_xor_si128(my, msign);
		mz = _mm_xor_si128(mz, msign);
		mw = _mm_xor_si128(mw, msign);

		mx = _mm_madd_epi16(mx, mweights);
		my = _mm_madd_epi16(my, mweights);
		mz = _mm_madd_epi16(mz, mweights);
		mw = _mm_madd_epi16(mw, mweights);

		mx = _mm_add_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));
		my = _mm_add_epi32(my, _mm_shuffle_epi32(my, _MM_SHUFFLE(2, 3, 0, 1)));
		mz = _mm_add_epi32(mz, _mm_shuffle_epi32(mz, _MM_SHUFFLE(2, 3, 0, 1)));
		mw = _mm_add_epi32(mw, _mm_shuffle_epi32(mw, _MM_SHUFFLE(2, 3, 0, 1)));

		mx = _mm_min_epi32(mx, my);
		mz = _mm_min_epi32(mz, mw);
		mx = _mm_min_epi32(mx, mz);
		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, mx);

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}
#endif

	return _mm_cvtsi128_si32(merrorBlock);
}

template<int shuffle>
static INLINED int ComputeOpaqueSubsetError3Pair(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
{
	__m128i merrorBlock = _mm_setzero_si128();

#if defined(OPTION_AVX2)
	const __m256i vweights = _mm256_broadcastsi128_si256(mweights);

	const __m256i vhalf = _mm256_set1_epi16(32);
	const __m256i vsign = _mm256_set1_epi16(-0x8000);
	const __m128i mfix2 = _mm_add_epi32(mfix, mfix);

	mc = _mm_shuffle_epi32(mc, shuffle);
	mc = _mm_packus_epi16(mc, mc);
	__m256i vc = _mm256_broadcastsi128_si256(mc);

	__m256i vt = *(const __m256i*)gTableInterpolate3GR_U8;

	vt = _mm256_maddubs_epi16(vc, vt);

	vt = _mm256_add_epi16(vt, vhalf);

	vt = _mm256_srli_epi16(vt, 6);

	__m256i vtx = _mm256_permute4x64_epi64(vt, 0x44);
	__m256i vty = _mm256_permute4x64_epi64(vt, 0xEE);

	int k = static_cast<int>(area.Count);
	const __m256i* p = (const __m256i*)area.DataMask_I16;

	while ((k -= 2) >= 0)
	{
		__m256i vpacked = _mm256_load_si256(p);
		__m256i vpixel = _mm256_shufflelo_epi16(vpacked, shuffle);
		vpixel = _mm256_unpacklo_epi64(vpixel, vpixel);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix2);

		__m256i vx = _mm256_sub_epi16(vpixel, vtx);
		__m256i vy = _mm256_sub_epi16(vpixel, vty);

		vx = _mm256_mullo_epi16(vx, vx);
		vy = _mm256_mullo_epi16(vy, vy);

		vx = _mm256_xor_si256(vx, vsign);
		vy = _mm256_xor_si256(vy, vsign);

		vx = _mm256_madd_epi16(vx, vweights);
		vy = _mm256_madd_epi16(vy, vweights);

		vx = _mm256_min_epi32(vx, vy);
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_castsi256_si128(vx));
		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_extracti128_si256(vx, 1));

		p++;

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}

	if (k & 1)
	{
		__m128i mpacked = _mm_load_si128((const __m128i*)p);
		__m128i mpixel = _mm_shufflelo_epi16(mpacked, shuffle);
		__m256i vpixel = _mm256_broadcastq_epi64(mpixel);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m256i vx = _mm256_sub_epi16(vpixel, vt);

		vx = _mm256_mullo_epi16(vx, vx);

		vx = _mm256_xor_si256(vx, vsign);

		vx = _mm256_madd_epi16(vx, vweights);

		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, _mm_min_epi32(_mm256_extracti128_si256(vx, 1), _mm256_castsi256_si128(vx)));
	}
#else
	const __m128i mhalf = _mm_set1_epi16(32);
	const __m128i msign = _mm_set1_epi16(-0x8000);

	mc = _mm_shuffle_epi32(mc, shuffle);
	mc = _mm_packus_epi16(mc, mc);

	__m128i mtx = gTableInterpolate3GR_U8[0];
	__m128i mty = gTableInterpolate3GR_U8[1];

	mtx = _mm_maddubs_epi16(mc, mtx);
	mty = _mm_maddubs_epi16(mc, mty);

	mtx = _mm_add_epi16(mtx, mhalf);
	mty = _mm_add_epi16(mty, mhalf);

	mtx = _mm_srli_epi16(mtx, 6);
	mty = _mm_srli_epi16(mty, 6);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_shufflelo_epi16(mpacked, shuffle);
		mpixel = _mm_unpacklo_epi64(mpixel, mpixel);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m128i mx = _mm_sub_epi16(mpixel, mtx);
		__m128i my = _mm_sub_epi16(mpixel, mty);

		mx = _mm_mullo_epi16(mx, mx);
		my = _mm_mullo_epi16(my, my);

		mx = _mm_xor_si128(mx, msign);
		my = _mm_xor_si128(my, msign);

		mx = _mm_madd_epi16(mx, mweights);
		my = _mm_madd_epi16(my, mweights);

		mx = _mm_min_epi32(mx, my);
		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));
		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, mx);

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}
#endif

	return _mm_cvtsi128_si32(merrorBlock);
}

static INLINED int ComputeOpaqueSubsetTable3(const Area& area, __m128i mc, uint64_t& indices) noexcept
{
	const __m128i mhalf = _mm_set1_epi16(32);

	mc = _mm_packus_epi16(mc, mc);

	__m128i mtx = gTableInterpolate3_U8[0];
	__m128i mty = gTableInterpolate3_U8[1];
	__m128i mtz = gTableInterpolate3_U8[2];
	__m128i mtw = gTableInterpolate3_U8[3];

	mtx = _mm_maddubs_epi16(mc, mtx);
	mty = _mm_maddubs_epi16(mc, mty);
	mtz = _mm_maddubs_epi16(mc, mtz);
	mtw = _mm_maddubs_epi16(mc, mtw);

	mtx = _mm_add_epi16(mtx, mhalf);
	mty = _mm_add_epi16(mty, mhalf);
	mtz = _mm_add_epi16(mtz, mhalf);
	mtw = _mm_add_epi16(mtw, mhalf);

	mtx = _mm_srli_epi16(mtx, 6);
	mty = _mm_srli_epi16(mty, 6);
	mtz = _mm_srli_epi16(mtz, 6);
	mtw = _mm_srli_epi16(mtw, 6);

	const __m128i mopaque = _mm_set_epi16(0, 0, 0, 255, 0, 0, 0, 255);
	mtx = _mm_or_si128(mtx, mopaque);
	mty = _mm_or_si128(mty, mopaque);
	mtz = _mm_or_si128(mtz, mopaque);
	mtw = _mm_or_si128(mtw, mopaque);

	Modulations state;
	_mm_store_si128((__m128i*)&state.Values_I16[0], mtx);
	_mm_store_si128((__m128i*)&state.Values_I16[2], mty);
	_mm_store_si128((__m128i*)&state.Values_I16[4], mtz);
	_mm_store_si128((__m128i*)&state.Values_I16[6], mtw);

	int error = ComputeSubsetTable(area, gWeightsAGRB, gFixWeightsAGRB, state, 8);

	for (size_t i = 0, n = area.Count; i < n; i++)
	{
		uint64_t index = static_cast<uint32_t>(state.Best[i]);
		indices |= index << (area.Indices[i] * 3);
	}

	return error;
}
