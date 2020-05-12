#pragma once

#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"

static INLINED int ComputeOpaqueSubsetError2(const Area& area, __m128i mc, const __m128i mwater) noexcept
{
	__m128i merrorBlock = _mm_setzero_si128();

	const __m128i mweights = gWeightsGRB;
	const __m128i mfix = gFixWeightsGRB;

#if defined(OPTION_AVX2)
	const __m256i vweights = _mm256_broadcastq_epi64(mweights);

	const __m256i vhalf = _mm256_set1_epi16(32);
	const __m256i vsign = _mm256_set1_epi16(-0x8000);
	const __m128i mfix4 = _mm_slli_epi32(mfix, 2);

	mc = _mm_packus_epi16(mc, mc);
	__m256i vc = _mm256_broadcastq_epi64(mc);

	__m256i vt = *(const __m256i*)gTableInterpolate2_U8;

	vt = _mm256_maddubs_epi16(vc, vt);

	vt = _mm256_add_epi16(vt, vhalf);

	vt = _mm256_srli_epi16(vt, 6);

	__m256i vtx = _mm256_permute4x64_epi64(vt, 0x44);
	__m256i vty = _mm256_permute4x64_epi64(vt, 0xEE);

	int k = static_cast<int>(area.Count);
	const __m256i* p = (const __m256i*)area.DataMask_I16;

	while ((k -= 4) >= 0)
	{
		__m256i vpacked0 = _mm256_load_si256(&p[0]);
		__m256i vpacked1 = _mm256_load_si256(&p[1]);
		__m256i vpixel0 = _mm256_unpacklo_epi64(vpacked0, vpacked0);
		__m256i vpixel1 = _mm256_unpacklo_epi64(vpacked1, vpacked1);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix4);

		__m256i vx = _mm256_sub_epi16(vpixel0, vtx);
		__m256i vy = _mm256_sub_epi16(vpixel0, vty);
		__m256i vz = _mm256_sub_epi16(vpixel1, vtx);
		__m256i vw = _mm256_sub_epi16(vpixel1, vty);

		vx = _mm256_mullo_epi16(vx, vx);
		vy = _mm256_mullo_epi16(vy, vy);
		vz = _mm256_mullo_epi16(vz, vz);
		vw = _mm256_mullo_epi16(vw, vw);

		vx = _mm256_xor_si256(vx, vsign);
		vy = _mm256_xor_si256(vy, vsign);
		vz = _mm256_xor_si256(vz, vsign);
		vw = _mm256_xor_si256(vw, vsign);

		vx = _mm256_madd_epi16(vx, vweights);
		vy = _mm256_madd_epi16(vy, vweights);
		vz = _mm256_madd_epi16(vz, vweights);
		vw = _mm256_madd_epi16(vw, vweights);

		vx = _mm256_add_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));
		vy = _mm256_add_epi32(vy, _mm256_shuffle_epi32(vy, _MM_SHUFFLE(2, 3, 0, 1)));
		vz = _mm256_add_epi32(vz, _mm256_shuffle_epi32(vz, _MM_SHUFFLE(2, 3, 0, 1)));
		vw = _mm256_add_epi32(vw, _mm256_shuffle_epi32(vw, _MM_SHUFFLE(2, 3, 0, 1)));

		vx = _mm256_min_epi32(vx, vy);
		vz = _mm256_min_epi32(vz, vw);
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));
		vz = _mm256_min_epi32(vz, _mm256_shuffle_epi32(vz, _MM_SHUFFLE(1, 0, 3, 2)));

		vx = _mm256_add_epi32(vx, vz);

		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_castsi256_si128(vx));
		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_extracti128_si256(vx, 1));

		p += 2;

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}

	if (k & 2)
	{
		__m256i vpacked = _mm256_load_si256(p);
		__m256i vpixel = _mm256_unpacklo_epi64(vpacked, vpacked);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);
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

		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_castsi256_si128(vx));
		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_extracti128_si256(vx, 1));

		p++;
	}

	if (k & 1)
	{
		__m128i mpacked = _mm_load_si128((const __m128i*)p);
		__m256i vpixel = _mm256_broadcastq_epi64(mpacked);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m256i vx = _mm256_sub_epi16(vpixel, vt);

		vx = _mm256_mullo_epi16(vx, vx);

		vx = _mm256_xor_si256(vx, vsign);

		vx = _mm256_madd_epi16(vx, vweights);

		vx = _mm256_add_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));

		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, _mm_min_epi32(_mm256_extracti128_si256(vx, 1), _mm256_castsi256_si128(vx)));
	}
#else
	const __m128i mhalf = _mm_set1_epi16(32);
	const __m128i msign = _mm_set1_epi16(-0x8000);

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
#endif

	return _mm_cvtsi128_si32(merrorBlock);
}

template<int shuffle>
static INLINED int ComputeOpaqueSubsetError2Pair(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
{
	__m128i merrorBlock = _mm_setzero_si128();

#if defined(OPTION_AVX2)
	const __m256i vweights = _mm256_broadcastq_epi64(mweights);

	const __m256i vhalf = _mm256_set1_epi16(32);
	const __m256i vsign = _mm256_set1_epi16(-0x8000);
	const __m128i mfix4 = _mm_slli_epi32(mfix, 2);

	mc = _mm_shuffle_epi32(mc, shuffle);
	mc = _mm_packus_epi16(mc, mc);
	__m256i vc = _mm256_broadcastq_epi64(mc);

	__m256i vt = *(const __m256i*)gTableInterpolate2GR_U8;

	vt = _mm256_maddubs_epi16(vc, vt);

	vt = _mm256_add_epi16(vt, vhalf);

	vt = _mm256_srli_epi16(vt, 6);

	int k = static_cast<int>(area.Count);
	const __m256i* p = (const __m256i*)area.DataMask_I16;

	while ((k -= 4) >= 0)
	{
		__m256i vpacked0 = _mm256_load_si256(&p[0]);
		__m256i vpacked1 = _mm256_load_si256(&p[1]);
		__m256i vpixel0 = _mm256_shufflelo_epi16(vpacked0, shuffle);
		__m256i vpixel1 = _mm256_shufflelo_epi16(vpacked1, shuffle);
		vpixel0 = _mm256_unpacklo_epi64(vpixel0, vpixel0);
		vpixel1 = _mm256_unpacklo_epi64(vpixel1, vpixel1);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix4);

		__m256i vx = _mm256_sub_epi16(vpixel0, vt);
		__m256i vy = _mm256_sub_epi16(vpixel1, vt);

		vx = _mm256_mullo_epi16(vx, vx);
		vy = _mm256_mullo_epi16(vy, vy);

		vx = _mm256_xor_si256(vx, vsign);
		vy = _mm256_xor_si256(vy, vsign);

		vx = _mm256_madd_epi16(vx, vweights);
		vy = _mm256_madd_epi16(vy, vweights);

		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));
		vy = _mm256_min_epi32(vy, _mm256_shuffle_epi32(vy, _MM_SHUFFLE(2, 3, 0, 1)));
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));
		vy = _mm256_min_epi32(vy, _mm256_shuffle_epi32(vy, _MM_SHUFFLE(1, 0, 3, 2)));

		vx = _mm256_add_epi32(vx, vy);

		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_castsi256_si128(vx));
		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_extracti128_si256(vx, 1));

		p += 2;

		if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
			break;
	}

	if (k & 2)
	{
		__m256i vpacked = _mm256_load_si256(p);
		__m256i vpixel = _mm256_shufflelo_epi16(vpacked, shuffle);
		vpixel = _mm256_unpacklo_epi64(vpixel, vpixel);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);
		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m256i vx = _mm256_sub_epi16(vpixel, vt);

		vx = _mm256_mullo_epi16(vx, vx);

		vx = _mm256_xor_si256(vx, vsign);

		vx = _mm256_madd_epi16(vx, vweights);

		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(2, 3, 0, 1)));
		vx = _mm256_min_epi32(vx, _mm256_shuffle_epi32(vx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_castsi256_si128(vx));
		merrorBlock = _mm_add_epi32(merrorBlock, _mm256_extracti128_si256(vx, 1));

		p++;
	}

	if (k & 1)
	{
		__m128i mpacked = _mm_load_si128((const __m128i*)p);
		__m128i mpixel = _mm_shufflelo_epi16(mpacked, shuffle);
		mpixel = _mm_unpacklo_epi64(mpixel, mpixel);

		merrorBlock = _mm_add_epi32(merrorBlock, mfix);

		__m128i mx = _mm_sub_epi16(mpixel, _mm256_castsi256_si128(vt));

		mx = _mm_mullo_epi16(mx, mx);

		mx = _mm_xor_si128(mx, _mm256_castsi256_si128(vsign));

		mx = _mm_madd_epi16(mx, _mm256_castsi256_si128(vweights));

		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));
		mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

		merrorBlock = _mm_add_epi32(merrorBlock, mx);
	}
#else
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
#endif

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
