
#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"
#include "Bc7Pca.h"

#include "SnippetDecompressIndexedSubset.h"
#include "SnippetInsertRemoveZeroBit.h"
#include "SnippetLevelsMinimum.h"
#include "SnippetLevelsBuffer.h"

// https://docs.microsoft.com/en-us/windows/desktop/direct3d11/bc7-format-mode-reference#mode-7

namespace Mode7 {

	constexpr int LevelsCapacity = 20;

#if defined(OPTION_COUNTERS)
	static std::atomic_int gComputeSubsetError2, gComputeSubsetError2AG, gComputeSubsetError2AR, gComputeSubsetError2AGR, gComputeSubsetError2AGB;
#endif

	static ALWAYS_INLINED int Max(int x, int y) noexcept
	{
		return (x > y) ? x : y;
	}

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept
	{
		uint64_t data0 = *(const uint64_t*)&input[0];
		uint64_t data1 = *(const uint64_t*)&input[8];

		data0 >>= 8;

		size_t partitionIndex = data0 & 0x3F; data0 >>= 6;

		__m128i mc0 = _mm_setzero_si128();
		__m128i mc1 = _mm_setzero_si128();

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 4); data0 >>= 5;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 5); data0 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data0), 4); data0 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data0), 5); data0 >>= 5;

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 2); data0 >>= 5;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 3); data0 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data0), 2); data0 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data0), 3); data0 >>= 5;

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 6); data0 >>= 5;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 7);
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data1), 6); data1 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data1), 7); data1 >>= 5;

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data1), 0); data1 >>= 5;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data1), 1); data1 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data1), 0); data1 >>= 5;
		mc1 = _mm_insert_epi16(mc1, static_cast<int>(data1), 1); data1 >>= 5;

		const __m128i m5 = _mm_set1_epi16(0x1F);
		mc0 = _mm_and_si128(mc0, m5);
		mc1 = _mm_and_si128(mc1, m5);

		mc0 = _mm_or_si128(_mm_slli_epi16(mc0, 3), _mm_srli_epi16(mc0, 3));
		mc1 = _mm_or_si128(_mm_slli_epi16(mc1, 3), _mm_srli_epi16(mc1, 3));

		int pbits0 = (static_cast<int>(data1 & 1) << 2) + (static_cast<int>(data1 & 2) << 17); data1 >>= 2;
		int pbits1 = (static_cast<int>(data1 & 1) << 2) + (static_cast<int>(data1 & 2) << 17); data1 >>= 2;
		mc0 = _mm_or_si128(mc0, _mm_shuffle_epi32(_mm_cvtsi32_si128(pbits0), 0));
		mc1 = _mm_or_si128(mc1, _mm_shuffle_epi32(_mm_cvtsi32_si128(pbits1), 0));

		data1 = InsertZeroBit(data1, gTableShrinked22[partitionIndex] * 2);
		data1 = InsertZeroBit(data1, 1);

		DecompressIndexedSubset<2>(mc0, gTableSelection12[partitionIndex], (int*)output.ImageRows_U8, data1);
		DecompressIndexedSubset<2>(mc1, gTableSelection22[partitionIndex], (int*)output.ImageRows_U8, data1);

		output.BestParameter = partitionIndex;
		output.BestMode = 7;
	}

	static INLINED void ComposeBlock(uint8_t output[16], __m128i mc0, __m128i mc1, uint64_t indices, const size_t partitionIndex) noexcept
	{
		uint64_t data1 = RemoveZeroBit(RemoveZeroBit(indices, 1), gTableShrinked22[partitionIndex] * 2);

		mc0 = _mm_srli_epi16(mc0, 2);
		mc1 = _mm_srli_epi16(mc1, 2);

		data1 <<= 1; data1 |= _mm_extract_epi16(mc1, 1) & 1;
		data1 <<= 1; data1 |= _mm_extract_epi16(mc1, 0) & 1;
		data1 <<= 1; data1 |= _mm_extract_epi16(mc0, 1) & 1;
		data1 <<= 1; data1 |= _mm_extract_epi16(mc0, 0) & 1;

		mc0 = _mm_srli_epi16(mc0, 1);
		mc1 = _mm_srli_epi16(mc1, 1);

		data1 <<= 5; data1 |= _mm_extract_epi16(mc1, 1);
		data1 <<= 5; data1 |= _mm_extract_epi16(mc1, 0);
		data1 <<= 5; data1 |= _mm_extract_epi16(mc0, 1);
		data1 <<= 5; data1 |= _mm_extract_epi16(mc0, 0);

		data1 <<= 5; data1 |= _mm_extract_epi16(mc1, 7);
		data1 <<= 5; data1 |= _mm_extract_epi16(mc1, 6);
		uint64_t data0 = _mm_extract_epi16(mc0, 7);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc0, 6);

		data0 <<= 5; data0 |= _mm_extract_epi16(mc1, 3);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc1, 2);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc0, 3);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc0, 2);

		data0 <<= 5; data0 |= _mm_extract_epi16(mc1, 5);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc1, 4);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc0, 5);
		data0 <<= 5; data0 |= _mm_extract_epi16(mc0, 4);

		data0 <<= 6; data0 |= partitionIndex;

		data0 <<= 8; data0 |= 1 << 7;

		*(uint64_t*)&output[0] = data0;
		*(uint64_t*)&output[8] = data1;
	}

	static INLINED int ComputeSubsetTransparentError2(const Area& area, const int alpha) noexcept
	{
		int error = static_cast<int>(area.Count - area.Active);
		if (error)
		{
			error *= kAlpha;
			int v = gTableDeltas2_Value8[0][alpha];
			error *= v * v;
		}

		return error;
	}

	static INLINED int ComputeSubsetError2(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
	{
		__m128i merrorBlock = _mm_setzero_si128();

#if defined(OPTION_AVX2)
		const __m256i vweights = _mm256_broadcastsi128_si256(mweights);

		const __m256i vhalf = _mm256_set1_epi16(32);
		const __m256i vsign = _mm256_set1_epi16(-0x8000);
		const __m128i mfix4 = _mm_slli_epi32(mfix, 2);

		mc = _mm_packus_epi16(mc, mc);
		__m256i vc = _mm256_broadcastsi128_si256(mc);

		__m256i vt = *(const __m256i*)gTableInterpolate2_U8;

		vt = _mm256_maddubs_epi16(vc, vt);

		vt = _mm256_add_epi16(vt, vhalf);

		vt = _mm256_srli_epi16(vt, 6);

		__m256i vtx = _mm256_permute4x64_epi64(vt, 0x44);
		__m256i vty = _mm256_permute4x64_epi64(vt, 0xEE);

		int k = static_cast<int>(area.Active);
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

		for (size_t i = 0, n = area.Active; i < n; i++)
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
	static INLINED int ComputeSubsetError2Pair(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
	{
		__m128i merrorBlock = _mm_setzero_si128();

#if defined(OPTION_AVX2)
		const __m256i vweights = _mm256_broadcastsi128_si256(mweights);

		const __m256i vhalf = _mm256_set1_epi16(32);
		const __m256i vsign = _mm256_set1_epi16(-0x8000);
		const __m128i mfix4 = _mm_slli_epi32(mfix, 2);

		mc = _mm_shuffle_epi32(mc, shuffle);
		mc = _mm_packus_epi16(mc, mc);
		__m256i vc = _mm256_broadcastsi128_si256(mc);

		__m256i vt = *(const __m256i*)gTableInterpolate2GR_U8;

		vt = _mm256_maddubs_epi16(vc, vt);

		vt = _mm256_add_epi16(vt, vhalf);

		vt = _mm256_srli_epi16(vt, 6);

		int k = static_cast<int>(area.Active);
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

		for (size_t i = 0, n = area.Active; i < n; i++)
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

	static NOTINLINED BlockError ComputeSubsetTable2(const Area& area, __m128i mc, uint64_t& indices) noexcept
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

		Modulations state;
		_mm_store_si128((__m128i*)&state.Values_I16[0], mtx);
		_mm_store_si128((__m128i*)&state.Values_I16[2], mty);

		int error = ComputeSubsetTable(area, gWeightsAGRB, gFixWeightsAGRB, state, 4);

		for (size_t i = 0, n = area.Count; i < n; i++)
		{
			uint64_t index = static_cast<uint32_t>(state.Best[i]);
			indices |= index << (area.Indices[i] << 1);
		}

		int errorAlpha = 0;
		for (size_t i = 0, n = area.Count; i < n; i++)
		{
			int da = *(const uint16_t*)&state.Values_I16[state.Best[i]] - *(const uint16_t*)&area.DataMask_I16[i];

			errorAlpha += da * da;
		}

		return BlockError(errorAlpha * kAlpha, error);
	}

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept
	{
		const size_t partitionIndex = input.BestParameter;

		Area& area1 = GetArea(input.Area12[partitionIndex], input.LazyArea12[partitionIndex], input, gTableSelection12[partitionIndex]);
		Area& area2 = GetArea(input.Area22[partitionIndex], input.LazyArea22[partitionIndex], input, gTableSelection22[partitionIndex]);

		__m128i mc0 = input.BestColor0;
		__m128i mc1 = input.BestColor1;

		uint64_t indices = 0;
		input.Error = ComputeSubsetTable2(area1, mc0, indices);
		input.Error.Add(ComputeSubsetTable2(area2, mc1, indices));

		AreaReduceTable2(area1, mc0, indices);
		AreaReduceTable2(area2, mc1, indices);

		ComposeBlock(output, mc0, mc1, indices, partitionIndex);
	}

	static INLINED int CompressSubsetFast(const Area& area, __m128i& mc, int water) noexcept
	{
		mc = area.Bounds_U16;

		const __m128i mp = _mm_set1_epi16(4);
		__m128i mpbits = _mm_shuffle_epi32(_mm_and_si128(mp, mc), 0);

		const __m128i mh5 = _mm_set1_epi16(0xF8);
		mc = _mm_and_si128(mc, mh5);

		mc = _mm_or_si128(mc, _mm_srli_epi16(mc, 6));

		mc = _mm_or_si128(mc, mpbits);

#if defined(OPTION_COUNTERS)
		gComputeSubsetError2++;
#endif

		const int ea = ComputeSubsetTransparentError2(area, _mm_extract_epi16(_mm_packus_epi16(mc, mc), 0));
		if (ea >= water)
			return water;

		return ComputeSubsetError2(area, mc, gWeightsAGRB, gFixWeightsAGRB, _mm_cvtsi32_si128(water - ea)) + ea;
	}

	void CompressBlockFast(Cell& input) noexcept
	{
		for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
		{
			__m128i mc0 = _mm_setzero_si128();
			__m128i mc1 = _mm_setzero_si128();

			int error = 0;
			if (error < input.Error.Total)
			{
				Area& area1 = GetArea(input.Area12[partitionIndex], input.LazyArea12[partitionIndex], input, gTableSelection12[partitionIndex]);

				error += CompressSubsetFast(area1, mc0, input.Error.Total - error);

				if (error < input.Error.Total)
				{
					Area& area2 = GetArea(input.Area22[partitionIndex], input.LazyArea22[partitionIndex], input, gTableSelection22[partitionIndex]);

					error += CompressSubsetFast(area2, mc1, input.Error.Total - error);

					if (input.Error.Total > error)
					{
						input.Error.Total = error;

						input.BestColor0 = mc0;
						input.BestColor1 = mc1;
						input.BestParameter = partitionIndex;
						input.BestMode = 7;

						if (error <= 0)
							return;
					}
				}
			}
		}
	}

	struct Estimation
	{
		int ch0, ch1, ch2, ch3;
	};

	class Subset final
	{
	public:
		LevelsBuffer<LevelsCapacity> ch0, ch1, ch2, ch3;

		ALWAYS_INLINED Subset() noexcept = default;

		template<int pbits>
		INLINED bool InitLevels(const Area& area, const int water, const Estimation& estimation) noexcept
		{
			if (area.IsOpaque)
			{
				ch0.SetZeroError(0xFFFF);
			}
			else
			{
				ch0.ComputeChannelLevelsReduced<5, pbits, false, gTableDeltas2_Value6>(area, 0, kAlpha, water - estimation.ch1 - estimation.ch2 - estimation.ch3);
			}
			int min0 = ch0.MinErr;
			if (min0 >= water)
				return false;

			ch1.ComputeChannelLevelsReduced<5, pbits, true, gTableDeltas2_Value6>(area, 1, kGreen, water - min0 - estimation.ch2 - estimation.ch3);
			int min1 = ch1.MinErr;
			if (min0 + min1 >= water)
				return false;

			ch2.ComputeChannelLevelsReduced<5, pbits, true, gTableDeltas2_Value6>(area, 2, kRed, water - min0 - min1 - estimation.ch3);
			int min2 = ch2.MinErr;
			if (min0 + min1 + min2 >= water)
				return false;

			ch3.ComputeChannelLevelsReduced<5, pbits, true, gTableDeltas2_Value6>(area, 3, kBlue, water - min0 - min1 - min2);
			int min3 = ch3.MinErr;
			if (min0 + min1 + min2 + min3 >= water)
				return false;

			return true;
		}

		INLINED int TryVariants(const Area& area, __m128i& best_color, int water) noexcept
		{
			int min0 = ch0.MinErr;
			int min1 = ch1.MinErr;
			int min2 = ch2.MinErr;
			int min3 = ch3.MinErr;
			if (min0 + min1 + min2 + min3 >= water)
				return water;

			int n0 = ch0.Count;
			int n1 = ch1.Count;
			int n2 = ch2.Count;
			int n3 = ch3.Count;

			int memAR[LevelsCapacity];
			int memAGB[LevelsCapacity];

			for (int i0 = 0; i0 < n0; i0++)
			{
				int e0 = ch0.Err[i0].Error;
				if (e0 + min1 + min2 + min3 >= water)
					break;

				int c0 = ch0.Err[i0].Color;

				const int ea = ComputeSubsetTransparentError2(area, c0);

				for (int i = 0; i < n2; i++)
				{
					memAR[i] = -1;
				}

				for (int i1 = 0; i1 < n1; i1++)
				{
					int e1 = ch1.Err[i1].Error + e0;
					if (e1 + min2 + min3 >= water)
						break;

					int c1 = ch1.Err[i1].Color;

					if (!area.IsOpaque)
					{
						__m128i mc = _mm_setzero_si128();
						mc = _mm_insert_epi16(mc, c0, 0);
						mc = _mm_insert_epi16(mc, c1, 1);
						mc = _mm_cvtepu8_epi16(mc);

#if defined(OPTION_COUNTERS)
						gComputeSubsetError2AG++;
#endif
						e1 = ComputeSubsetError2Pair<_MM_SHUFFLE(1, 0, 1, 0)>(area, mc, gWeightsAGAG, gFixWeightsAG, _mm_cvtsi32_si128(water - ea - min2 - min3)) + ea;
						if (e1 + min2 + min3 >= water)
							continue;
					}

					for (int i = 0; i < n3; i++)
					{
						memAGB[i] = -1;
					}

					for (int i2 = 0; i2 < n2; i2++)
					{
						int e2 = ch2.Err[i2].Error + e1;
						if (e2 + min3 >= water)
							break;

						int c2 = ch2.Err[i2].Color;

						if (!area.IsOpaque)
						{
							int ear = memAR[i2];
							if (ear < 0)
							{
								__m128i mc = _mm_setzero_si128();
								mc = _mm_insert_epi16(mc, c0, 0);
								mc = _mm_insert_epi16(mc, c2, 2);
								mc = _mm_cvtepu8_epi16(mc);

#if defined(OPTION_COUNTERS)
								gComputeSubsetError2AR++;
#endif
								ear = ComputeSubsetError2Pair<_MM_SHUFFLE(2, 0, 2, 0)>(area, mc, gWeightsARAR, gFixWeightsAR, _mm_cvtsi32_si128(water - ea - min1 - min3)) + ea;
								memAR[i2] = ear;
							}
							if (ear + min1 + min3 >= water)
								continue;
						}

						{
							__m128i mc = _mm_setzero_si128();
							mc = _mm_insert_epi16(mc, c0, 0);
							mc = _mm_insert_epi16(mc, c1, 1);
							mc = _mm_insert_epi16(mc, c2, 2);
							mc = _mm_cvtepu8_epi16(mc);

#if defined(OPTION_COUNTERS)
							gComputeSubsetError2AGR++;
#endif
							e2 = ComputeSubsetError2(area, mc, gWeightsAGR, gFixWeightsAGR, _mm_cvtsi32_si128(water - ea - min3)) + ea;
							if (e2 + min3 >= water)
								continue;
						}

						for (int i3 = 0; i3 < n3; i3++)
						{
							int e3 = ch3.Err[i3].Error + e2;
							if (e3 >= water)
								break;

							int c3 = ch3.Err[i3].Color;

							int eagb = memAGB[i3];
							if (eagb < 0)
							{
								__m128i mc = _mm_setzero_si128();
								mc = _mm_insert_epi16(mc, c0, 0);
								mc = _mm_insert_epi16(mc, c1, 1);
								mc = _mm_insert_epi16(mc, c3, 3);
								mc = _mm_cvtepu8_epi16(mc);

#if defined(OPTION_COUNTERS)
								gComputeSubsetError2AGB++;
#endif
								eagb = ComputeSubsetError2(area, mc, gWeightsAGB, gFixWeightsAGB, _mm_cvtsi32_si128(water - ea - min2)) + ea;
								memAGB[i3] = eagb;
							}
							if (eagb + min2 >= water)
								continue;

							__m128i mc = _mm_setzero_si128();
							mc = _mm_insert_epi16(mc, c0, 0);
							mc = _mm_insert_epi16(mc, c1, 1);
							mc = _mm_insert_epi16(mc, c2, 2);
							mc = _mm_insert_epi16(mc, c3, 3);
							mc = _mm_cvtepu8_epi16(mc);

#if defined(OPTION_COUNTERS)
							gComputeSubsetError2++;
#endif
							int err = ComputeSubsetError2(area, mc, gWeightsAGRB, gFixWeightsAGRB, _mm_cvtsi32_si128(water - ea)) + ea;

							if (water > err)
							{
								water = err;

								best_color = mc;
							}
						}
					}
				}
			}

			return water;
		}
	};

	class Subsets final
	{
	public:
		Subset subset3;
		Subset subset0;
		Subset subset1;
		Subset subset2;

		bool valid3 = false;
		bool valid0 = false;
		bool valid1 = false;
		bool valid2 = false;

		ALWAYS_INLINED Subsets() noexcept = default;

		INLINED bool InitLevels(const Area& area, const int water, const Estimation& estimation) noexcept
		{
			valid3 = subset3.InitLevels<0x0101>(area, water, estimation);

			if (!area.IsOpaque)
			{
				valid0 = subset0.InitLevels<0>(area, water, estimation);
				valid1 = subset1.InitLevels<1>(area, water, estimation);
				valid2 = subset2.InitLevels<1 << 8>(area, water, estimation);
			}

			return valid3 || valid0 || valid1 || valid2;
		}

		INLINED int TryVariants(const Area& area, __m128i& best_color, int water) noexcept
		{
			if (valid3)
			{
				water = subset3.TryVariants(area, best_color, water);
			}

			if (valid0)
			{
				water = subset0.TryVariants(area, best_color, water);
			}

			if (valid1)
			{
				water = subset1.TryVariants(area, best_color, water);
			}

			if (valid2)
			{
				water = subset2.TryVariants(area, best_color, water);
			}

			return water;
		}
	};

	void CompressBlock(Cell& input) noexcept
	{
		const size_t partitionIndex = input.BestParameter;

		__m128i mc0 = input.BestColor0;
		__m128i mc1 = input.BestColor1;

		const Estimation estimation{ 0, 0, 0, 0 };

		int error = 0;
		{
			Area& area1 = GetArea(input.Area12[partitionIndex], input.LazyArea12[partitionIndex], input, gTableSelection12[partitionIndex]);

#if defined(OPTION_COUNTERS)
			gComputeSubsetError2++;
#endif
			int water1 = ComputeSubsetError2(area1, mc0, gWeightsAGRB, gFixWeightsAGRB, _mm_cvtsi32_si128(kBlockMaximalAlphaError + kBlockMaximalColorError));
			water1 += ComputeSubsetTransparentError2(area1, _mm_extract_epi16(_mm_packus_epi16(mc0, mc0), 0));
			if (water1)
			{
				Subsets subsets1;
				if (subsets1.InitLevels(area1, water1, estimation))
				{
					error += subsets1.TryVariants(area1, mc0, water1);
				}
				else
				{
					error += water1;
				}
			}

			Area& area2 = GetArea(input.Area22[partitionIndex], input.LazyArea22[partitionIndex], input, gTableSelection22[partitionIndex]);

#if defined(OPTION_COUNTERS)
			gComputeSubsetError2++;
#endif
			int water2 = ComputeSubsetError2(area2, mc1, gWeightsAGRB, gFixWeightsAGRB, _mm_cvtsi32_si128(kBlockMaximalAlphaError + kBlockMaximalColorError));
			water2 += ComputeSubsetTransparentError2(area2, _mm_extract_epi16(_mm_packus_epi16(mc1, mc1), 0));
			if (water2)
			{
				Subsets subsets2;
				if (subsets2.InitLevels(area2, water2, estimation))
				{
					error += subsets2.TryVariants(area2, mc1, water2);
				}
				else
				{
					error += water2;
				}
			}

			if (input.Error.Total > error)
			{
				input.Error.Total = error;

				input.BestColor0 = mc0;
				input.BestColor1 = mc1;
				//input.BestParameter = partitionIndex;
				//input.BestMode = 7;
			}
		}
	}

	static INLINED int EstimateBest(Area& area) noexcept
	{
#if defined(OPTION_PCA)
		if (area.IsOpaque)
		{
			return AreaGetBestPca3(area);
		}

		return PrincipalComponentAnalysis4(area);
#else
		(void)area;
		return 0;
#endif
	}

	static INLINED int EstimateLevels(const Area& area, const int water, Estimation& estimation) noexcept
	{
		int error = 0;
		if (error < water)
		{
			if (!area.Active)
			{
				estimation.ch0 = 0;
				estimation.ch1 = 0;
				estimation.ch2 = 0;
				estimation.ch3 = 0;
				return 0;
			}

			if (!area.IsOpaque)
			{
				int level0 = LevelsMinimum::EstimateChannelLevelsReduced<6, false, gTableDeltas2_Value6, gTableCuts2_Value6>(area, 0, kAlpha, water - error);
				estimation.ch0 = level0;
				error += level0;
			}
			else
			{
				estimation.ch0 = 0;
			}

			if (error < water)
			{
				int level1 = LevelsMinimum::EstimateChannelLevelsReduced<6, true, gTableDeltas2_Value6, gTableCuts2_Value6>(area, 1, kGreen, water - error);
				estimation.ch1 = level1;
				error += level1;

				if (error < water)
				{
					int level2 = LevelsMinimum::EstimateChannelLevelsReduced<6, true, gTableDeltas2_Value6, gTableCuts2_Value6>(area, 2, kRed, water - error);
					estimation.ch2 = level2;
					error += level2;

					if (error < water)
					{
						int level3 = LevelsMinimum::EstimateChannelLevelsReduced<6, true, gTableDeltas2_Value6, gTableCuts2_Value6>(area, 3, kBlue, water - error);
						estimation.ch3 = level3;
						error += level3;

						if (error < water)
						{
							return error;
						}
					}
				}
			}
		}

		return water;
	}

	void CompressBlockFull(Cell& input) noexcept
	{
		Node order[64]{};
		int lines1[64];
		int lines2[64];

		size_t partitionsCount = 0;
		for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
		{
			if ((input.PersonalMode == 7) && (input.PersonalParameter == partitionIndex))
				continue;

			Area& area1 = GetArea(input.Area12[partitionIndex], input.LazyArea12[partitionIndex], input, gTableSelection12[partitionIndex]);

			const int water1 = input.Error.Total;
			int line1 = EstimateBest(area1);
			if (line1 < water1)
			{
				Area& area2 = GetArea(input.Area22[partitionIndex], input.LazyArea22[partitionIndex], input, gTableSelection22[partitionIndex]);

				const int water2 = water1 - line1;
				int line2 = EstimateBest(area2);
				if (line2 < water2)
				{
					lines1[partitionIndex] = line1;
					lines2[partitionIndex] = line2;

					order[partitionsCount++].Init(line1 + line2, static_cast<int>(partitionIndex));
				}
			}
		}

		{
			Node order2[64];
			if (radix_sort(order, order2, partitionsCount) == order2)
			{
				memcpy(order, order2, sizeof(Node) * partitionsCount);
			}
		}

		for (size_t i = 0; i < partitionsCount; i++)
		{
			if (order[i].Error < input.Error.Total)
			{
				const size_t partitionIndex = order[i].Color;

				__m128i mc0 = _mm_setzero_si128();
				__m128i mc1 = _mm_setzero_si128();

				Area& area1 = GetArea(input.Area12[partitionIndex], input.LazyArea12[partitionIndex], input, gTableSelection12[partitionIndex]);
				Area& area2 = GetArea(input.Area22[partitionIndex], input.LazyArea22[partitionIndex], input, gTableSelection22[partitionIndex]);

				int line1 = lines1[partitionIndex];
				int line2 = lines2[partitionIndex];

				Estimation estimations1;
				int water1 = input.Error.Total - line2;
				line1 = Max(line1, EstimateLevels(area1, water1, estimations1));
				if (line1 < water1)
				{
					Estimation estimations2;
					int water2 = water1 - line1;
					line2 = Max(line2, EstimateLevels(area2, water2, estimations2));
					if (line2 < water2)
					{
						water1 = input.Error.Total - line2;
						Subsets subsets1;
						if (subsets1.InitLevels(area1, water1, estimations1))
						{
							int error = subsets1.TryVariants(area1, mc0, water1);
							if (error < water1)
							{
								water2 = input.Error.Total - error;
								Subsets subsets2;
								if (subsets2.InitLevels(area2, water2, estimations2))
								{
									error += subsets2.TryVariants(area2, mc1, water2);

									if (input.Error.Total > error)
									{
										input.Error.Total = error;

										input.BestColor0 = mc0;
										input.BestColor1 = mc1;
										input.BestParameter = partitionIndex;
										input.BestMode = 7;

										if (error <= 0)
											return;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	void PrintCounters() noexcept
	{
#if defined(OPTION_COUNTERS)
		PRINTF("[Mode 7]\tAG2 = %i, AR2 = %i, AGR2 = %i, AGB2 = %i, AGRB2 = %i",
			gComputeSubsetError2AG.load(), gComputeSubsetError2AR.load(),
			gComputeSubsetError2AGR.load(), gComputeSubsetError2AGB.load(),
			gComputeSubsetError2.load());
#endif
	}

} // namespace Mode7
