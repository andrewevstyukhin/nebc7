
#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"

#include "SnippetInsertRemoveZeroBit.h"
#include "SnippetHorizontalSum4.h"
#include "SnippetLevelsBuffer.h"

// https://docs.microsoft.com/en-us/windows/desktop/direct3d11/bc7-format-mode-reference#mode-6

namespace Mode6 {

#if defined(OPTION_COUNTERS)
	static std::atomic_int gComputeSubsetError4, gComputeSubsetError4AG, gComputeSubsetError4AR, gComputeSubsetError4AGR, gComputeSubsetError4AGB;
#endif

	static INLINED void DecompressSubset(__m128i mc, int* output, uint64_t data) noexcept
	{
		const __m128i mhalf = _mm_set1_epi16(32);

		mc = _mm_packus_epi16(mc, mc);

		for (size_t i = 0; i < 16; i++)
		{
			__m128i mratio = _mm_loadl_epi64((const __m128i*)&((const uint64_t*)gTableInterpolate4_U8)[data & 0xF]); data >>= 4;

			__m128i mx = _mm_maddubs_epi16(mc, mratio);
			mx = _mm_add_epi16(mx, mhalf);
			mx = _mm_srli_epi16(mx, 6);

			output[i] = _mm_cvtsi128_si32(_mm_packus_epi16(mx, mx));
		}
	}

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept
	{
		uint64_t data0 = *(const uint64_t*)&input[0];
		uint64_t data1 = *(const uint64_t*)&input[8];

		data0 >>= 7;

		__m128i mc0 = _mm_setzero_si128();

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 4); data0 >>= 7;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 5); data0 >>= 7;

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 2); data0 >>= 7;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 3); data0 >>= 7;

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 6); data0 >>= 7;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 7); data0 >>= 7;

		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 0); data0 >>= 7;
		mc0 = _mm_insert_epi16(mc0, static_cast<int>(data0), 1); data0 >>= 7;

		const __m128i m7 = _mm_set1_epi16(0x7F);
		mc0 = _mm_and_si128(mc0, m7);

		mc0 = _mm_add_epi16(mc0, mc0);

		int pbits = static_cast<int>(data0) | (static_cast<int>(data1 & 1) << 16); data1 >>= 1;
		mc0 = _mm_or_si128(mc0, _mm_shuffle_epi32(_mm_cvtsi32_si128(pbits), 0));

		data1 = InsertZeroBit(data1, 3);

		DecompressSubset(mc0, (int*)output.ImageRows_U8, data1);

		output.BestParameter = 0;
		output.BestMode = 6;
	}

	static INLINED void ComposeBlock(uint8_t output[16], __m128i mc0, uint64_t indices) noexcept
	{
		uint64_t data1 = RemoveZeroBit(indices, 3);

		data1 <<= 1; data1 |= _mm_extract_epi16(mc0, 1) & 1;

		uint64_t data0 = _mm_extract_epi16(mc0, 0) & 1;

		mc0 = _mm_srli_epi16(mc0, 1);

		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 1);
		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 0);

		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 7);
		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 6);

		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 3);
		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 2);

		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 5);
		data0 <<= 7; data0 |= _mm_extract_epi16(mc0, 4);

		data0 <<= 7; data0 |= 1 << 6;

		*(uint64_t*)&output[0] = data0;
		*(uint64_t*)&output[8] = data1;
	}

	static INLINED int ComputeSubsetTransparentError4(const Area& area, const int alpha) noexcept
	{
		int error = static_cast<int>(area.Count - area.Active);
		if (error)
		{
			error *= kAlpha;
			error *= gTableLevels4_U16[0][alpha];
		}

		return error;
	}

	static INLINED int ComputeSubsetError4(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
	{
		__m128i merrorBlock = _mm_setzero_si128();

		const __m128i mhalf = _mm_set1_epi16(32);
		const __m128i msign = _mm_set1_epi16(-0x8000);

		mc = _mm_packus_epi16(mc, mc);

		__m128i mtx = gTableInterpolate4_U8[0];
		__m128i mty = gTableInterpolate4_U8[1];
		__m128i mtz = gTableInterpolate4_U8[2];
		__m128i mtw = gTableInterpolate4_U8[3];
		__m128i mrx = gTableInterpolate4_U8[4];
		__m128i mry = gTableInterpolate4_U8[5];
		__m128i mrz = gTableInterpolate4_U8[6];
		__m128i mrw = gTableInterpolate4_U8[7];

		mtx = _mm_maddubs_epi16(mc, mtx);
		mty = _mm_maddubs_epi16(mc, mty);
		mtz = _mm_maddubs_epi16(mc, mtz);
		mtw = _mm_maddubs_epi16(mc, mtw);
		mrx = _mm_maddubs_epi16(mc, mrx);
		mry = _mm_maddubs_epi16(mc, mry);
		mrz = _mm_maddubs_epi16(mc, mrz);
		mrw = _mm_maddubs_epi16(mc, mrw);

		mtx = _mm_add_epi16(mtx, mhalf);
		mty = _mm_add_epi16(mty, mhalf);
		mtz = _mm_add_epi16(mtz, mhalf);
		mtw = _mm_add_epi16(mtw, mhalf);
		mrx = _mm_add_epi16(mrx, mhalf);
		mry = _mm_add_epi16(mry, mhalf);
		mrz = _mm_add_epi16(mrz, mhalf);
		mrw = _mm_add_epi16(mrw, mhalf);

		mtx = _mm_srli_epi16(mtx, 6);
		mty = _mm_srli_epi16(mty, 6);
		mtz = _mm_srli_epi16(mtz, 6);
		mtw = _mm_srli_epi16(mtw, 6);
		mrx = _mm_srli_epi16(mrx, 6);
		mry = _mm_srli_epi16(mry, 6);
		mrz = _mm_srli_epi16(mrz, 6);
		mrw = _mm_srli_epi16(mrw, 6);

		for (size_t i = 0, n = area.Active; i < n; i++)
		{
			__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
			__m128i mpixel = _mm_unpacklo_epi64(mpacked, mpacked);
			__m128i mmask = _mm_unpackhi_epi64(mpacked, mpacked);

			__m128i mx = _mm_sub_epi16(mpixel, mtx);
			__m128i my = _mm_sub_epi16(mpixel, mty);
			__m128i mz = _mm_sub_epi16(mpixel, mtz);
			__m128i mw = _mm_sub_epi16(mpixel, mtw);

			mx = _mm_mullo_epi16(mx, mx);
			my = _mm_mullo_epi16(my, my);
			mz = _mm_mullo_epi16(mz, mz);
			mw = _mm_mullo_epi16(mw, mw);

			mx = _mm_and_si128(mx, mmask);
			my = _mm_and_si128(my, mmask);
			mz = _mm_and_si128(mz, mmask);
			mw = _mm_and_si128(mw, mmask);

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
			__m128i mx0 = _mm_min_epi32(mx, mz);

			merrorBlock = _mm_add_epi32(merrorBlock, mfix);

			mx = _mm_sub_epi16(mpixel, mrx);
			my = _mm_sub_epi16(mpixel, mry);
			mz = _mm_sub_epi16(mpixel, mrz);
			mw = _mm_sub_epi16(mpixel, mrw);

			mx = _mm_mullo_epi16(mx, mx);
			my = _mm_mullo_epi16(my, my);
			mz = _mm_mullo_epi16(mz, mz);
			mw = _mm_mullo_epi16(mw, mw);

			mx = _mm_and_si128(mx, mmask);
			my = _mm_and_si128(my, mmask);
			mz = _mm_and_si128(mz, mmask);
			mw = _mm_and_si128(mw, mmask);

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
			mx = _mm_min_epi32(mx, mx0);
			mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

			merrorBlock = _mm_add_epi32(merrorBlock, mx);

			if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
				break;
		}

		return _mm_cvtsi128_si32(merrorBlock);
	}

	template<int shuffle>
	static INLINED int ComputeSubsetError4Pair(const Area& area, __m128i mc, const __m128i mweights, const __m128i mfix, const __m128i mwater) noexcept
	{
		__m128i merrorBlock = _mm_setzero_si128();

		const __m128i mhalf = _mm_set1_epi16(32);
		const __m128i msign = _mm_set1_epi16(-0x8000);

		mc = _mm_shuffle_epi32(mc, shuffle);
		mc = _mm_packus_epi16(mc, mc);

		__m128i mtx = gTableInterpolate4GR_U8[0];
		__m128i mty = gTableInterpolate4GR_U8[1];
		__m128i mtz = gTableInterpolate4GR_U8[2];
		__m128i mtw = gTableInterpolate4GR_U8[3];

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

		for (size_t i = 0, n = area.Active; i < n; i++)
		{
			__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
			__m128i mpixel = _mm_shufflelo_epi16(mpacked, shuffle);
			__m128i mmask = _mm_shufflehi_epi16(mpacked, shuffle);
			mpixel = _mm_unpacklo_epi64(mpixel, mpixel);
			mmask = _mm_unpackhi_epi64(mmask, mmask);

			merrorBlock = _mm_add_epi32(merrorBlock, mfix);

			__m128i mx = _mm_sub_epi16(mpixel, mtx);
			__m128i my = _mm_sub_epi16(mpixel, mty);
			__m128i mz = _mm_sub_epi16(mpixel, mtz);
			__m128i mw = _mm_sub_epi16(mpixel, mtw);

			mx = _mm_mullo_epi16(mx, mx);
			my = _mm_mullo_epi16(my, my);
			mz = _mm_mullo_epi16(mz, mz);
			mw = _mm_mullo_epi16(mw, mw);

			mx = _mm_and_si128(mx, mmask);
			my = _mm_and_si128(my, mmask);
			mz = _mm_and_si128(mz, mmask);
			mw = _mm_and_si128(mw, mmask);

			mx = _mm_xor_si128(mx, msign);
			my = _mm_xor_si128(my, msign);
			mz = _mm_xor_si128(mz, msign);
			mw = _mm_xor_si128(mw, msign);

			mx = _mm_madd_epi16(mx, mweights);
			my = _mm_madd_epi16(my, mweights);
			mz = _mm_madd_epi16(mz, mweights);
			mw = _mm_madd_epi16(mw, mweights);

			mx = _mm_min_epi32(mx, my);
			mz = _mm_min_epi32(mz, mw);
			mx = _mm_min_epi32(mx, mz);
			mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));
			mx = _mm_min_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(1, 0, 3, 2)));

			merrorBlock = _mm_add_epi32(merrorBlock, mx);

			if (!(_mm_movemask_epi8(_mm_cmpgt_epi32(mwater, merrorBlock)) & 0xF))
				break;
		}

		return _mm_cvtsi128_si32(merrorBlock);
	}

	static INLINED BlockError ComputeSubsetTable4(const Area& area, __m128i mc, uint64_t& indices) noexcept
	{
		const __m128i mhalf = _mm_set1_epi16(32);

		mc = _mm_packus_epi16(mc, mc);

		__m128i mtx = gTableInterpolate4_U8[0];
		__m128i mty = gTableInterpolate4_U8[1];
		__m128i mtz = gTableInterpolate4_U8[2];
		__m128i mtw = gTableInterpolate4_U8[3];
		__m128i mrx = gTableInterpolate4_U8[4];
		__m128i mry = gTableInterpolate4_U8[5];
		__m128i mrz = gTableInterpolate4_U8[6];
		__m128i mrw = gTableInterpolate4_U8[7];

		mtx = _mm_maddubs_epi16(mc, mtx);
		mty = _mm_maddubs_epi16(mc, mty);
		mtz = _mm_maddubs_epi16(mc, mtz);
		mtw = _mm_maddubs_epi16(mc, mtw);
		mrx = _mm_maddubs_epi16(mc, mrx);
		mry = _mm_maddubs_epi16(mc, mry);
		mrz = _mm_maddubs_epi16(mc, mrz);
		mrw = _mm_maddubs_epi16(mc, mrw);

		mtx = _mm_add_epi16(mtx, mhalf);
		mty = _mm_add_epi16(mty, mhalf);
		mtz = _mm_add_epi16(mtz, mhalf);
		mtw = _mm_add_epi16(mtw, mhalf);
		mrx = _mm_add_epi16(mrx, mhalf);
		mry = _mm_add_epi16(mry, mhalf);
		mrz = _mm_add_epi16(mrz, mhalf);
		mrw = _mm_add_epi16(mrw, mhalf);

		mtx = _mm_srli_epi16(mtx, 6);
		mty = _mm_srli_epi16(mty, 6);
		mtz = _mm_srli_epi16(mtz, 6);
		mtw = _mm_srli_epi16(mtw, 6);
		mrx = _mm_srli_epi16(mrx, 6);
		mry = _mm_srli_epi16(mry, 6);
		mrz = _mm_srli_epi16(mrz, 6);
		mrw = _mm_srli_epi16(mrw, 6);

		Modulations state;
		_mm_store_si128((__m128i*)&state.Values_I16[0], mtx);
		_mm_store_si128((__m128i*)&state.Values_I16[2], mty);
		_mm_store_si128((__m128i*)&state.Values_I16[4], mtz);
		_mm_store_si128((__m128i*)&state.Values_I16[6], mtw);
		_mm_store_si128((__m128i*)&state.Values_I16[8], mrx);
		_mm_store_si128((__m128i*)&state.Values_I16[10], mry);
		_mm_store_si128((__m128i*)&state.Values_I16[12], mrz);
		_mm_store_si128((__m128i*)&state.Values_I16[14], mrw);

		int error = ComputeSubsetTable(area, gWeightsAGRB, gFixWeightsAGRB, state, 16);

		for (size_t i = 0; i < 16; i++)
		{
			uint64_t index = static_cast<uint32_t>(state.Best[i]);
			indices |= index << (area.Indices[i] << 2);
		}

		int errorAlpha = 0;
		for (size_t i = 0; i < 16; i++)
		{
			int da = *(const uint16_t*)&state.Values_I16[state.Best[i]] - *(const uint16_t*)&area.DataMask_I16[i];

			errorAlpha += da * da;
		}

		return BlockError(errorAlpha * kAlpha, error);
	}

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept
	{
		Area& area = input.Area1;

		__m128i mc = input.BestColor0;

		uint64_t indices = 0;
		input.Error = ComputeSubsetTable4(area, mc, indices);
		AreaReduceTable4(mc, indices);

		ComposeBlock(output, mc, indices);
	}

	static INLINED int CompressSubsetFast(const Area& area, __m128i& mc, int water) noexcept
	{
		mc = area.Bounds_U16;

		const __m128i m1 = _mm_set1_epi16(1);
		__m128i mpbits = _mm_shuffle_epi32(_mm_and_si128(m1, mc), 0);

		const __m128i mh7 = _mm_set1_epi16(0xFE);
		mc = _mm_and_si128(mc, mh7);

		mc = _mm_or_si128(mc, mpbits);

#if defined(OPTION_COUNTERS)
		gComputeSubsetError4++;
#endif

		const int ea = ComputeSubsetTransparentError4(area, _mm_extract_epi16(_mm_packus_epi16(mc, mc), 0));
		if (ea >= water)
			return water;

		return ComputeSubsetError4(area, mc, gWeightsAGRB, gFixWeightsAGRB, _mm_cvtsi32_si128(water - ea)) + ea;
	}

	void CompressBlockFast(Cell& input) noexcept
	{
		__m128i mc = _mm_setzero_si128();

		int error = 0;
		if (error < input.Error.Total)
		{
			Area& area = input.Area1;

			error += CompressSubsetFast(area, mc, input.Error.Total - error);
		}

		if (input.Error.Total > error)
		{
			input.Error.Total = error;

			input.BestColor0 = mc;
			input.BestMode = 6;
		}
	}

	class Subset final
	{
	public:
		LevelsBuffer<20> ch0, ch1, ch2, ch3;

		INLINED Subset() noexcept = default;

		template<int pbits>
		INLINED bool InitLevels(const Area& area, const int water) noexcept
		{
			if (area.IsOpaque)
			{
				ch0.SetZeroError(0xFFFF);
			}
			else
			{
				ch0.ComputeChannelLevelsReduced<7, pbits, false, gTableLevels4_U16>(area, 0, kAlpha, water);
			}
			int min0 = ch0.MinErr;
			if (min0 >= water)
				return false;

			ch1.ComputeChannelLevelsReduced<7, pbits, true, gTableLevels4_U16>(area, 1, kGreen, water - min0);
			int min1 = ch1.MinErr;
			if (min0 + min1 >= water)
				return false;

			ch2.ComputeChannelLevelsReduced<7, pbits, true, gTableLevels4_U16>(area, 2, kRed, water - min0 - min1);
			int min2 = ch2.MinErr;
			if (min0 + min1 + min2 >= water)
				return false;

			ch3.ComputeChannelLevelsReduced<7, pbits, true, gTableLevels4_U16>(area, 3, kBlue, water - min0 - min1 - min2);
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

			int memAR[20];
			int memAGB[20];

			for (int i0 = 0; i0 < n0; i0++)
			{
				int e0 = ch0.Err[i0].Error;
				if (e0 + min1 + min2 + min3 >= water)
					break;

				int c0 = ch0.Err[i0].Color;

				const int ea = ComputeSubsetTransparentError4(area, c0);

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
						gComputeSubsetError4AG++;
#endif
						e1 = ComputeSubsetError4Pair<_MM_SHUFFLE(1, 0, 1, 0)>(area, mc, gWeightsAGAG, gFixWeightsAG, _mm_cvtsi32_si128(water - ea - min2 - min3)) + ea;
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
								gComputeSubsetError4AR++;
#endif
								ear = ComputeSubsetError4Pair<_MM_SHUFFLE(2, 0, 2, 0)>(area, mc, gWeightsARAR, gFixWeightsAR, _mm_cvtsi32_si128(water - ea - min1 - min3)) + ea;
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
							gComputeSubsetError4AGR++;
#endif
							e2 = ComputeSubsetError4(area, mc, gWeightsAGR, gFixWeightsAGR, _mm_cvtsi32_si128(water - ea - min3)) + ea;
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
								gComputeSubsetError4AGB++;
#endif
								eagb = ComputeSubsetError4(area, mc, gWeightsAGB, gFixWeightsAGB, _mm_cvtsi32_si128(water - ea - min2)) + ea;
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
							gComputeSubsetError4++;
#endif
							int err = ComputeSubsetError4(area, mc, gWeightsAGRB, gFixWeightsAGRB, _mm_cvtsi32_si128(water - ea)) + ea;

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

	static INLINED int CompressSubset(const Area& area, __m128i& mc, int water)
	{
		Subset subset3;
		if (subset3.InitLevels<0x0101>(area, water))
		{
			water = subset3.TryVariants(area, mc, water);
		}

		if (!area.IsOpaque)
		{
			Subset subset0;
			if (subset0.InitLevels<0>(area, water))
			{
				water = subset0.TryVariants(area, mc, water);
			}

			Subset subset1;
			if (subset1.InitLevels<1>(area, water))
			{
				water = subset1.TryVariants(area, mc, water);
			}

			Subset subset2;
			if (subset2.InitLevels<1 << 8>(area, water))
			{
				water = subset2.TryVariants(area, mc, water);
			}
		}

		return water;
	}

	void CompressBlock(Cell& input) noexcept
	{
		__m128i mc = input.BestColor0;

		int error = 0;
		{
			Area& area = input.Area1;

			error += CompressSubset(area, mc, input.Error.Total - error);
		}

		if (input.Error.Total > error)
		{
			input.Error.Total = error;

			input.BestColor0 = mc;
			//input.BestMode = 6;
		}
	}

	void CompressBlockFull(Cell& input) noexcept
	{
		if (input.PersonalMode == 6)
			return;

		__m128i mc = _mm_setzero_si128();

		int error = 0;
		if (error < input.Error.Total)
		{
			Area& area = input.Area1;

			error += CompressSubset(area, mc, input.Error.Total - error);
		}

		if (input.Error.Total > error)
		{
			input.Error.Total = error;

			input.BestColor0 = mc;
			input.BestMode = 6;
		}
	}

	void PrintCounters() noexcept
	{
#if defined(OPTION_COUNTERS)
		PRINTF("[Mode 6]\tAG4 = %i, AR4 = %i, AGR4 = %i, AGB4 = %i, AGRB4 = %i",
			gComputeSubsetError4AG.load(), gComputeSubsetError4AR.load(),
			gComputeSubsetError4AGR.load(), gComputeSubsetError4AGB.load(),
			gComputeSubsetError4.load());
#endif
	}

} // namespace Mode6
