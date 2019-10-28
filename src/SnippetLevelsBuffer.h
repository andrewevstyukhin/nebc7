#pragma once

#include "pch.h"
#include "Bc7Core.h"

#if defined(OPTION_COUNTERS)
inline std::atomic_int gEstimateFull, gEstimateShort;

inline std::atomic_int gLevels[17];
#endif

template<int MaxSize>
class LevelsBuffer final
{
public:
	int MinErr = kBlockMaximalAlphaError + kBlockMaximalColorError;
	int Count = 0;

	Node Err[MaxSize];

	static INLINED int Min(int x, int y) noexcept
	{
		return (x < y) ? x : y;
	}

	static INLINED int Max(int x, int y) noexcept
	{
		return (x > y) ? x : y;
	}

	INLINED LevelsBuffer() noexcept = default;

	INLINED void SetZeroError(int color) noexcept
	{
		MinErr = 0;
		Count = 1;

		Err[0].Error = 0;
		Err[0].Color = color;
	}

	template<int bits, int pbits, bool transparent, const uint16_t table[0x100][(1 << (bits + int(pbits >= 0))) * (1 << (bits + int(pbits >= 0)))]>
	NOTINLINED void ComputeChannelLevelsReduced(const Area& area, const size_t offset, const int weight, const int water, const int QuerySize = MaxSize) noexcept
	{
		const uint16_t* values[16];

		size_t count;
		if constexpr (transparent)
		{
			count = area.Active;

			if (!count)
			{
				if constexpr (pbits < 0)
				{
					SetZeroError(0);
				}
				else
				{
					SetZeroError(pbits << (7 - bits));
				}
				return;
			}
		}
		else
		{
			count = area.Count;
		}

		for (size_t i = 0; i < count; i++)
		{
			size_t value = ((const uint16_t*)&area.DataMask_I16[i])[offset];

			values[i] = table[value];
		}

		int top = (water + weight - 1) / weight;
		if (!top)
			return;

		int d = _mm_cvtsi128_si32(_mm_cvttps_epi32(_mm_sqrt_ss(_mm_cvtepi32_ps(_mm_cvtsi32_si128(top)))));
		d -= int(d * d >= top);

		constexpr int tailmask = (1 << (8 - bits)) - 1;
		constexpr int shift = bits + int(pbits >= 0);

		int L = ((const short*)&area.MinMax_U16)[offset + offset + 0];
		int H = ((const short*)&area.MinMax_U16)[offset + offset + 1];

		const bool reverse = (L != ((const short*)&area.Bounds_U16)[offset + offset + 0]);

		if constexpr ((pbits < 0) || ((pbits >> 8) == (pbits & 1)))
		{
			if (L == H)
			{
				if ((((L >> shift) ^ L) & ((1 << (8 - shift)) - 1)) == 0)
				{
					if constexpr (pbits < 0)
					{
						SetZeroError((H << 8) + L);
						return;
					}
					else if (((L >> (8 - shift)) & 1) == (pbits & 1))
					{
						SetZeroError((H << 8) + L);
						return;
					}
				}
			}
		}

#if defined(OPTION_COUNTERS)
		gLevels[count]++;
#endif

		int LH = Min(L + d, 255);
		int HL = Max(H - d, 0);

		LH = Min(H - (H >> shift) + tailmask, LH - (LH >> shift)) & ~tailmask;
		HL = Max(L - (L >> shift), HL - (HL >> shift) + tailmask) & ~tailmask;

		int HH = 0x100;

		LH >>= 8 - shift;
		HL >>= 8 - shift;
		HH >>= 8 - shift;

		Node nodes1[(1 << shift) * (1 << shift)];
		Node nodes2[(1 << shift) * (1 << shift)];
		Node* nodesPtr = nodes1;

		const __m128i mtop = _mm_shuffle_epi32(_mm_cvtsi32_si128(top), 0);

		if (top <= 0xFFFF)
		{
#if defined(OPTION_COUNTERS)
			gEstimateShort++;
#endif

			if constexpr (pbits < 0)
			{
				for (int iH = HL; iH < HH; iH++)
				{
					int cH = (iH << shift);

#if defined(OPTION_AVX512)
					int hL = LH + cH - 24;
					int iL = cH;
					while (iL <= hL)
					{
						Estimate32Short(nodesPtr, values, count, iL, mtop);

						iL += 32;
					}
					hL += 16;
					if (iL <= hL)
					{
						Estimate16Short(nodesPtr, values, count, iL, mtop);

						iL += 16;
					}
					hL += 8;
					if (iL <= hL)
					{
						Estimate8Short(nodesPtr, values, count, iL, mtop);
					}
#elif defined(OPTION_AVX2)
					int hL = LH + cH - 8;
					int iL = cH;
					while (iL <= hL)
					{
						Estimate16Short(nodesPtr, values, count, iL, mtop);

						iL += 16;
					}
					hL += 8;
					if (iL <= hL)
					{
						Estimate8Short(nodesPtr, values, count, iL, mtop);
					}
#else
					for (int iL = cH, hL = LH + cH; iL <= hL; iL += 8)
					{
						Estimate8Short(nodesPtr, values, count, iL, mtop);
					}
#endif
				}
			}
			else
			{
				const int pH = reverse ? (pbits & 1) : (pbits >> 8);
				const int pL = reverse ? (pbits >> 8) : (pbits & 1);

				HL += pH;

				for (int iH = HL; iH < HH; iH += 2)
				{
					int cH = (iH << shift);

					int hL = LH + cH - 24;
					int iL = cH;
					while (iL <= hL)
					{
						Estimate16PShort(nodesPtr, values, count, iL, mtop, pL);

						iL += 32;
					}
					hL += 16;
					if (iL <= hL)
					{
						Estimate8PShort(nodesPtr, values, count, iL, mtop, pL);

						iL += 16;
					}
					hL += 8;
					if (iL <= hL)
					{
						Estimate4PShort(nodesPtr, values, count, iL, mtop, pL);
					}
				}
			}
		}
		else
		{
#if defined(OPTION_COUNTERS)
			gEstimateFull++;
#endif

			if constexpr (pbits < 0)
			{
				for (int iH = HL; iH < HH; iH++)
				{
					int cH = (iH << shift);

					for (int iL = cH, hL = LH + cH; iL <= hL; iL += 8)
					{
						Estimate8Full(nodesPtr, values, count, iL, mtop);
					}
				}
			}
			else
			{
				const int pH = reverse ? (pbits & 1) : (pbits >> 8);
				const int pL = reverse ? (pbits >> 8) : (pbits & 1);

				HL += pH;

				for (int iH = HL; iH < HH; iH += 2)
				{
					int cH = (iH << shift);

					for (int iL = cH, hL = LH + cH; iL <= hL; iL += 8)
					{
						Estimate4PFull(nodesPtr, values, count, iL, mtop, pL);
					}
				}
			}
		}

		Node* sorted = nodes1;
		if (int nodesCount = int(nodesPtr - nodes1); nodesCount)
		{
			if (QuerySize == 1)
			{
				Node best = sorted[0];
				for (int i = 1; i < nodesCount; i++)
				{
					Node item = sorted[i];
					if (best.Error > item.Error)
					{
						best = item;
					}
				}
				sorted[0] = best;

				MinErr = best.Error * weight;
				Count = 1;
			}
			else
			{
				sorted = radix_sort(nodes1, nodes2, (uint32_t)nodesCount);

				MinErr = sorted[0].Error * weight;
				Count = Min(nodesCount, QuerySize);
			}
		}

		if (reverse)
		{
			for (int i = 0, n = Count; i < n; i++)
			{
				int e = sorted[i].Error;
				int c = sorted[i].Color;

				{
					int cL = (c >> shift) << (8 - shift);
					int cH = (c & ((1 << shift) - 1)) << (8 - shift);

					c = ((cH + (cH >> shift)) << 8) + cL + (cL >> shift);
				}

				Err[i].Error = e * weight;
				Err[i].Color = c;
			}
		}
		else
		{
			for (int i = 0, n = Count; i < n; i++)
			{
				int e = sorted[i].Error;
				int c = sorted[i].Color;

				if constexpr (shift != 8)
				{
					int cH = (c >> shift) << (8 - shift);
					int cL = (c & ((1 << shift) - 1)) << (8 - shift);

					c = ((cH + (cH >> shift)) << 8) + cL + (cL >> shift);
				}

				Err[i].Error = e * weight;
				Err[i].Color = c;
			}
		}
	}

private:
#if defined(OPTION_AVX512)

	static INLINED void Store8N(Node*& nodesPtr, __m128i msum, uint32_t flags, const int c) noexcept
	{
		for (int i = 0; i < 8; i++)
		{
			nodesPtr->Init((int)_mm_cvtsi128_si64(msum) & 0xFFFF, c + i);
			nodesPtr += flags & 1;

			msum = _mm_alignr_epi8(msum, msum, 2);
			flags >>= 1;
		}
	}

#endif

	static INLINED void Store8(Node*& nodesPtr, __m128i msum, int flags, const int c) noexcept
	{
		for (int i = 0; i < 8; i++)
		{
			nodesPtr->Init((int)_mm_cvtsi128_si64(msum) & 0xFFFF, c + i);
			nodesPtr += flags & 1;

			msum = _mm_alignr_epi8(msum, msum, 2);
			flags >>= 2;
		}
	}

	static INLINED void Store4(Node*& nodesPtr, __m128i msum, int flags, const int c) noexcept
	{
		for (int i = 0; i < 4; i++)
		{
			nodesPtr->Init((int)_mm_cvtsi128_si64(msum), c + i);
			nodesPtr += flags & 1;

			msum = _mm_shuffle_epi32(msum, _MM_SHUFFLE(0, 3, 2, 1));
			flags >>= 4;
		}
	}

#if defined(OPTION_AVX512)

	static INLINED void Store4PN(Node*& nodesPtr, __m128i msum, int flags, const int c) noexcept
	{
		for (int i = 0; i < 4; i++)
		{
			nodesPtr->Init((int)_mm_cvtsi128_si64(msum), c + i + i);
			nodesPtr += flags & 1;

			msum = _mm_shuffle_epi32(msum, _MM_SHUFFLE(0, 3, 2, 1));
			flags >>= 1;
		}
	}

#endif

	static INLINED void Store4P(Node*& nodesPtr, __m128i msum, int flags, const int c) noexcept
	{
		for (int i = 0; i < 4; i++)
		{
			nodesPtr->Init((int)_mm_cvtsi128_si64(msum), c + i + i);
			nodesPtr += flags & 1;

			msum = _mm_shuffle_epi32(msum, _MM_SHUFFLE(0, 3, 2, 1));
			flags >>= 4;
		}
	}

#if defined(OPTION_AVX512)

	static INLINED void Estimate32Short(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop) noexcept
	{
		const __m512i wwater = _mm512_maskz_broadcastw_epi16(kFullMask32, mtop);

		__m512i wsum = _mm512_setzero_si512();
		uint32_t flags = ~0ui32;

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m512i* p = (const __m512i*)&value[c];

			__m512i wadd = _mm512_maskz_load_epi64(kFullMask8, p);

			wsum = _mm512_maskz_adds_epu16(kFullMask32, wsum, wadd);

			flags = _mm512_mask_cmp_epu16_mask(kFullMask32, wwater, wsum, _MM_CMPINT_GT);
			if (!flags)
				return;
		}

		Store8N(nodesPtr, _mm512_castsi512_si128(wsum), flags, c);
		Store8N(nodesPtr, _mm512_maskz_extracti32x4_epi32(kFullMask8, wsum, 1), flags >> 8, c + 8);
		Store8N(nodesPtr, _mm512_maskz_extracti32x4_epi32(kFullMask8, wsum, 2), flags >> 16, c + 16);
		Store8N(nodesPtr, _mm512_maskz_extracti32x4_epi32(kFullMask8, wsum, 3), flags >> 24, c + 24);
	}

	static INLINED void Estimate16Short(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop) noexcept
	{
		const __m256i vwater = _mm256_broadcastw_epi16(mtop);

		__m256i vsum = _mm256_setzero_si256();
		int flags = 0xFFFF;

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m256i* p = (const __m256i*)&value[c];

			__m256i vadd = _mm256_load_si256(p);

			vsum = _mm256_adds_epu16(vsum, vadd);

			flags = _mm256_cmp_epu16_mask(vwater, vsum, _MM_CMPINT_GT);
			if (!flags)
				return;
		}

		Store8N(nodesPtr, _mm256_castsi256_si128(vsum), flags, c);
		Store8N(nodesPtr, _mm256_extracti128_si256(vsum, 1), flags >> 8, c + 8);
	}

#elif defined(OPTION_AVX2)

	static INLINED void Estimate16Short(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop) noexcept
	{
		const __m256i vsign = _mm256_set1_epi16(-0x8000);
		const __m256i vwater = _mm256_xor_si256(_mm256_broadcastw_epi16(mtop), vsign);

		__m256i vsum = _mm256_setzero_si256();
		int flags = -1;

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m256i* p = (const __m256i*)&value[c];

			__m256i vadd = _mm256_load_si256(p);

			vsum = _mm256_adds_epu16(vsum, vadd);

			flags = _mm256_movemask_epi8(_mm256_cmpgt_epi16(vwater, _mm256_xor_si256(vsum, vsign)));
			if (!flags)
				return;
		}

		Store8(nodesPtr, _mm256_castsi256_si128(vsum), flags, c);
		Store8(nodesPtr, _mm256_extracti128_si256(vsum, 1), flags >> 16, c + 8);
	}

#endif

	static INLINED void Estimate8Short(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop) noexcept
	{
		const __m128i msign = _mm_set1_epi16(-0x8000);
		const __m128i mwater = _mm_xor_si128(_mm_packus_epi32(mtop, mtop), msign);

		__m128i msum = _mm_setzero_si128();
		int flags = 0xFFFF;

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m128i* p = (const __m128i*)&value[c];

			__m128i madd = _mm_load_si128(p);

			msum = _mm_adds_epu16(msum, madd);

			flags = _mm_movemask_epi8(_mm_cmpgt_epi16(mwater, _mm_xor_si128(msum, msign)));
			if (!flags)
				return;
		}

		Store8(nodesPtr, msum, flags, c);
	}

#if defined(OPTION_AVX512)

	static INLINED void Estimate16PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m512i wsum = _mm512_setzero_si512();

		const __m512i wmask = pL ?
			_mm512_set_epi8(
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2,
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2,
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2,
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm512_set_epi8(
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0,
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0,
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0,
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m512i* p = (const __m512i*)&value[c];

			__m512i wadd = _mm512_maskz_load_epi64(kFullMask8, p);

			wsum = _mm512_maskz_adds_epu16(kFullMask32, wsum, wadd);
		}

		const __m512i wtop = _mm512_maskz_broadcastd_epi32(kFullMask16, mtop);

		wsum = _mm512_maskz_shuffle_epi8(kFullMask64, wsum, wmask);

		int flags = _mm512_mask_cmp_epi32_mask(kFullMask16, wtop, wsum, _MM_CMPINT_GT);

		if (flags & 0xF)
		{
			Store4PN(nodesPtr, _mm512_castsi512_si128(wsum), flags, c + pL);
		}

		if (flags & 0xF0)
		{
			Store4PN(nodesPtr, _mm512_maskz_extracti32x4_epi32(kFullMask8, wsum, 1), flags >> 4, c + pL + 8);
		}

		if (flags & 0xF00)
		{
			Store4PN(nodesPtr, _mm512_maskz_extracti32x4_epi32(kFullMask8, wsum, 2), flags >> 8, c + pL + 16);
		}

		if (flags & 0xF000)
		{
			Store4PN(nodesPtr, _mm512_maskz_extracti32x4_epi32(kFullMask8, wsum, 3), flags >> 12, c + pL + 24);
		}
	}

	static INLINED void Estimate8PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m256i vsum = _mm256_setzero_si256();

		const __m256i vmask = pL ?
			_mm256_set_epi8(
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2,
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm256_set_epi8(
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0,
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m256i* p = (const __m256i*)&value[c];

			__m256i vadd = _mm256_load_si256(p);

			vsum = _mm256_adds_epu16(vsum, vadd);
		}

		__m256i vtop = _mm256_broadcastd_epi32(mtop);

		vsum = _mm256_shuffle_epi8(vsum, vmask);

		int flags = _mm256_movemask_epi8(_mm256_cmpgt_epi32(vtop, vsum));

		if (flags & 0xFFFF)
		{
			Store4P(nodesPtr, _mm256_castsi256_si128(vsum), flags, c + pL);
		}

		if (flags >> 16)
		{
			Store4P(nodesPtr, _mm256_extracti128_si256(vsum, 1), flags >> 16, c + pL + 8);
		}
	}

#elif defined(OPTION_AVX2)

	static INLINED void Estimate16PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m256i vsum0 = _mm256_setzero_si256();
		__m256i vsum1 = _mm256_setzero_si256();

		const __m256i vmask = pL ?
			_mm256_set_epi8(
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2,
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm256_set_epi8(
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0,
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m256i* p = (const __m256i*)&value[c];

			__m256i vadd0 = _mm256_load_si256(&p[0]);
			__m256i vadd1 = _mm256_load_si256(&p[1]);

			vsum0 = _mm256_adds_epu16(vsum0, vadd0);
			vsum1 = _mm256_adds_epu16(vsum1, vadd1);
		}

		__m256i vtop = _mm256_broadcastd_epi32(mtop);

		vsum0 = _mm256_shuffle_epi8(vsum0, vmask);
		vsum1 = _mm256_shuffle_epi8(vsum1, vmask);

		int flags0 = _mm256_movemask_epi8(_mm256_cmpgt_epi32(vtop, vsum0));
		int flags1 = _mm256_movemask_epi8(_mm256_cmpgt_epi32(vtop, vsum1));

		if (flags0 & 0xFFFF)
		{
			Store4P(nodesPtr, _mm256_castsi256_si128(vsum0), flags0, c + pL);
		}

		if (flags0 >> 16)
		{
			Store4P(nodesPtr, _mm256_extracti128_si256(vsum0, 1), flags0 >> 16, c + pL + 8);
		}

		if (flags1 & 0xFFFF)
		{
			Store4P(nodesPtr, _mm256_castsi256_si128(vsum1), flags1, c + pL + 16);
		}

		if (flags1 >> 16)
		{
			Store4P(nodesPtr, _mm256_extracti128_si256(vsum1, 1), flags1 >> 16, c + pL + 24);
		}
	}

	static INLINED void Estimate8PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m256i vsum = _mm256_setzero_si256();

		const __m256i vmask = pL ?
			_mm256_set_epi8(
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2,
				-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm256_set_epi8(
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0,
				-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m256i* p = (const __m256i*)&value[c];

			__m256i vadd = _mm256_load_si256(p);

			vsum = _mm256_adds_epu16(vsum, vadd);
		}

		__m256i vtop = _mm256_broadcastd_epi32(mtop);

		vsum = _mm256_shuffle_epi8(vsum, vmask);

		int flags = _mm256_movemask_epi8(_mm256_cmpgt_epi32(vtop, vsum));

		if (flags & 0xFFFF)
		{
			Store4P(nodesPtr, _mm256_castsi256_si128(vsum), flags, c + pL);
		}

		if (flags >> 16)
		{
			Store4P(nodesPtr, _mm256_extracti128_si256(vsum, 1), flags >> 16, c + pL + 8);
		}
	}

#else

	static INLINED void Estimate16PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m128i msum0 = _mm_setzero_si128();
		__m128i msum1 = _mm_setzero_si128();
		__m128i msum2 = _mm_setzero_si128();
		__m128i msum3 = _mm_setzero_si128();

		const __m128i mmask = pL ?
			_mm_set_epi8(-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm_set_epi8(-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m128i* p = (const __m128i*)&value[c];

			__m128i madd0 = _mm_load_si128(&p[0]);
			__m128i madd1 = _mm_load_si128(&p[1]);
			__m128i madd2 = _mm_load_si128(&p[2]);
			__m128i madd3 = _mm_load_si128(&p[3]);

			msum0 = _mm_adds_epu16(msum0, madd0);
			msum1 = _mm_adds_epu16(msum1, madd1);
			msum2 = _mm_adds_epu16(msum2, madd2);
			msum3 = _mm_adds_epu16(msum3, madd3);
		}

		msum0 = _mm_shuffle_epi8(msum0, mmask);
		msum1 = _mm_shuffle_epi8(msum1, mmask);
		msum2 = _mm_shuffle_epi8(msum2, mmask);
		msum3 = _mm_shuffle_epi8(msum3, mmask);

		int flags0 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum0));
		int flags1 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum1));

		if (flags0)
		{
			Store4P(nodesPtr, msum0, flags0, c + pL);
		}

		if (flags1)
		{
			Store4P(nodesPtr, msum1, flags1, c + pL + 8);
		}

		int flags2 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum2));
		int flags3 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum3));

		if (flags2)
		{
			Store4P(nodesPtr, msum2, flags2, c + pL + 16);
		}

		if (flags3)
		{
			Store4P(nodesPtr, msum3, flags3, c + pL + 24);
		}
	}

	static INLINED void Estimate8PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m128i msum0 = _mm_setzero_si128();
		__m128i msum1 = _mm_setzero_si128();

		const __m128i mmask = pL ?
			_mm_set_epi8(-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm_set_epi8(-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m128i* p = (const __m128i*)&value[c];

			__m128i madd0 = _mm_load_si128(&p[0]);
			__m128i madd1 = _mm_load_si128(&p[1]);

			msum0 = _mm_adds_epu16(msum0, madd0);
			msum1 = _mm_adds_epu16(msum1, madd1);
		}

		msum0 = _mm_shuffle_epi8(msum0, mmask);
		msum1 = _mm_shuffle_epi8(msum1, mmask);

		int flags0 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum0));
		int flags1 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum1));

		if (flags0)
		{
			Store4P(nodesPtr, msum0, flags0, c + pL);
		}

		if (flags1)
		{
			Store4P(nodesPtr, msum1, flags1, c + pL + 8);
		}
	}

#endif

	static INLINED void Estimate4PShort(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m128i msum = _mm_setzero_si128();

		const __m128i mmask = pL ?
			_mm_set_epi8(-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm_set_epi8(-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m128i* p = (const __m128i*)&value[c];

			__m128i madd = _mm_load_si128(p);

			msum = _mm_adds_epu16(msum, madd);
		}

		msum = _mm_shuffle_epi8(msum, mmask);

		int flags = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum));
		if (flags)
		{
			Store4P(nodesPtr, msum, flags, c + pL);
		}
	}

	static INLINED void Estimate8Full(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop) noexcept
	{
		__m128i msum0 = _mm_setzero_si128();
		__m128i msum1 = _mm_setzero_si128();

		__m128i mzero = _mm_setzero_si128();

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m128i* p = (const __m128i*)&value[c];

			__m128i madd = _mm_load_si128(p);

			msum0 = _mm_add_epi32(msum0, _mm_unpacklo_epi16(madd, mzero));
			msum1 = _mm_add_epi32(msum1, _mm_unpackhi_epi16(madd, mzero));
		}

		int flags0 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum0));
		int flags1 = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum1));

		if (flags0)
		{
			Store4(nodesPtr, msum0, flags0, c);
		}

		if (flags1)
		{
			Store4(nodesPtr, msum1, flags1, c + 4);
		}
	}

	static INLINED void Estimate4PFull(Node*& nodesPtr, const uint16_t* values[16], const size_t count, const int c, const __m128i mtop, const int pL) noexcept
	{
		__m128i msum = _mm_setzero_si128();

		const __m128i mmask = pL ?
			_mm_set_epi8(-0x80, -0x80, 15, 14, -0x80, -0x80, 11, 10, -0x80, -0x80, 7, 6, -0x80, -0x80, 3, 2) :
			_mm_set_epi8(-0x80, -0x80, 13, 12, -0x80, -0x80, 9, 8, -0x80, -0x80, 5, 4, -0x80, -0x80, 1, 0);

		for (size_t i = 0; i < count; i++)
		{
			auto value = values[i];

			const __m128i* p = (const __m128i*)&value[c];

			__m128i madd = _mm_load_si128(p);

			msum = _mm_add_epi32(msum, _mm_shuffle_epi8(madd, mmask));
		}

		int flags = _mm_movemask_epi8(_mm_cmpgt_epi32(mtop, msum));
		if (flags)
		{
			Store4P(nodesPtr, msum, flags, c + pL);
		}
	}
};
