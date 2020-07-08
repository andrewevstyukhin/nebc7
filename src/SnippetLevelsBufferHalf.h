#pragma once

#include "pch.h"
#include "Bc7Core.h"

#if defined(OPTION_COUNTERS)
inline std::atomic_int gEstimateHalf;

extern std::atomic_int gLevels[17];
#endif

template<int MaxSize>
class LevelsBufferHalf final
{
public:
	int MinErr = kBlockMaximalAlphaError + kBlockMaximalColorError;
	int Count = 0;

	Node Err[MaxSize];

	static ALWAYS_INLINED int Min(int x, int y) noexcept
	{
		return (x < y) ? x : y;
	}

	static ALWAYS_INLINED int Max(int x, int y) noexcept
	{
		return (x > y) ? x : y;
	}

	ALWAYS_INLINED LevelsBufferHalf() noexcept = default;

	ALWAYS_INLINED void SetZeroError(int color) noexcept
	{
		MinErr = 0;
		Count = 1;

		Err[0].Error = 0;
		Err[0].Color = color;
	}

	template<int bits, int pbits, bool transparent, const uint8_t table[0x100][1 << 1 * (2 * (bits + 1) - 1)]>
	NOTINLINED void ComputeChannelLevelsReduced(const Area& area, const size_t offset, const int weight, const int water) noexcept
	{
		const uint8_t* values[16];

		size_t count;
		if constexpr (transparent)
		{
			count = area.Active;

			if (!count)
			{
				SetZeroError(pbits << (7 - bits));
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
		d <<= kDenoise;
		d += (1 << kDenoise) - 1;

		constexpr int tailmask = (1 << (8 - bits)) - 1;
		constexpr int shift = bits + 1;

		int L = ((const short*)&area.MinMax_U16)[offset + offset + 0];
		int H = ((const short*)&area.MinMax_U16)[offset + offset + 1];

		const bool reverse = (L != ((const short*)&area.Bounds_U16)[offset + offset + 0]);

		if constexpr (((pbits >> 8) == (pbits & 1)))
		{
			if (L == H)
			{
				if ((((L >> shift) ^ L) & ((1 << (8 - shift)) - 1)) == 0)
				{
					if (((L >> (8 - shift)) & 1) == (pbits & 1))
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

		top = Min(top, (0xF * 0xF) * 16 + 1);

#if defined(OPTION_COUNTERS)
		gEstimateHalf++;
#endif

		NodeShort nodes1[(1 << shift) * (1 << shift)];
		NodeShort nodes2[(1 << shift) * (1 << shift)];
		NodeShort* nodesPtr = nodes1;

		const __m128i mtop = _mm_shuffle_epi32(_mm_shufflelo_epi16(_mm_cvtsi32_si128(top), 0), 0);

		{
			const int pH = reverse ? (pbits & 1) : (pbits >> 8);
			const int pL = reverse ? (pbits >> 8) : (pbits & 1);

			HL += pH;

			for (int iH = HL; iH < HH; iH += 2)
			{
				int cH = (iH << shift);

				Estimate16PStep8Short(nodesPtr, values, count, cH, LH + cH, mtop, pL);
			}
		}

		NodeShort* sorted = nodes1;
		if (int nodesCount = int(nodesPtr - nodes1); nodesCount)
		{
			sorted = radix_sort(nodes1, nodes2, (uint32_t)nodesCount);

			MinErr = (sorted[0].ColorError >> 16) * weight;
			Count = Min(nodesCount, MaxSize);
		}

		if (reverse)
		{
			for (int i = 0, n = Count; i < n; i++)
			{
				uint32_t ce = sorted[i].ColorError;
				int e = ce >> 16;
				int c = ce & 0xFFFFu;

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
				uint32_t ce = sorted[i].ColorError;
				int e = ce >> 16;
				int c = ce & 0xFFFFu;

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

	static ALWAYS_INLINED void Store8P(NodeShort*& nodesPtr, __m128i msum, int flags, const int c) noexcept
	{
		for (int i = 0; i < 8; i++)
		{
			nodesPtr->Init(_mm_cvtsi128_si64(msum), c + i + i);
			nodesPtr += flags & 1;

			msum = _mm_alignr_epi8(msum, msum, 2);
			flags >>= 2;
		}
	}

#if defined(OPTION_AVX2)

	static INLINED void Estimate16PStep8Short(NodeShort*& nodesPtr, const uint8_t* values[16], const size_t count, const int cH, int hL, const __m128i mtop, const int pL) noexcept
	{
		const __m256i vtop = _mm256_broadcastw_epi16(mtop);

		const __m128i mshift = _mm_cvtsi32_si128(pL << 2);
		const __m128i mmask = _mm_set1_epi8(0xF);

		int c = cH;
		do
		{
			__m256i vsum = _mm256_setzero_si256();

			for (size_t i = 0; i < count; i++)
			{
				auto value = values[i];

				const __m128i* p = (const __m128i*)&value[c >> 1];

				__m128i mdelta = _mm_load_si128(p);
				mdelta = _mm_and_si128(_mm_srl_epi16(mdelta, mshift), mmask);

				__m256i vadd = _mm256_cvtepu8_epi16(mdelta);

				vadd = _mm256_mullo_epi16(vadd, vadd);

				vsum = _mm256_add_epi16(vsum, vadd);
			}

			int flags = _mm256_movemask_epi8(_mm256_cmpgt_epi16(vtop, vsum));

			flags &= static_cast<int>(~0u >> Max(0, (c | 31) - (hL | 7)));

			if (flags & 0xFFFF)
			{
				Store8P(nodesPtr, _mm256_castsi256_si128(vsum), flags, c + pL);
			}

			if (flags & 0xFFFF0000)
			{
				Store8P(nodesPtr, _mm256_extracti128_si256(vsum, 1), flags >> 16, c + pL + 16);
			}

			c += 32;

		} while (c <= hL);
	}

#else

	static INLINED void Estimate16PStep8Short(NodeShort*& nodesPtr, const uint8_t* values[16], const size_t count, const int cH, int hL, const __m128i mtop, const int pL) noexcept
	{
		const __m128i mshift = _mm_cvtsi32_si128(pL << 2);
		const __m128i mmask = _mm_set1_epi8(0xF);

		int c = cH;
		do
		{
			__m128i msum = _mm_setzero_si128();

			for (size_t i = 0; i < count; i++)
			{
				auto value = values[i];

				const __m128i* p = (const __m128i*)&value[c >> 1];

				__m128i mdelta = _mm_loadl_epi64(p);
				mdelta = _mm_and_si128(_mm_srl_epi16(mdelta, mshift), mmask);

				__m128i madd = _mm_cvtepu8_epi16(mdelta);

				madd = _mm_mullo_epi16(madd, madd);

				msum = _mm_add_epi16(msum, madd);
			}

			int flags = _mm_movemask_epi8(_mm_cmpgt_epi16(mtop, msum));

			flags &= (0xFFFF >> Max(0, (c | 15) - (hL | 7)));

			if (flags)
			{
				Store8P(nodesPtr, msum, flags, c + pL);
			}

			c += 16;

		} while (c <= hL);
	}

#endif
};
