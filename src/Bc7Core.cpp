
#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"
#include "Bc7Pca.h"
#include "Metrics.h"

#if defined(OPTION_COUNTERS)
#include "SnippetLevelsBuffer.h"
#endif

#if defined(OPTION_COUNTERS)
static std::atomic_int gCounterModes[8];
static std::atomic_int gCompressAlready, gCompress, gCompressBad;
#endif

bool gDoDraft = false;
bool gDoFast = false;
bool gDoNormal = false;
bool gDoSlow = false;

static INLINED int ComputeOpaqueAlphaError(const Area& area) noexcept
{
	int error = 0;

	if (!area.IsOpaque)
	{
		for (size_t i = 0, n = area.Count; i < n; i++)
		{
			int da = *(const short*)&area.DataMask_I16[i] - 255;
			error += da * da;
		}

		error *= kAlpha;
	}

	return error;
}

static INLINED void ComputeColorCovariances(const Area& area,
	int64_t& covGR, int64_t& covRB, int64_t& covBG) noexcept
{
	__m128i msum = _mm_setzero_si128();
	__m128i msum2 = _mm_setzero_si128();

	size_t count = area.Count;

	for (size_t i = 0; i < count; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_cvtepu16_epi32(mpacked);

		msum = _mm_add_epi16(msum, mpixel);
		msum2 = _mm_add_epi32(msum2, _mm_mullo_epi16(mpixel, _mm_shuffle_epi32(mpixel, _MM_SHUFFLE(1, 3, 2, 0))));
	}

	__m128i mcount = _mm_shuffle_epi32(_mm_cvtsi64_si128(static_cast<int64_t>(count)), 0);

	msum2 = _mm_sub_epi32(_mm_mullo_epi32(msum2, mcount), _mm_mullo_epi32(msum, _mm_shuffle_epi32(msum, _MM_SHUFFLE(1, 3, 2, 0))));

	covGR = _mm_extract_epi32(msum2, 1);
	covRB = _mm_extract_epi32(msum2, 2);
	covBG = _mm_extract_epi32(msum2, 3);
}

static INLINED void ComputeAlphaColorCovariances(const Area& area,
	int64_t& covAG, int64_t& covAR, int64_t& covAB, int64_t& covGR, int64_t& covRB, int64_t& covBG) noexcept
{
	__m128i msum = _mm_setzero_si128();
	__m128i msumA = _mm_setzero_si128();
	__m128i msum2 = _mm_setzero_si128();

	size_t count = area.Active;

	for (size_t i = 0; i < count; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_cvtepu16_epi32(mpacked);

		msum = _mm_add_epi16(msum, mpixel);
		msumA = _mm_add_epi32(msumA, _mm_mullo_epi16(mpixel, _mm_shuffle_epi32(mpixel, 0)));
		msum2 = _mm_add_epi32(msum2, _mm_mullo_epi16(mpixel, _mm_shuffle_epi32(mpixel, _MM_SHUFFLE(1, 3, 2, 0))));
	}

	__m128i mcount = _mm_shuffle_epi32(_mm_cvtsi64_si128(static_cast<int64_t>(count)), 0);

	msumA = _mm_sub_epi32(_mm_mullo_epi32(msumA, mcount), _mm_mullo_epi32(msum, _mm_shuffle_epi32(msum, 0)));
	msum2 = _mm_sub_epi32(_mm_mullo_epi32(msum2, mcount), _mm_mullo_epi32(msum, _mm_shuffle_epi32(msum, _MM_SHUFFLE(1, 3, 2, 0))));

	covAG = _mm_extract_epi32(msumA, 1);
	covAR = _mm_extract_epi32(msumA, 2);
	covAB = _mm_extract_epi32(msumA, 3);

	covGR = _mm_extract_epi32(msum2, 1);
	covRB = _mm_extract_epi32(msum2, 2);
	covBG = _mm_extract_epi32(msum2, 3);
}

static INLINED void MakeCell(Cell& input) noexcept
{
	input.BestColor0 = _mm_setzero_si128();
	input.BestColor1 = _mm_setzero_si128();
	input.BestColor2 = _mm_setzero_si128();
	//input.BestParameter = 0;
	//input.BestMode = 8;

	for (size_t i = 0; i < 16; i++)
	{
		__m128i mc = _mm_cvtepu8_epi16(_mm_cvtsi32_si128(((int*)input.ImageRows_U8)[i]));
		__m128i mmask = _mm_cvtepi8_epi16(_mm_cvtsi32_si128(((int*)input.MaskRows_S8)[i]));

		mc = _mm_and_si128(mc, mmask);

		input.DataMask_I16[i] = _mm_unpacklo_epi64(mc, mmask);
	}

	for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
	{
		input.LazyArea12[partitionIndex] = true;
	}
	for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
	{
		input.LazyArea22[partitionIndex] = true;
	}

	for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
	{
		input.LazyArea13[partitionIndex] = true;
	}
	for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
	{
		input.LazyArea23[partitionIndex] = true;
	}
	for (size_t partitionIndex = 0; partitionIndex < 64; partitionIndex++)
	{
		input.LazyArea33[partitionIndex] = true;
	}
}

static NOTINLINED void MakeAreaFromCell(Area& area, const Cell& cell, const size_t count, uint64_t indices) noexcept
{
	// Initialize Indices

	area.Count = static_cast<uint32_t>(count);

	area.ZeroIndex = static_cast<uint8_t>(indices & 0xF);

	{
		uint8_t TransparentIndices[16];

		uint8_t* pOpaque = area.Indices;
		uint8_t* pTransparent = TransparentIndices;

		for (size_t i = 0; i < count; i++)
		{
			const size_t index = indices & 0xF;
			indices >>= 4;

			size_t opaque = ((const uint16_t*)&cell.DataMask_I16[index])[4 + 1] & 1;

			*(opaque ? pOpaque : pTransparent) = static_cast<uint8_t>(index);

			pTransparent++;
			pOpaque += opaque;
			pTransparent -= opaque;
		}

		area.Active = static_cast<uint32_t>(pOpaque - area.Indices);

		uint8_t* p = TransparentIndices;
		while (p != pTransparent)
		{
			*pOpaque++ = *p++;
		}
	}

	// Min & Max

	__m128i m0 = _mm_set1_epi16(255);
	__m128i m1 = _mm_setzero_si128();

	for (size_t i = 0; i < count; i++)
	{
		size_t index = area.Indices[i];

		__m128i mpacked = _mm_load_si128(&cell.DataMask_I16[index]);
		__m128i mc = _mm_unpacklo_epi64(mpacked, mpacked);
		__m128i mmask = _mm_unpackhi_epi64(mpacked, mpacked);

		m0 = _mm_blendv_epi8(m0, _mm_min_epi16(mc, m0), mmask);
		m1 = _mm_max_epi16(m1, mc);

		area.DataMask_I16[i] = mpacked;
	}

	m0 = _mm_andnot_si128(_mm_cmplt_epi16(m1, m0), m0);

	__m128i mbounds = _mm_unpacklo_epi16(m0, m1);

	_mm_store_si128(&area.MinMax_U16, mbounds);

	// Flags

	area.IsOpaque = ((_mm_extract_epi16(mbounds, 0) & _mm_extract_epi16(mbounds, 1)) == 255);

	area.BestPca3 = -1;

	// Orient channels

	if (area.IsOpaque)
	{
		int64_t covGR = 0, covRB = 0, covBG = 0;
		ComputeColorCovariances(area, covGR, covRB, covBG);

		for (;;)
		{
			bool changes = false;

			int64_t b = covBG * kGreen + covRB * kRed;
			if (b < 0)
			{
				mbounds = _mm_shufflehi_epi16(mbounds, _MM_SHUFFLE(2, 3, 1, 0));

				covBG = -covBG;
				covRB = -covRB;

				changes = true;
			}

			int64_t r = covGR * kGreen + covRB * kBlue;
			if (r < 0)
			{
				mbounds = _mm_shufflehi_epi16(mbounds, _MM_SHUFFLE(3, 2, 0, 1));

				covGR = -covGR;
				covRB = -covRB;

				changes = true;
			}

			if (!changes)
				break;
		}
	}
	else
	{
		int64_t covAG = 0, covAR = 0, covAB = 0, covGR = 0, covRB = 0, covBG = 0;
		ComputeAlphaColorCovariances(area, covAG, covAR, covAB, covGR, covRB, covBG);

		for (;;)
		{
			bool changes = false;

			int64_t b = covAB * kAlpha + covBG * kGreen + covRB * kRed;
			if (b < 0)
			{
				mbounds = _mm_shufflehi_epi16(mbounds, _MM_SHUFFLE(2, 3, 1, 0));

				covAB = -covAB;
				covBG = -covBG;
				covRB = -covRB;

				changes = true;
			}

			int64_t r = covAR * kAlpha + covGR * kGreen + covRB * kBlue;
			if (r < 0)
			{
				mbounds = _mm_shufflehi_epi16(mbounds, _MM_SHUFFLE(3, 2, 0, 1));

				covAR = -covAR;
				covGR = -covGR;
				covRB = -covRB;

				changes = true;
			}

			int64_t g = covAG * kAlpha + covGR * kRed + covBG * kBlue;
			if (g < 0)
			{
				mbounds = _mm_shufflelo_epi16(mbounds, _MM_SHUFFLE(2, 3, 1, 0));

				covAG = -covAG;
				covGR = -covGR;
				covBG = -covBG;

				changes = true;
			}

			if (!changes)
				break;
		}
	}

	_mm_store_si128(&area.Bounds_U16, mbounds);
}

Area& GetArea(Area& area, bool& lazy, const Cell& cell, const uint64_t indices) noexcept
{
	if (lazy)
	{
		lazy = false;

		MakeAreaFromCell(area, cell, indices & 0xF, indices >> 4);
	}

	return area;
}

int AreaGetBestPca3(Area& area) noexcept
{
#if defined(OPTION_PCA)
	int best = area.BestPca3;
	if (best < 0)
	{
		best = PrincipalComponentAnalysis3(area);
		area.BestPca3 = best;
	}
#else
	(void)area;
	int best = 0;
#endif

	return best;
}

void AreaReduceTable2(const Area& area, __m128i& mc, uint64_t& indices) noexcept
{
	const int zerobit = area.ZeroIndex * 2 + 1;

	if (indices & (1ui64 << zerobit))
	{
		uint64_t mask = 0;

		for (size_t i = 0, n = area.Count; i < n; i++)
		{
			mask |= 3ui64 << (area.Indices[i] << 1);
		}

		indices ^= mask;

		mc = _mm_shufflelo_epi16(mc, _MM_SHUFFLE(2, 3, 0, 1));
		mc = _mm_shufflehi_epi16(mc, _MM_SHUFFLE(2, 3, 0, 1));
	}
}

void AreaReduceTable3(const Area& area, __m128i& mc, uint64_t& indices) noexcept
{
	const int zerobit = area.ZeroIndex * 3 + 2;

	if (indices & (1ui64 << zerobit))
	{
		uint64_t mask = 0;

		for (size_t i = 0, n = area.Count; i < n; i++)
		{
			mask |= 7ui64 << (area.Indices[i] * 3);
		}

		indices ^= mask;

		mc = _mm_shufflelo_epi16(mc, _MM_SHUFFLE(2, 3, 0, 1));
		mc = _mm_shufflehi_epi16(mc, _MM_SHUFFLE(2, 3, 0, 1));
	}
}

void AreaReduceTable4(__m128i& mc, uint64_t& indices) noexcept
{
	const int zerobit = 0 * 4 + 3;

	if (indices & (1ui64 << zerobit))
	{
		indices = ~indices;

		mc = _mm_shufflelo_epi16(mc, _MM_SHUFFLE(2, 3, 0, 1));
		mc = _mm_shufflehi_epi16(mc, _MM_SHUFFLE(2, 3, 0, 1));
	}
}

NOTINLINED Node* radix_sort(Node* input, Node* work, size_t N) noexcept
{
	constexpr size_t radix = 7;
	constexpr size_t bucketCount = 1 << radix;
	constexpr size_t bucketMask = bucketCount - 1;

	alignas(32) uint32_t counts[bucketCount];

	Node* A = input;
	Node* B = work;

	int any_error = 0;
	for (size_t i = 0; i < N; i++)
	{
		any_error |= A[i].Error;
	}

	for (size_t shift = 0; shift < 32; shift += radix)
	{
		if (((any_error >> shift) & bucketMask) == 0)
			continue;

		for (size_t i = 0; i < bucketCount; i++)
		{
			counts[i] = 0;
		}

		for (size_t i = 0; i < N; i++)
		{
			size_t value = (uint32_t)A[i].Error;
			counts[(value >> shift) & bucketMask]++;
		}

		uint32_t total = 0;
		for (size_t i = 0; i < bucketCount; i++)
		{
			uint32_t oldCount = counts[i];
			counts[i] = total;
			total += oldCount;
		}

		for (size_t i = 0; i < N; i++)
		{
			Node val = A[i];
			uint32_t p = counts[((size_t)(uint32_t)val.Error >> shift) & bucketMask]++;
			B[p] = val;
		}

		Node* C = A; A = B; B = C;
	}

	return A;
}

NOTINLINED int ComputeSubsetTable(const Area& area, const __m128i mweights, const __m128i mfix, Modulations& state, const int M) noexcept
{
	int good = (1 << M) - 1;
	{
		const int m = M >> 1;

		for (int i = 1; i < m; i++)
		{
			good ^= int(state.Values_I16[i - 1] == state.Values_I16[i]) << i;
		}

		for (int i = M - 2; i >= m; i--)
		{
			good ^= int(state.Values_I16[i + 1] == state.Values_I16[i]) << i;
		}

		good &= ~(int(state.Values_I16[m - 1] == state.Values_I16[m]) << m);

		good ^= int(state.Values_I16[0] == state.Values_I16[M - 1]) << (M - 1);
	}

	int errorBlock = 0;

	const __m128i msign = _mm_set1_epi16(-0x8000);

	for (size_t i = 0, n = area.Active; i < n; i++)
	{
		__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
		__m128i mpixel = _mm_unpacklo_epi64(mpacked, mpacked);
		__m128i mmask = _mm_unpackhi_epi64(mpacked, mpacked);

		int errors[16];

		int bottom = kBlockMaximalAlphaError + kBlockMaximalColorError;

		for (int j = 0; j < M; j++)
		{
			__m128i mx = _mm_loadl_epi64((const __m128i*)&state.Values_I16[j]);

			mx = _mm_sub_epi16(mx, mpixel);

			mx = _mm_mullo_epi16(mx, mx);

			mx = _mm_and_si128(mx, mmask);

			mx = _mm_xor_si128(mx, msign);

			mx = _mm_madd_epi16(mx, mweights);

			mx = _mm_add_epi32(mx, _mm_shuffle_epi32(mx, _MM_SHUFFLE(2, 3, 0, 1)));

			mx = _mm_add_epi32(mx, mfix);

			int error = (int)_mm_cvtsi128_si64(mx);

			errors[j] = error;

			if (bottom > error)
			{
				bottom = error;
			}
		}

		errorBlock += bottom;

		int way = 0;

		for (int j = 0; j < M; j++)
		{
			way |= int(errors[j] == bottom) << j;
		}

		way &= good;

		way |= (1 << M);

		state.Ways[i] = (way & (way - 1) & good) ? way : 0;

		int k = 0;
		while ((way & (1 << k)) == 0) k++;
		state.Loops[i] = k;
	}

	//

	if (area.Active < area.Count)
	{
		good = (1 << M) - 1;
		{
			const int m = M >> 1;

			for (int i = 1; i < m; i++)
			{
				good ^= int(*(const uint16_t*)&state.Values_I16[i - 1] == *(const uint16_t*)&state.Values_I16[i]) << i;
			}

			for (int i = M - 2; i >= m; i--)
			{
				good ^= int(*(const uint16_t*)&state.Values_I16[i + 1] == *(const uint16_t*)&state.Values_I16[i]) << i;
			}

			good &= ~(int(*(const uint16_t*)&state.Values_I16[m - 1] == *(const uint16_t*)&state.Values_I16[m]) << m);

			good ^= int(*(const uint16_t*)&state.Values_I16[0] == *(const uint16_t*)&state.Values_I16[M - 1]) << (M - 1);
		}

		int errors[16];

		int bottom = kBlockMaximalAlphaError + kBlockMaximalColorError;

		for (int j = 0; j < M; j++)
		{
			int x = *(const uint16_t*)&state.Values_I16[j];

			int error = x * x;

			errors[j] = error;

			if (bottom > error)
			{
				bottom = error;
			}
		}

		errorBlock += bottom * _mm_extract_epi16(mweights, 0) * int(area.Count - area.Active);

		int way = 0;

		for (int j = 0; j < M; j++)
		{
			way |= int(errors[j] == bottom) << j;
		}

		way &= good;

		way |= (1 << M);

		int k = 0;
		while ((way & (1 << k)) == 0) k++;

		way = (way & (way - 1) & good) ? way : 0;

		for (size_t i = area.Active, n = area.Count; i < n; i++)
		{
			state.Ways[i] = way;

			state.Loops[i] = k;
		}
	}

	//

	double best = -(kAlpha + kColor + 0.1);

	for (;;)
	{
		SSIM_INIT();

		for (size_t i = 0, n = area.Count; i < n; i++)
		{
			__m128i mpacked = _mm_load_si128(&area.DataMask_I16[i]);
			__m128i mb = _mm_cvtepu16_epi32(mpacked);
			__m128i mmask = _mm_unpackhi_epi64(mpacked, mpacked);

			__m128i mvalue = _mm_loadl_epi64((const __m128i*)&state.Values_I16[state.Loops[i]]);
			__m128i mt = _mm_cvtepu16_epi32(_mm_and_si128(mvalue, mmask));

			SSIM_UPDATE(mt, mb);
		}

		// SSIM_CLOSE(PowerOfTwo(area.Count))
		const int count = static_cast<int>(area.Count);
		if (count > 8)
		{
			sab = _mm_slli_epi32(sab, 4 + 1);
			saa_sbb = _mm_slli_epi32(saa_sbb, 4);
		}
		else if (count > 4)
		{
			sab = _mm_slli_epi32(sab, 3 + 1);
			saa_sbb = _mm_slli_epi32(saa_sbb, 3);
		}
		else if (count > 2)
		{
			sab = _mm_slli_epi32(sab, 2 + 1);
			saa_sbb = _mm_slli_epi32(saa_sbb, 2);
		}
		else if (count > 1)
		{
			sab = _mm_slli_epi32(sab, 1 + 1);
			saa_sbb = _mm_slli_epi32(saa_sbb, 1);
		}
		else
		{
			sab = _mm_slli_epi32(sab, 0 + 1);
		}
		__m128i sasb = _mm_mullo_epi32(sa, sb);
		sasb = _mm_add_epi32(sasb, sasb);
		__m128i sasa_sbsb = _mm_add_epi32(_mm_mullo_epi32(sa, sa), _mm_mullo_epi32(sb, sb));
		sab = _mm_sub_epi32(sab, sasb);
		saa_sbb = _mm_sub_epi32(saa_sbb, sasa_sbsb);

		SSIM_FINAL(mssim_ga, gSsim16k1L, gSsim16k2L);
		SSIM_OTHER();
		SSIM_FINAL(mssim_br, gSsim16k1L, gSsim16k2L);

		double ssim = _mm_cvtsd_f64(mssim_ga) * _mm_extract_epi16(mweights, 0);
		ssim += _mm_cvtsd_f64(_mm_unpackhi_pd(mssim_ga, mssim_ga)) * _mm_extract_epi16(mweights, 1);
		ssim += _mm_cvtsd_f64(mssim_br) * _mm_extract_epi16(mweights, 2);
		ssim += _mm_cvtsd_f64(_mm_unpackhi_pd(mssim_br, mssim_br)) * _mm_extract_epi16(mweights, 3);

		if (best < ssim)
		{
			best = ssim;

			memcpy(state.Best, state.Loops, sizeof(state.Best));
		}

		int i = 0;
		for (;; )
		{
			if (int way = state.Ways[i])
			{
				int k = state.Loops[i];
				do { k++; } while ((way & (1 << k)) == 0);
				if (k < M)
				{
					state.Loops[i] = k;
					break;
				}

				k = 0;
				while ((way & (1 << k)) == 0) k++;
				state.Loops[i] = k;
			}

			if (++i >= count)
				return errorBlock;
		}
	}
}

void DecompressBlock(uint8_t input[16], Cell& output) noexcept
{
	const int encoded_mode = input[0];
	if (encoded_mode & 0xF)
	{
		if (encoded_mode & 3)
		{
			if (encoded_mode & 1)
			{
				Mode0::DecompressBlock(input, output);
			}
			else
			{
				Mode1::DecompressBlock(input, output);
			}
		}
		else
		{
			if (encoded_mode & 4)
			{
				Mode2::DecompressBlock(input, output);
			}
			else
			{
				Mode3::DecompressBlock(input, output);
			}
		}
	}
	else if (encoded_mode & 0xF0)
	{
		if (encoded_mode & 0x30)
		{
			if (encoded_mode & 0x10)
			{
				Mode4::DecompressBlock(input, output);
			}
			else
			{
				Mode5::DecompressBlock(input, output);
			}
		}
		else
		{
			if (encoded_mode & 0x40)
			{
				Mode6::DecompressBlock(input, output);
			}
			else
			{
				Mode7::DecompressBlock(input, output);
			}
		}
	}
	else
	{
		__m128i mzero = _mm_setzero_si128();

		_mm_store_si128(&output.ImageRows_U8[0], mzero);
		_mm_store_si128(&output.ImageRows_U8[1], mzero);
		_mm_store_si128(&output.ImageRows_U8[2], mzero);
		_mm_store_si128(&output.ImageRows_U8[3], mzero);

		output.BestMode = 8;
	}
}

static void CompressBlockFastTransparent(Cell& input) noexcept
{
	if (input.Error.Total > 0)
	{
		Mode4::CompressBlockFast(input);

		if (input.Error.Total > 0)
		{
			Mode6::CompressBlockFast(input);

			if (input.Error.Total > 0)
			{
				Mode5::CompressBlockFast(input);

				if (input.Error.Total > 0)
				{
					Mode7::CompressBlockFast(input);
				}
			}
		}
	}
}

static void CompressBlockFastOpaque(Cell& input) noexcept
{
	// DetectGlitches
	if (*(const short*)&input.Area1.MinMax_U16 <= (255 - 16))
		return;

	input.OpaqueAlphaError = ComputeOpaqueAlphaError(input.Area1);

	if (input.Error.Total > input.OpaqueAlphaError)
	{
		Mode1::CompressBlockFast(input);

		if (input.Error.Total > input.OpaqueAlphaError)
		{
			Mode3::CompressBlockFast(input);

			if (input.Error.Total > input.OpaqueAlphaError)
			{
				Mode2::CompressBlockFast(input);

				if (input.Error.Total > input.OpaqueAlphaError)
				{
					Mode0::CompressBlockFast(input);
				}
			}
		}
	}
}

void CompressBlock(uint8_t output[16], Cell& input) noexcept
{
	if (!output[0])
	{
		*(uint64_t*)&output[0] = 1 << 6;
		*(uint64_t*)&output[8] = 0;
	}

#if defined(OPTION_SELFCHECK)
	uint8_t saved_output[16];
	memcpy(saved_output, output, sizeof(saved_output));
#endif

	Cell temp;
	DecompressBlock(output, temp);
	input.BestParameter = temp.BestParameter;
	input.BestMode = temp.BestMode;

	input.Error = CompareBlocks(input, temp);
	if (input.Error.Total > 0)
	{
		MakeCell(input);

		MakeAreaFromCell(input.Area1, input, 16, gTableSelection11);

		if (gDoDraft)
		{
			int water = input.Error.Total;

			if (input.Area1.IsOpaque)
			{
				CompressBlockFastOpaque(input);
				CompressBlockFastTransparent(input);
			}
			else
			{
				CompressBlockFastTransparent(input);
				CompressBlockFastOpaque(input);
			}

#if defined(OPTION_COUNTERS)
			gCompress++;
#endif

			if (gDoFast)
			{
				input.PersonalParameter = input.BestParameter;
				input.PersonalMode = input.BestMode;

				switch (input.BestMode)
				{
				case 0:
					if (input.Error.Total > input.OpaqueAlphaError)
					{
						Mode0::CompressBlock(input);
					}
					break;

				case 1:
					if (input.Error.Total > input.OpaqueAlphaError)
					{
						Mode1::CompressBlock(input);
					}
					break;

				case 2:
					if (input.Error.Total > input.OpaqueAlphaError)
					{
						Mode2::CompressBlock(input);
					}
					break;

				case 3:
					if (input.Error.Total > input.OpaqueAlphaError)
					{
						Mode3::CompressBlock(input);
					}
					break;

				case 4:
					if (input.Error.Total > 0)
					{
						Mode4::CompressBlock(input);
					}
					break;

				case 5:
					if (input.Error.Total > 0)
					{
						Mode5::CompressBlock(input);
					}
					break;

				case 6:
					if (input.Error.Total > 0)
					{
						Mode6::CompressBlock(input);
					}
					break;

				case 7:
					if (input.Error.Total > 0)
					{
						Mode7::CompressBlock(input);
					}
					break;
				}

				if (input.Area1.IsOpaque)
				{
					if (input.Error.Total > input.OpaqueAlphaError)
					{
						Mode2::CompressBlockFull(input);

						if (input.Error.Total > input.OpaqueAlphaError)
						{
							Mode0::CompressBlockFull(input);

							if (gDoSlow)
							{
								if (input.Error.Total > input.OpaqueAlphaError)
								{
									Mode1::CompressBlockFull(input);

									if (input.Error.Total > input.OpaqueAlphaError)
									{
										Mode3::CompressBlockFull(input);
									}
								}
							}
						}
					}
				}
				else
				{
					if (input.Error.Total > 0)
					{
						Mode5::CompressBlockFull(input);

						if (input.Error.Total > 0)
						{
							Mode4::CompressBlockFull(input);

							if (gDoNormal)
							{
								if (input.Error.Total > 0)
								{
									Mode7::CompressBlockFull(input);

									if (gDoSlow)
									{
										if (input.Error.Total > 0)
										{
											Mode6::CompressBlockFull(input);
										}
									}
								}
							}
						}
					}
				}
			}

			if (water > input.Error.Total)
			{
				switch (input.BestMode)
				{
				case 0:
					Mode0::FinalPackBlock(output, input);
					break;

				case 1:
					Mode1::FinalPackBlock(output, input);
					break;

				case 2:
					Mode2::FinalPackBlock(output, input);
					break;

				case 3:
					Mode3::FinalPackBlock(output, input);
					break;

				case 4:
					Mode4::FinalPackBlock(output, input);
					break;

				case 5:
					Mode5::FinalPackBlock(output, input);
					break;

				case 6:
					Mode6::FinalPackBlock(output, input);
					break;

				case 7:
					Mode7::FinalPackBlock(output, input);
					break;
				}

				DecompressBlock(output, temp);
			}
		}

#if defined(OPTION_COUNTERS)
		if (DetectGlitches(input, temp))
		{
			gCompressBad++;
		}
#endif

#if defined(OPTION_SELFCHECK)
		DecompressBlock(output, temp);
		auto e = CompareBlocks(input, temp);
		if ((e.Alpha != input.Error.Alpha) || (e.Total != input.Error.Total))
		{
			__debugbreak();
			memcpy(output, saved_output, sizeof(saved_output));
		}
#endif

		input.Quality = CompareBlocksSSIM(input, temp);
	}
	else
	{
#if defined(OPTION_COUNTERS)
		gCompressAlready++;
#endif

		input.Quality = BlockSSIM(1.0, 1.0);
	}

#if defined(OPTION_COUNTERS)
	gCounterModes[input.BestMode]++;
#endif
}

void CompressStatistics()
{
#if defined(OPTION_COUNTERS)
	PRINTF("[Transparent]\tM4 = %i, M6 = %i, M5 = %i, M7 = %i",
		gCounterModes[4].load(), gCounterModes[6].load(), gCounterModes[5].load(), gCounterModes[7].load());

	PRINTF("[Opaque]\tM1 = %i, M3 = %i, M2 = %i, M0 = %i",
		gCounterModes[1].load(), gCounterModes[3].load(), gCounterModes[2].load(), gCounterModes[0].load());

	PRINTF("[Compress]\tAlready = %i, Compress = %i, Bad = %i",
		gCompressAlready.load(), gCompress.load(), gCompressBad.load());

	Mode1::PrintCounters();
	Mode3::PrintCounters();
	Mode2::PrintCounters();
	Mode0::PrintCounters();

	Mode4::PrintCounters();
	Mode6::PrintCounters();
	Mode5::PrintCounters();
	Mode7::PrintCounters();

	PRINTF("[Estimate]\tFull = %i, Short = %i",
		gEstimateFull.load(), gEstimateShort.load());

	PRINTF("\t\t[1] = %i, [2] = %i, [3] = %i, [4] = %i",
		gLevels[1].load(), gLevels[2].load(), gLevels[3].load(), gLevels[4].load());
	PRINTF("\t\t[5] = %i, [6] = %i, [7] = %i, [8] = %i",
		gLevels[5].load(), gLevels[6].load(), gLevels[7].load(), gLevels[8].load());
	PRINTF("\t\t[9] = %i, [10] = %i, [11] = %i, [12] = %i",
		gLevels[9].load(), gLevels[10].load(), gLevels[11].load(), gLevels[12].load());
	PRINTF("\t\t[13] = %i, [14] = %i, [15] = %i, [16] = %i",
		gLevels[13].load(), gLevels[14].load(), gLevels[15].load(), gLevels[16].load());
#endif
}
