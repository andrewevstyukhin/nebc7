#pragma once

#include "pch.h"

constexpr int kBlockMaximalAlphaError = 16 * (255 * 255) * kAlpha + 1;
constexpr int kBlockMaximalColorError = 16 * (255 * 255) * kColor + 1;

struct BlockError
{
	int Alpha, Total;

	BlockError(int alpha, int total) noexcept
		: Alpha(alpha)
		, Total(total)
	{
	}

	void Add(const BlockError& other) noexcept
	{
		Alpha += other.Alpha;
		Total += other.Total;
	}
};

struct BlockSSIM
{
	double Alpha, Color;

	BlockSSIM(double alpha, double color) noexcept
		: Alpha(alpha)
		, Color(color)
	{
	}
};

// A,G,R,B
struct alignas(64) Area
{
	__m128i DataMask_I16[16];

	uint8_t Indices[16];
	uint32_t Count, Active;

	bool IsOpaque;

	uint8_t ZeroIndex;

	int BestPca3;

	__m128i MinMax_U16;
	__m128i Bounds_U16;
};

// A,G,R,B
struct alignas(64) Cell
{
	__m128i BestColor0, BestColor1, BestColor2;
	uint64_t BestParameter;
	uint32_t BestMode;

	int OpaqueAlphaError;

	BlockError Error = BlockError(0, 0);
	BlockSSIM Quality = BlockSSIM(1.0, 1.0);

	uint64_t PersonalParameter;
	uint32_t PersonalMode;

	size_t unused[3];

	__m128i DataMask_I16[16];

	__m128i ImageRows_U8[4];
	__m128i MaskRows_S8[4];

	//

	Area Area1;

	bool LazyArea12[0x40], LazyArea22[0x40];
	Area Area12[0x40], Area22[0x40];

	bool LazyArea13[0x40], LazyArea23[0x40], LazyArea33[0x40];
	Area Area13[0x40], Area23[0x40], Area33[0x40];
};

struct alignas(8) Node
{
	int Error, Color;

	INLINED void Init(int error, int color) noexcept
	{
		Error = error;
		Color = color;
	}
};

// A,G,R,B
struct alignas(64) Modulations
{
	uint64_t Values_I16[16];

	int Ways[16];

	int Loops[16];

	int Best[16];
};

Area& GetArea(Area& area, bool& lazy, const Cell& cell, const uint64_t indices) noexcept;

int AreaGetBestPca3(Area& area) noexcept;

void AreaReduceTable2(const Area& area, __m128i& mc, uint64_t& indices) noexcept;
void AreaReduceTable3(const Area& area, __m128i& mc, uint64_t& indices) noexcept;
void AreaReduceTable4(__m128i& mc, uint64_t& indices) noexcept;

NOTINLINED Node* radix_sort(Node* input, Node* work, size_t N) noexcept;

NOTINLINED int ComputeSubsetTable(const Area& area, const __m128i mweights, const __m128i mfix, Modulations& state, const int M) noexcept;

namespace Mode0 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode0

namespace Mode1 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode1

namespace Mode2 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode2

namespace Mode3 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode3

namespace Mode4 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode4

namespace Mode5 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode5

namespace Mode6 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode6

namespace Mode7 {

	void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

	void FinalPackBlock(uint8_t output[16], Cell& input) noexcept;

	void CompressBlockFast(Cell& input) noexcept;

	void CompressBlock(Cell& input) noexcept;

	void CompressBlockFull(Cell& input) noexcept;

	void PrintCounters() noexcept;

} // namespace Mode7

void DecompressBlock(uint8_t input[16], Cell& output) noexcept;

void CompressBlock(uint8_t output[16], Cell& input) noexcept;

void CompressStatistics();

extern bool gDoDraft, gDoFast, gDoNormal, gDoSlow;
