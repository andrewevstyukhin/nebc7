#pragma once

#include "pch.h"
#include "Bc7Core.h"

struct WorkerItem
{
	uint8_t* _Output;
	uint8_t* _Cell;
	uint8_t* _Mask;

	WorkerItem()
	{
	}

	WorkerItem(uint8_t* output, uint8_t* cell, uint8_t* mask)
		: _Output(output)
		, _Cell(cell)
		, _Mask(mask)
	{
	}
};

using PBlockKernel = void(*)(const WorkerItem* begin, const WorkerItem* end, int stride, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim) noexcept;

void DecompressKernel(const WorkerItem* begin, const WorkerItem* end, int stride, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim) noexcept;
void CompressKernel(const WorkerItem* begin, const WorkerItem* end, int stride, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim) noexcept;

void ProcessTexture(uint8_t* dst, uint8_t* src_bgra, uint8_t* mask_agrb, int stride, int src_w, int src_h, PBlockKernel blockKernel, size_t block_size, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim);

void ShowBadBlocks(const uint8_t* src_bgra, const uint8_t* dst_bgra, uint8_t* mask_agrb, int stride, int src_w, int src_h) noexcept;
