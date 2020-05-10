#pragma once

#include "pch.h"

static ALWAYS_INLINED uint64_t InsertZeroBit(uint64_t value, int index) noexcept
{
	uint64_t high = (~0ui64 << index) & value;

	return (value ^ high) + (high << 1);
}

static ALWAYS_INLINED uint64_t RemoveZeroBit(uint64_t value, int index) noexcept
{
	uint64_t high = (~0ui64 << index) & value;

	return (value ^ high) + (high >> 1);
}
