// Bc7Compress.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "Bc7Core.h"
#include "Bc7Tables.h"
#include "Bc7Pca.h"
#include "IO.h"
#include "Worker.h"

#include <chrono>
#include <string>
#include <vector>

static ALWAYS_INLINED int Max(int x, int y) noexcept
{
	return (x > y) ? x : y;
}

static INLINED void FullMask(uint8_t* mask_agrb, int stride, int src_w, int src_h) noexcept
{
	for (int y = 0; y < src_h; y++)
	{
		int* w = (int*)&mask_agrb[y * stride];

		for (int x = 0; x < src_w; x++)
		{
			w[x] = -1;
		}
	}
}

static INLINED void ComputeAlphaMaskWithOutline(uint8_t* mask_agrb, const uint8_t* src_bgra, int stride, int src_w, int src_h, int radius)
{
	if (radius < 0)
		return;

	const int full_w = radius + src_w + radius;

	// row buffer to accumulate Y sum
	int* buf = new int[full_w + size_t(15)];
	memset(buf, 0, full_w * sizeof(int));

	// add kernel rows except last
	for (int y = 0; y < radius; y++)
	{
		const uint32_t* __restrict r = reinterpret_cast<const uint32_t*>(&src_bgra[y * stride]);

		for (int i = 0; i < src_w; i++)
		{
			const uint8_t v = uint8_t((r[i] & 0xFF000000u) != 0);
			buf[radius + i] += v;
		}
	}

	for (int y = 0; y < src_h; y++)
	{
		// add last kernel row
		if (y + radius < src_h)
		{
			const uint32_t* __restrict r = reinterpret_cast<const uint32_t*>(&src_bgra[(y + radius) * stride]);

			for (int i = 0; i < src_w; i++)
			{
				const uint8_t v = uint8_t((r[i] & 0xFF000000u) != 0);
				buf[radius + i] += v;
			}
		}

		int32_t* __restrict p = reinterpret_cast<int32_t*>(&mask_agrb[y * stride]);

		// filter current row
		int sum = 0;
		for (int i = 0; i < radius + radius; i++)
		{
			sum -= buf[i];
		}
		for (int i = 0; i < src_w; i++)
		{
			// include last kernel pixel
			sum -= buf[radius + radius + i];

			const int32_t v = (sum < 0) ? -1 : 0;

			// remove first kernel pixel
			sum += buf[i];

			p[i] = v | 0xFF;
		}

		// subtract first kernel row
		if (y >= radius)
		{
			const uint32_t* __restrict r = reinterpret_cast<const uint32_t*>(&src_bgra[(y - radius) * stride]);

			for (int i = 0; i < src_w; i++)
			{
				const uint8_t v = uint8_t((r[i] & 0xFF000000u) != 0);
				buf[radius + i] -= v;
			}
		}
	}

	delete[] buf;
}

static void PackTexture(const IBc7Core& bc7Core, uint8_t* dst_bc7, uint8_t* src_bgra, uint8_t* mask_agrb, int stride, int src_w, int src_h, PBlockKernel blockKernel, size_t block_size, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim)
{
	auto start = std::chrono::high_resolution_clock::now();

	ProcessTexture(dst_bc7, src_bgra, mask_agrb, stride, src_w, src_h, blockKernel, block_size, pErrorAlpha, pErrorColor, pssim);

	auto finish = std::chrono::high_resolution_clock::now();

	if (blockKernel == bc7Core.pCompress)
	{
		int pixels = src_h * src_w;

		int span = Max((int)std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count(), 1);

		int kpx_s = pixels / span;

		PRINTF("    Compressed %d blocks, elapsed %i ms, throughput %d.%03d Mpx/s", pixels >> 4, span, kpx_s / 1000, kpx_s % 1000);

#if defined(OPTION_COUNTERS)
		CompressStatistics();
#endif
	}
}

static INLINED void VisualizePartitionsGRB(uint8_t* dst_bc7, int size)
{
	for (int i = 0; i < size; i += 16)
	{
		uint64_t data0 = *(const uint64_t*)&dst_bc7[i + 0];
		uint64_t data1 = *(const uint64_t*)&dst_bc7[i + 8];

		if (data0 & 0xF)
		{
			if (data0 & 3)
			{
				if (data0 & 1)
				{
					data0 &= (1u << (4 + 1)) - 1u;

					data0 |=
						(0xFuLL << 29) + // G0
						(0xFuLL << 13); // R2

					data1 = 
						(0xFuLL << 5); // B4
				}
				else
				{
					data0 &= (1u << (6 + 2)) - 1u;

					data0 |=
						(0x3FuLL << 32) + // G0
						(0x3FuLL << 20); // R2

					data1 = 0;
				}
			}
			else
			{
				if (data0 & 4)
				{
					data0 &= (1u << (6 + 3)) - 1u;

					data0 |=
						(0x1FuLL << 39) + // G0
						(0x1FuLL << 19); // R2

					data1 =
						(0x1FuLL << 25); // B4
				}
				else
				{
					data0 &= (1u << (6 + 4)) - 1u;

					data0 |=
						(0x7FuLL << 38) + // G0
						(0x7FuLL << 24); // R2

					data1 = 0;
				}
			}
		}
		else if (data0 & 0xF0)
		{
			if (data0 & 0x30)
			{
				if (data0 & 0x10)
				{
					data0 &= (1u << 5) - 1u;

					data0 |=
						(0x1FuLL << 18) + // G0
						(0x1FuLL << 8) + // R0
						(0x3FuLL << 38); // A0

					data1 = 0;
				}
				else
				{
					data0 &= (1u << 6) - 1u;

					data0 |=
						(0x7FuLL << 22) + // G0
						(0x7FuLL << 8) + // R0
						(0xFFuLL << 50); // A0

					data1 = 0;
				}
			}
			else
			{
				if (data0 & 0x40)
				{
					if ((data0 == 0x40) && (data1 == 0))
						continue;

					data0 &= (1u << 7) - 1u;

					// data1 = InsertZeroBit(data1 >> 1, 3);
					data1 = (data1 & ~0xFuLL) + ((data1 & 0xFuLL) >> 1);

					if (((data1 >> 2) ^ data1) & 0x3333333333333333uLL)
					{
						data0 |=
							(0x50uLL << 21) + // G0
							(0x50uLL << 7) + // R0
							(0x50uLL << 35) + // B0
							(0x7FuLL << 49); // A0
					}
					else
					{
						data0 |=
							(0x30uLL << 21) + // G0
							(0x30uLL << 7) + // R0
							(0x30uLL << 35) + // B0
							(0x7FuLL << 49); // A0
					}

					data1 = 0;
				}
				else
				{
					data0 &= (1u << (6 + 8)) - 1u;

					data0 |=
						(0x1FuLL << 34) + // G0
						(0x1FuLL << 24); // R2

					data1 =
						(0x1FuLL << 10) + // A0
						(0x1FuLL << 20); // A2
				}
			}
		}

		*(uint64_t*)&dst_bc7[i + 0] = data0;
		*(uint64_t*)&dst_bc7[i + 8] = data1;
	}
}

int Bc7MainWithArgs(const IBc7Core& bc7Core, const std::vector<std::string>& args)
{
	bool doDraft = true;
	bool doNormal = true;
	bool doSlow = false;

	bool flip = true;
	bool mask = true;
	int border = 1;
	bool linearData = false;

	const char* src_name = nullptr;
	const char* dst_name = nullptr;
	const char* result_name = nullptr;
	const char* partitions_name = nullptr;
	const char* bad_name = nullptr;

	for (int i = 0, n = (int)args.size(); i < n; i++)
	{
		const char* arg = args[i].c_str();

		if (arg[0] == '/')
		{
			if (strcmp(arg, "/compare") == 0)
			{
				doDraft = false;
				doNormal = false;
				doSlow = false;
				continue;
			}
			else if (strcmp(arg, "/draft") == 0)
			{
				doDraft = true;
				doNormal = false;
				doSlow = false;
				continue;
			}
			else if (strcmp(arg, "/slow") == 0)
			{
				doDraft = true;
				doNormal = true;
				doSlow = true;
				continue;
			}
			else if (strcmp(arg, "/noflip") == 0)
			{
				flip = false;
				continue;
			}
			else if (strcmp(arg, "/nomask") == 0)
			{
				mask = false;
				continue;
			}
			else if (strcmp(arg, "/retina") == 0)
			{
				border = 2;
				continue;
			}
			else if (strcmp(arg, "/data") == 0)
			{
				mask = false;
				linearData = true;
				continue;
			}
			else if (strcmp(arg, "/debug") == 0)
			{
				if (++i < n)
				{
					result_name = args[i].c_str();
				}
				continue;
			}
			else if (strcmp(arg, "/map") == 0)
			{
				if (++i < n)
				{
					partitions_name = args[i].c_str();
				}
				continue;
			}
			else if (strcmp(arg, "/bad") == 0)
			{
				if (++i < n)
				{
					bad_name = args[i].c_str();
				}
				continue;
			}
#if defined(WIN32)
			else
			{
				PRINTF("Unknown %s", arg);
				continue;
			}
#endif
		}

		if (src_name == nullptr)
		{
			src_name = arg;
		}
		else if (dst_name == nullptr)
		{
			dst_name = arg;
		}
		else
		{
			PRINTF("Error: %s", arg);
			return 1;
		}
	}

	if (!src_name)
	{
		PRINTF("No input");
		return 1;
	}

	uint8_t* src_image_bgra;
	int src_image_w, src_image_h;

	if (!ReadImage(src_name, src_image_bgra, src_image_w, src_image_h, flip))
	{
		PRINTF("Problem with image %s", src_name);
		return 1;
	}

	PRINTF("Loaded %s", src_name);

	int src_texture_w = (Max(4, src_image_w) + 3) & ~3;
	int src_texture_h = (Max(4, src_image_h) + 3) & ~3;

	if (Max(src_texture_w, src_texture_h) > 16384)
	{
		PRINTF("Huge image %s", src_name);
		delete[] src_image_bgra;
		return 1;
	}

	int c = 4;
	int src_image_stride = src_image_w * c;
	int src_texture_stride = src_texture_w * c;

	uint8_t* src_texture_bgra = new uint8_t[src_texture_h * src_texture_stride];

	for (int i = 0; i < src_image_h; i++)
	{
		memcpy(&src_texture_bgra[i * src_texture_stride], &src_image_bgra[i * src_image_stride], src_image_stride);

		for (int j = src_image_stride; j < src_texture_stride; j += c)
		{
			memcpy(&src_texture_bgra[i * src_texture_stride + j], &src_image_bgra[i * src_image_stride + src_image_stride - c], c);
		}
	}

	for (int i = src_image_h; i < src_texture_h; i++)
	{
		memcpy(&src_texture_bgra[i * src_texture_stride], &src_texture_bgra[(src_image_h - 1) * src_texture_stride], src_texture_stride);
	}

	PRINTF("  Image %dx%d, Texture %dx%d", src_image_w, src_image_h, src_texture_w, src_texture_h);

	delete[] src_image_bgra, src_image_bgra = nullptr;

	uint8_t* dst_texture_bgra = new uint8_t[src_texture_h * src_texture_stride];

	memcpy(dst_texture_bgra, src_texture_bgra, src_texture_h* src_texture_stride);

	int Size = src_texture_h * src_texture_w;

	// https://www.khronos.org/opengles/sdk/tools/KTX/file_format_spec/
	uint32_t head[16 + 7 + 1];
	head[0] = 0x58544BABu; // identifier
	head[1] = 0xBB313120u; // identifier
	head[2] = 0x0A1A0A0Du; // identifier
	head[3] = 0x04030201u; // endianness
	head[4] = 0; // glType
	head[5] = 1; // glTypeSize
	head[6] = 0; // glFormat
	head[7] = 0x8E8C; // glInternalFormat = COMPRESSED_RGBA_BPTC_UNORM_ARB
	head[8] = 0x1908; // glBaseInternalFormat = GL_RGBA
	head[9] = static_cast<uint32_t>(src_image_w); // pixelWidth
	head[10] = static_cast<uint32_t>(src_image_h); // pixelHeight
	head[11] = 0; // pixelDepth
	head[12] = 0; // numberOfArrayElements
	head[13] = 1; // numberOfFaces
	head[14] = 1; // numberOfMipmapLevels
	head[15] = 0x1C; // bytesOfKeyValueData
	head[16] = 0x17; // keyAndValueByteSize
	head[17] = 0x4F58544Bu; // "KTXOrientation"
	head[18] = 0x6E656972u;
	head[19] = 0x69746174u;
	head[20] = 0x53006E6Fu; // "S=r,T=u"
	head[21] = 0x542C723Du;
	head[22] = flip ? 0x00753Du : 0x00643Du;
	head[23] = static_cast<uint32_t>(Size); // imageSize

	bc7Core.pInitTables(doDraft, doNormal, doSlow, linearData);

	if ((dst_name != nullptr) && dst_name[0])
	{
		uint8_t* mask_agrb = new uint8_t[src_texture_h * src_texture_stride];

		if (mask)
		{
			ComputeAlphaMaskWithOutline(mask_agrb, src_texture_bgra, src_texture_stride, src_texture_w, src_texture_h, border);
		}
		else
		{
			FullMask(mask_agrb, src_texture_stride, src_texture_w, src_texture_h);
		}

		uint8_t* dst_bc7 = new uint8_t[Size];
		memset(dst_bc7, 0, Size);

		LoadBc7(dst_name, sizeof(head), dst_bc7, Size);

		int64_t mse_alpha = 0;
		int64_t mse_color = 0;
		BlockSSIM ssim = BlockSSIM(0, 0);
		PackTexture(bc7Core, dst_bc7, src_texture_bgra, mask_agrb, src_texture_stride, src_texture_w, src_texture_h, bc7Core.pCompress, 16, mse_alpha, mse_color, ssim);

		int pixels = src_texture_h * src_texture_w;

		if (mse_alpha > 0)
		{
			const double weightAlpha = bc7Core.pGetWeightAlpha();

			PRINTF("      SubTexture A qMSE = %.1f, qPSNR = %f, SSIM_4x4 = %.8f",
				(1.0 / weightAlpha) * mse_alpha / pixels,
				10.0 * log((255.0 * 255.0) * weightAlpha * pixels / mse_alpha) / log(10.0),
				ssim.Alpha * 16.0 / pixels);
		}
		else
		{
			PRINTF("      Whole A");
		}

		if (mse_color > 0)
		{
			const double weightColor = bc7Core.pGetWeightColor();

			PRINTF("      SubTexture RGB qMSE = %.1f, qPSNR = %f, wSSIM_4x4 = %.8f",
				(1.0 / weightColor) * mse_color / pixels,
				10.0 * log((255.0 * 255.0) * weightColor * pixels / mse_color) / log(10.0),
				ssim.Color * 16.0 / pixels);
		}
		else
		{
			PRINTF("      Whole RGB");
		}

		SaveBc7(dst_name, (const uint8_t*)head, sizeof(head), dst_bc7, Size);

#if !defined(OPTION_LIBRARY)
		PackTexture(bc7Core, dst_bc7, dst_texture_bgra, mask_agrb, src_texture_stride, src_texture_w, src_texture_h, bc7Core.pDecompress, 16, mse_alpha, mse_color, ssim);

		if ((bad_name != nullptr) && bad_name[0])
		{
			ShowBadBlocks(src_texture_bgra, dst_texture_bgra, mask_agrb, src_texture_stride, src_texture_w, src_texture_h);

			WriteImage(bad_name, mask_agrb, src_texture_w, src_texture_h, flip);
		}

		if ((partitions_name != nullptr) && partitions_name[0])
		{
			VisualizePartitionsGRB(dst_bc7, Size);

			PackTexture(bc7Core, dst_bc7, mask_agrb, dst_texture_bgra, src_texture_stride, src_texture_w, src_texture_h, bc7Core.pDecompress, 16, mse_alpha, mse_color, ssim);

			WriteImage(partitions_name, mask_agrb, src_texture_w, src_texture_h, flip);
		}
#else
		(void)bad_name;
		(void)partitions_name;
#endif

		delete[] dst_bc7;
		delete[] mask_agrb;
	}

	if ((result_name != nullptr) && result_name[0])
	{
#if !defined(OPTION_LIBRARY)
		WriteImage(result_name, dst_texture_bgra, src_texture_w, src_texture_h, flip);
#endif
	}

	delete[] dst_texture_bgra;
	delete[] src_texture_bgra;

	return 0;
}

#if !defined(OPTION_LIBRARY) && defined(WIN32)

bool GetBc7Core(void* bc7Core);

int __cdecl main(int argc, char* argv[])
{
	IBc7Core bc7Core{};
	if (!GetBc7Core(&bc7Core))
	{
		PRINTF("Unsupported CPU");
		return 2;
	}

	if (argc < 2)
	{
		PRINTF("Usage: Bc7Compress [/retina] [/nomask] [/noflip] [/data] src");
		PRINTF("                   [/draft | /slow] [dst.ktx]");
		PRINTF("                   [/debug result.png] [/map partitions.png] [/bad bad.png]");
		return 1;
	}

	std::vector<std::string> args;
	args.reserve(argc);

	for (int i = 1; i < argc; i++)
	{
		args.emplace_back(argv[i]);
	}

	return Bc7MainWithArgs(bc7Core, args);
}

#endif
