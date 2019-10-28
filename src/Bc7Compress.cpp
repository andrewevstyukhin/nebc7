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

static INLINED int Max(int x, int y) noexcept
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
	int full_w = 1 + radius + src_w + radius;
	int full_h = 1 + radius + src_h + radius;

	uint8_t* data = new uint8_t[full_h * full_w];
	memset(data, 0, full_h * full_w);

	for (int y = 0; y < src_h; y++)
	{
		const uint8_t* r = &src_bgra[y * stride + 3];
		uint8_t* w = &data[(y + radius + 1) * full_w + (radius + 1)];

		for (int x = 0; x < src_w; x++)
		{
			w[x] = uint8_t(r[x * 4] != 0);
		}
	}

	int* sum = new int[full_h * full_w];
	memset(sum, 0, full_h * full_w * sizeof(int));

	int from_py = (radius + 1) * full_w;
	int to_py = full_h * full_w;
	for (int py = from_py; py < to_py; py += full_w)
	{
		int prev_py = py - full_w;

		for (int x = radius + 1; x < full_w; x++)
		{
			sum[py + x] = sum[prev_py + x] - sum[prev_py + x - 1] + data[py + x] + sum[py + x - 1];
		}
	}

	int a = radius + radius + 1;
	for (int y = 0; y < src_h; y++)
	{
		int* w = (int*)&mask_agrb[y * stride];

		const int* rL = &sum[y * full_w];
		const int* rH = &sum[(y + a) * full_w];

		for (int x = 0; x < src_w; x++, rL++, rH++)
		{
			int v = rH[a] - *rH + *rL - rL[a];

			w[x] = (v != 0) ? -1 : 0xFF;
		}
	}

	delete[] sum, sum = nullptr;

	delete[] data, data = nullptr;
}

static void PackTexture(uint8_t* dst_bc7, uint8_t* src_bgra, uint8_t* mask_agrb, int stride, int src_w, int src_h, PBlockKernel blockKernel, size_t block_size, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim)
{
	auto start = std::chrono::high_resolution_clock::now();

	ProcessTexture(dst_bc7, src_bgra, mask_agrb, stride, src_w, src_h, blockKernel, block_size, pErrorAlpha, pErrorColor, pssim);

	auto finish = std::chrono::high_resolution_clock::now();

	if (blockKernel != DecompressKernel)
	{
		int pixels = src_h * src_w;

		int span = Max((int)std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count(), 1);

		int kpx_s = pixels / span;

		PRINTF("    Compressed %d blocks, elapsed %i ms, throughput %d.%03d Mpx/s", pixels >> 4, span, kpx_s / 1000, kpx_s % 1000);

#if defined(OPTION_COUNTERS)
		if (blockKernel == CompressKernel)
		{
			CompressStatistics();
		}
#endif
	}
}

int Bc7MainWithArgs(const std::vector<std::string>& args)
{
	gDoDraft = true;
	gDoFast = true;
	gDoNormal = true;
	gDoSlow = false;

	bool flip = true;
	bool mask = true;
	int border = 1;

	const char* src_name = nullptr;
	const char* dst_name = nullptr;
	const char* result_name = nullptr;

	for (int i = 0, n = (int)args.size(); i < n; i++)
	{
		const char* arg = args[i].c_str();

		if (arg[0] == '/')
		{
			if (strcmp(arg, "/compare") == 0)
			{
				gDoDraft = false;
				gDoFast = false;
				gDoNormal = false;
				gDoSlow = false;
				continue;
			}
			else if (strcmp(arg, "/draft") == 0)
			{
				gDoDraft = true;
				gDoFast = false;
				gDoNormal = false;
				gDoSlow = false;
				continue;
			}
			else if (strcmp(arg, "/fast") == 0)
			{
				gDoDraft = true;
				gDoFast = true;
				gDoNormal = false;
				gDoSlow = false;
				continue;
			}
			else if (strcmp(arg, "/normal") == 0)
			{
				gDoDraft = true;
				gDoFast = true;
				gDoNormal = true;
				gDoSlow = false;
				continue;
			}
			else if (strcmp(arg, "/slow") == 0)
			{
				gDoDraft = true;
				gDoFast = true;
				gDoNormal = true;
				gDoSlow = true;
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
			else if (strcmp(arg, "/debug") == 0)
			{
				if (++i < n)
				{
					result_name = args[i].c_str();
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

	InitInterpolation();
	InitShrinked();
	InitSelection();
	InitLevels();

#if defined(OPTION_PCA)
	InitPCA();
#endif

	memcpy(dst_texture_bgra, src_texture_bgra, src_texture_h * src_texture_stride);

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
		PackTexture(dst_bc7, src_texture_bgra, mask_agrb, src_texture_stride, src_texture_w, src_texture_h, &CompressKernel, 16, mse_alpha, mse_color, ssim);

		int pixels = src_texture_h * src_texture_w;

		if (mse_alpha > 0)
		{
			PRINTF("      SubTexture A MSE = %.1f, PSNR = %f, SSIM_4x4 = %.8f",
				(1.0 / kAlpha) * mse_alpha / pixels,
				10.0 * log((255.0 * 255.0) * kAlpha * pixels / mse_alpha) / log(10.0),
				ssim.Alpha * 16.0 / pixels);
		}
		else
		{
			PRINTF("      Exactly A");
		}

		if (mse_color > 0)
		{
#if defined(OPTION_LINEAR)
			PRINTF("      SubTexture RGB MSE = %.1f, PSNR = %f, SSIM_4x4 = %.8f",
				(1.0 / kColor) * mse_color / pixels,
				10.0 * log((255.0 * 255.0) * kColor * pixels / mse_color) / log(10.0),
				ssim.Color * 16.0 / pixels);
#else
			PRINTF("      SubTexture RGB wMSE = %.1f, wPSNR = %f, wSSIM_4x4 = %.8f",
				(1.0 / kColor) * mse_color / pixels,
				10.0 * log((255.0 * 255.0) * kColor * pixels / mse_color) / log(10.0),
				ssim.Color * 16.0 / pixels);
#endif
		}
		else
		{
			PRINTF("      Exactly RGB");
		}

		SaveBc7(dst_name, (const uint8_t*)head, sizeof(head), dst_bc7, Size);

		PackTexture(dst_bc7, dst_texture_bgra, mask_agrb, src_texture_stride, src_texture_w, src_texture_h, &DecompressKernel, 16, mse_alpha, mse_color, ssim);

		delete[] dst_bc7;
		delete[] mask_agrb;
	}

	if ((result_name != nullptr) && result_name[0])
	{
		WriteImage(result_name, dst_texture_bgra, src_texture_w, src_texture_h, flip);
	}

	delete[] dst_texture_bgra;
	delete[] src_texture_bgra;

	return 0;
}

#if defined(WIN32)

int __cdecl main(int argc, char* argv[])
{
	if (argc < 2)
	{
		PRINTF("Usage: Bc7Compress [/fast | /normal | /slow] [/retina] [/nomask] [/noflip] src");
		PRINTF("                   [dst.ktx] [/debug result.png]");
		return 1;
	}

	std::vector<std::string> args;
	args.reserve(argc);

	for (int i = 1; i < argc; i++)
	{
		args.emplace_back(argv[i]);
	}

	return Bc7MainWithArgs(args);
}

#endif
