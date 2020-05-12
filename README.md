# nebc7 - Nearly Exhaustive BC7 compressor

Nebc7 is an open-source software published as is.

Nebc7 converts specified RGBA image to BC7 format. It focuses on image quality and uses weighted PSNR metric with some assistance of SSIM.

Nebc7 always preserves opaque alpha for opaque blocks.

## Usage

The solution was tested on SSSE3, SSE4.1, AVX, AVX2, AVX512BW - capable CPUs for Win64 API only.

`Bc7Compress /nomask /noflip source.png destination.ktx [/debug result.png]`

I would recommend using AVX2 for the best performance. See Bc7Mode.h about settings.

## Example

Recompressing "BC7Ltest.png" (gained from https://code.google.com/archive/p/nvidia-texture-tools/downloads bc7_export.zip) on i7-6700 CPU:

    Bc7Compress.exe /slow /nomask /noflip BC7Ltest.png output.ktx /debug output.png
    Loaded BC7Ltest.png
      Image 152x152, Texture 152x152
        Compressed 1444 blocks, elapsed 26 ms, throughput 0.888 Mpx/s
          SubTexture A MSE = 0.0, PSNR = 73.986163, SSIM_4x4 = 0.99999923
          SubTexture RGB wMSE = 0.0, wPSNR = 62.172358, wSSIM_4x4 = 0.99999238
        Saved output.ktx
      Saved output.png

Compressing "frymire.png" (gained from https://github.com/castano/nvidia-texture-tools/blob/master/data/testsuite/waterloo/frymire.png):

    Bc7Compress.exe /nomask /noflip frymire.png frymire.ktx
    Loaded frymire.png
      Image 1118x1105, Texture 1120x1108
        Compressed 77560 blocks, elapsed 449 ms, throughput 2.763 Mpx/s
          Exactly A
          SubTexture RGB wMSE = 0.2, wPSNR = 55.181449, wSSIM_4x4 = 0.99980677
        Saved frymire.ktx

Compressing "frymire.png" in development mode:

    Bc7Compress.exe /draft /nomask /noflip frymire.png frymire.ktx
    Loaded frymire.png
      Image 1118x1105, Texture 1120x1108
        Compressed 77560 blocks, elapsed 141 ms, throughput 8.801 Mpx/s
          Exactly A
          SubTexture RGB wMSE = 0.4, wPSNR = 52.056761, wSSIM_4x4 = 0.99952034
        Saved frymire.ktx

Compressing "8192.png" (gained from https://bitbucket.org/wolfpld/etcpak/downloads/8192.png) in development mode:

    Bc7Compress.exe /draft /nomask /noflip 8192.png 8192.ktx
    Loaded 8192.png
      Image 8192x8192, Texture 8192x8192
        Compressed 4194304 blocks, elapsed 12377 ms, throughput 5.422 Mpx/s
          Exactly A
          SubTexture RGB wMSE = 0.4, wPSNR = 52.364416, wSSIM_4x4 = 0.99625929
        Saved 8192.ktx

## Copyright

Copyright (c) 2019-2020 Andrew Evstyukhin

Licensed under the MIT License.

All other trademarks are property of their respective owners.
