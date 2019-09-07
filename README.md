# nebc7 - Nearly Exhaustive BC7 compressor

Nebc7 is an experimental solution published as is.

Nebc7 converts specified RGBA image to BC7 format. It focuses on image quality and uses weighted PSNR metric with some assistance of SSIM.

Nebc7 always preserves opaque alpha for opaque blocks.

## Usage

The solution was tested on AVX-capable CPU for Win64 API only.

`Bc7Compress /nomask /noflip source.png destination.ktx [/debug result.png]`

I would recommend using AVX2 for the best performance. See Bc7Mode.h about settings.

## Example

Recompressing "BC7Ltest.png" (gained from https://code.google.com/archive/p/nvidia-texture-tools/downloads bc7_export.zip) on i7-6700 CPU:

    Bc7Compress.exe /slow /nomask /noflip BC7Ltest.png output.ktx /debug output.png
    Loaded BC7Ltest.png
      Image 152x152, Texture 152x152
        Compressed 1444 blocks, elapsed 39 ms, throughput 0.592 Mpx/s
          SubTexture A MSE = 0.0, PSNR = 73.986163, SSIM_4x4 = 0.99999923
          SubTexture RGB wMSE = 0.0, wPSNR = 62.172358, wSSIM_4x4 = 0.99999238
        Saved output.ktx
      Saved output.png

Compressing "frymire.png" (gained from https://github.com/castano/nvidia-texture-tools/blob/master/data/testsuite/waterloo/frymire.png):

    Bc7Compress.exe /nomask /noflip frymire.png frymire.ktx
    Loaded frymire.png
      Image 1118x1105, Texture 1120x1108
        Compressed 77560 blocks, elapsed 512 ms, throughput 2.423 Mpx/s
          Exactly A
          SubTexture RGB wMSE = 0.2, wPSNR = 55.181449, wSSIM_4x4 = 0.99980677
        Saved frymire.ktx

## Copyright

Copyright (c) 2019 Andrew Evstyukhin

Licensed under the MIT License.

All other trademarks are property of their respective owners.
