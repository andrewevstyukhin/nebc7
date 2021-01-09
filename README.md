# nebc7 - Nearly Exhaustive BC7 compressor

Nebc7 is an open-source software published as is.

Nebc7 converts specified RGBA image to BC7 format. It focuses on image quality and uses weighted PSNR metric with some assistance of SSIM.

Nebc7 always preserves opaque alpha for opaque blocks.

## Details

Nebc7 computes bounding box for each partition to choose the most suitable mode. Trivial min & max bounds swap to conform channel's covariations. Such mode is very fast and can be enabled by "/draft" command-line switch. Usually only -3dB below of an ideal encoding.

Hi-quality BC7 encodings can waste space for noise contained in source images. Error calculation function shifts away least bit of |delta| to allow simple denoising (of rounded pixels). While indices selection function treats full 8-bit values preserving gradients. Unseen random delta with magnitude 1 gives 48.13dB cap. Shifting trick forces unusual naming for metrics: qMSE and qPSNR.

It is suitable to use low entropy modes for opaque blocks. Mode 6 has excellent block structure, while modes 4&5 can discard single index. Nebc7 tries modes in special order and uses error stepping for only valuable switching. Random deltas with magnitude 3 are usually invisible in non-gradient areas and give 38.58dB cap. An interesting fact: many monitors can show only 6-bit true color. It is possible to render modes map specifying "/map ..." command-line parameter. Mode 6 is painted gray, 4&5 - yellow, other modes use green, red, blue for partitions 0, 1, 2 accordingly.

Nebc7 sorts nearly all possible endpoint values of each single channel and chooses some good of them. Then exhaustive search tries all selected combinations and reveals the best solution. Sorting provides a surprisingly fast convergence for such slow process.

Modes 7, 1, 3 are memory-bound because of large tables, they partially limited in default working mode. Slow modes can be fully activated by impractical "/slow" command-line switch.

For premultiplied alpha it is necessary to specify "/nomask" command-line option. While extruded RGBA images can highly benefit from masking. Switch "/retina" allows future artifact-free scaling by 0.5. Masking gives smaller compressed images and better borders, because masked pixels can have any value.

Many encoders use single metric (RMSE / MSE / PSNR) insensitive to a direction. While SSIM is unhandy for direct compression, it enhances correlation when encoding produces equal deltas and so SSIM overcomes dithering.

## Usage

The solution was tested on SSSE3, SSE4.1, AVX, AVX2, AVX-512BW - capable CPUs for Win64 API only.

`Bc7Compress /nomask /noflip source.png destination.ktx [/debug result.png]`

I would recommend using AVX2 or AVX-512 for the best performance. See Bc7Mode.h about settings.

## Customization

The kDenoise and kDenoiseStep constants in Bc7Mode.h define RDO capabilities:

1. The most performant and compact mode is defined by kDenoise = 1, kDenoiseStep = 3 \* 3.
2. Generally precise but much slower mode is defined by kDenoise = 0, kDenoiseStep = 0.
3. Especially for Moon shots (black backround with low light amplitude) I recommend kDenoise = 0, kDenoiseStep = 3 \* 3.

The constants kAlpha, kGreen, kRed, kBlue set weights for channels and depend on the nature of the data.

## Example

Recompressing "BC7Ltest.png" (gained from https://code.google.com/archive/p/nvidia-texture-tools/downloads bc7_export.zip) on i7-6700 CPU:

    Bc7Compress.exe /nomask /noflip BC7Ltest.png output.ktx /debug output.png
    Loaded BC7Ltest.png
      Image 152x152, Texture 152x152
        Compressed 1444 blocks, elapsed 11 ms, throughput 2.100 Mpx/s
          SubTexture A qMSE = 0.0, qPSNR = 69.026097, SSIM_4x4 = 0.99997870
          SubTexture RGB qMSE = 0.6, qPSNR = 50.685068, wSSIM_4x4 = 0.99964057
        Saved output.ktx
      Saved output.png

Compressing "frymire.png" (gained from https://github.com/castano/nvidia-texture-tools/blob/master/data/testsuite/waterloo/frymire.png):

    Bc7Compress.exe /nomask /noflip frymire.png frymire.ktx
    Loaded frymire.png
      Image 1118x1105, Texture 1120x1108
        Compressed 77560 blocks, elapsed 77 ms, throughput 16.116 Mpx/s
          Whole A
          SubTexture RGB qMSE = 0.5, qPSNR = 50.950326, wSSIM_4x4 = 0.97143024
        Saved frymire.ktx

Compressing "8192.png" (gained from https://bitbucket.org/wolfpld/etcpak/downloads/8192.png) in development mode:

    Bc7Compress.exe /draft /nomask /noflip 8192.png 8192.ktx
    Loaded 8192.png
      Image 8192x8192, Texture 8192x8192
        Compressed 4194304 blocks, elapsed 529 ms, throughput 126.859 Mpx/s
          Whole A
          SubTexture RGB qMSE = 0.1, qPSNR = 56.778563, wSSIM_4x4 = 0.99534311
        Saved 8192.ktx

## Copyright

Copyright (c) 2019-2020 Andrew Evstyukhin

Licensed under the MIT License.

All other trademarks are property of their respective owners.
