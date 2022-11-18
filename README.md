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

## Special features

Additional mode 6 with 2-bit index instead of 4-bit gives a significant reduction in size.

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

## Examples

Let's try the soup from https://cbloomrants.blogspot.com/2020/06/oodle-texture-sample-run.html on an i7-7800X CPU with AVX-512 enabled:

  A precise compression with kDenoise = 0, kDenoiseStep = 0

    Bc7Compress.exe /nomask /noflip mysoup1024.png mysoup1024-precise.ktx /debug mysoup1024-precise.png
      Image 1024x1024, Texture 1024x1024
        Compressed 65536 blocks, elapsed 1396 ms, throughput 0.751 Mpx/s
          Whole A
          SubTexture RGB qMSE = 1.6, qPSNR = 46.144434, wSSIM_4x4 = 0.99468383

    dssim-1.13.exe mysoup1024.png mysoup1024-precise.png
      0.00011307      mysoup1024-precise.png

    mysoup1024-precise.7z  973 652

  Standard BC7 modes separated by error stepping kDenoise = 1, kDenoiseStep = 3 \* 3:

    Bc7Compress.exe /nomask /noflip mysoup1024.png mysoup1024-stepping.ktx /debug mysoup1024-stepping.png
      Image 1024x1024, Texture 1024x1024
        Compressed 65536 blocks, elapsed 345 ms, throughput 3.039 Mpx/s
          Whole A
          SubTexture RGB qMSE = 1.6, qPSNR = 46.097026, wSSIM_4x4 = 0.99034363

    dssim-1.13.exe mysoup1024.png mysoup1024-stepping.png
      0.00023433      mysoup1024-stepping.png

    mysoup1024-stepping.7z  840 465

  Using additional RDO Mode6Index2 (Bc7Core.cpp, Mode6Index2::CompressBlockFast(input)):

    Bc7Compress.exe /nomask /noflip mysoup1024.png mysoup1024-mode6index2.ktx /debug mysoup1024-mode6index2.png
      Image 1024x1024, Texture 1024x1024
        Compressed 65536 blocks, elapsed 350 ms, throughput 2.995 Mpx/s
          Whole A
          SubTexture RGB qMSE = 2.1, qPSNR = 44.949510, wSSIM_4x4 = 0.98476649

    dssim-1.13.exe mysoup1024.png mysoup1024-mode6index2.png
      0.00045838      mysoup1024-mode6index2.png

    mysoup1024-mode6index2.7z  651 851

Then check all modes with BC7Ltest.png from https://code.google.com/archive/p/nvidia-texture-tools/downloads bc7_export.zip:

    Bc7Compress.exe /nomask /noflip BC7Ltest.png BC7Ltest.ktx /debug output.png
      Image 152x152, Texture 152x152
        Compressed 1444 blocks, elapsed 17 ms, throughput 1.359 Mpx/s
          SubTexture A qMSE = 0.0, qPSNR = 69.026097, SSIM_4x4 = 0.99997865
          SubTexture RGB qMSE = 0.6, qPSNR = 50.687018, wSSIM_4x4 = 0.99964438

    dssim-1.13.exe BC7Ltest.png output.png
      0.00003303      output.png

And finally compress an interesting image https://github.com/castano/image-datasets/blob/master/waterloo/frymire.png:

    Bc7Compress.exe /nomask /noflip frymire.png frymire.ktx /debug output.png
      Image 1118x1105, Texture 1120x1108
        Compressed 77560 blocks, elapsed 69 ms, throughput 17.984 Mpx/s
          Whole A
          SubTexture RGB qMSE = 0.6, qPSNR = 50.594453, wSSIM_4x4 = 0.97119160

    frymire.7z  313 550

For non-production purposes we can apply a "draft" mode on earlily saved https://bitbucket.org/wolfpld/etcpak/downloads/8192.png:

    Bc7Compress.exe /draft /nomask /noflip 8192.png 8192.ktx /debug output.png
      Image 8192x8192, Texture 8192x8192
        Compressed 4194304 blocks, elapsed 282 ms, throughput 237.974 Mpx/s
          Whole A
          SubTexture RGB qMSE = 0.9, qPSNR = 48.734921, wSSIM_4x4 = 0.98443459

    dssim-1.13.exe 8192.png output.png
      0.00085877      output.png

    8192.7z  31 309 975

## Copyright

Copyright (c) 2019-2021 Andrew Evstyukhin

Licensed under the MIT License.

All other trademarks are property of their respective owners.
