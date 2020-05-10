
#include "pch.h"
#include "Worker.h"
#include "Metrics.h"

#if defined(WIN32)
#include <windows.h>
#endif

#include <stdexcept>
#include <thread>
#include <atomic>

static constexpr int WorkerThreadStackSize = 2 * 1024 * 1024;

static ALWAYS_INLINED int Max(int x, int y) noexcept
{
	return (x > y) ? x : y;
}

class WorkerJob
{
public:
	WorkerJob* _Next;

protected:
	int _Count;
	WorkerItem _Items[0x100];

public:
	WorkerJob()
		: _Next(nullptr)
		, _Count(0)
	{
	}

	bool Add(const WorkerItem& item)
	{
		_Items[_Count++] = item;

		return _Count >= (int)(sizeof(_Items) / sizeof(_Items[0]));
	}

	WorkerItem* begin()
	{
		return _Items;
	}

	WorkerItem* end()
	{
		return _Items + _Count;
	}
};

class Worker
{
protected:
#if defined(WIN32)
	CRITICAL_SECTION _Sync;
	HANDLE _Done;
#else
	std::mutex _Sync;
#endif

	PBlockKernel _BlockKernel;
	int _Stride;

	WorkerJob* _First;
	WorkerJob* _Last;

	int64_t _errorAlpha;
	int64_t _errorColor;
	BlockSSIM _ssim;

	std::atomic_int _Running;

public:
	Worker()
		: _errorAlpha(0)
		, _errorColor(0)
		, _ssim(0, 0)
	{
#if defined(WIN32)
		if (!InitializeCriticalSectionAndSpinCount(&_Sync, 1000))
			throw std::runtime_error("init");

		_Done = CreateEvent(NULL, FALSE, FALSE, NULL);
		if (_Done == nullptr)
			throw std::runtime_error("init");
#endif

		_BlockKernel = nullptr;
		_Stride = 0;

		_First = nullptr;
		_Last = nullptr;

		_Running = 0;
	}

	~Worker()
	{
		for (WorkerJob* job; (job = _First) != nullptr;)
		{
			_First = job->_Next;

			delete job;
		}

		_Last = nullptr;

#if defined(WIN32)
		if (_Done != nullptr)
			CloseHandle(_Done), _Done = nullptr;

		DeleteCriticalSection(&_Sync);
#endif
	}

	void Lock()
	{
#if defined(WIN32)
		EnterCriticalSection(&_Sync);
#else
		_Sync.lock();
#endif
	}

	void UnLock()
	{
#if defined(WIN32)
		LeaveCriticalSection(&_Sync);
#else
		_Sync.unlock();
#endif
	}

	void Add(WorkerJob* job)
	{
		if (_Last)
			_Last->_Next = job;
		else
			_First = job;

		_Last = job;
	}

protected:
	WorkerJob* Take()
	{
		Lock();

		WorkerJob* job = _First;
		if (job)
		{
			_First = job->_Next;

			if (_First == nullptr)
				_Last = nullptr;
		}

		UnLock();

		return job;
	}

#if defined(WIN32)
	static DWORD WINAPI ThreadProc(LPVOID lpParameter)
#else
	static int ThreadProc(Worker* lpParameter)
#endif
	{
		Worker* worker = static_cast<Worker*>(lpParameter);

		int64_t errorAlpha = 0;
		int64_t errorColor = 0;
		BlockSSIM ssim = BlockSSIM(0, 0);

		for (WorkerJob* job; (job = worker->Take()) != nullptr;)
		{
			worker->_BlockKernel(job->begin(), job->end(), worker->_Stride, errorAlpha, errorColor, ssim);

			delete job;
		}

		worker->Lock();

		worker->_errorAlpha += errorAlpha;
		worker->_errorColor += errorColor;

		worker->_ssim.Alpha += ssim.Alpha;
		worker->_ssim.Color += ssim.Color;

		worker->UnLock();

		worker->_Running--;

#if defined(WIN32)
		if (worker->_Running <= 0)
		{
			SetEvent(worker->_Done);
		}
#endif

		return 0;
	}

public:
	void Run(PBlockKernel blockKernel, int stride, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim)
	{
		_BlockKernel = blockKernel;
		_Stride = stride;

		_errorAlpha = 0;
		_errorColor = 0;
		_ssim = BlockSSIM(0, 0);

		int n = Max(1, (int)std::thread::hardware_concurrency());
		_Running = n;

		for (int i = 0; i < n; i++)
		{
#if defined(WIN32)
			HANDLE hthread = CreateThread(NULL, WorkerThreadStackSize, ThreadProc, this, 0, NULL);
			if (hthread == nullptr)
				throw std::runtime_error("fork");
			CloseHandle(hthread);
#else
			//std::thread thread(std::bind(ThreadProc, this));
			boost::thread::attributes attrs;
			attrs.set_stack_size(WorkerThreadStackSize);
			boost::thread thread(attrs, boost::bind(ThreadProc, this));
			thread.detach();
#endif
		}

#if defined(WIN32)
		WaitForSingleObject(_Done, INFINITE);
#else
		for (;; )
		{
			std::this_thread::yield();

			if (_Running <= 0)
				break;
		}
#endif

		pErrorAlpha = _errorAlpha;
		pErrorColor = _errorColor;
		pssim = _ssim;
	}
};

static ALWAYS_INLINED __m128i ConvertBgraToAgrb(__m128i mc) noexcept
{
	const __m128i mrot = _mm_set_epi8(
		0 + 12, 2 + 12, 1 + 12, 3 + 12,
		0 + 8, 2 + 8, 1 + 8, 3 + 8,
		0 + 4, 2 + 4, 1 + 4, 3 + 4,
		0, 2, 1, 3);

	return _mm_shuffle_epi8(mc, mrot);
}

void DecompressKernel(const WorkerItem* begin, const WorkerItem* end, int stride, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim) noexcept
{
	(void)pErrorAlpha;
	(void)pErrorColor;
	(void)pssim;

	Cell output;

	for (auto it = begin; it != end; it++)
	{
		DecompressBlock(it->_Output, output);

		{
			const uint8_t* p = it->_Cell;

			_mm_storeu_si128((__m128i*)p, ConvertBgraToAgrb(output.ImageRows_U8[0]));

			p += stride;

			_mm_storeu_si128((__m128i*)p, ConvertBgraToAgrb(output.ImageRows_U8[1]));

			p += stride;

			_mm_storeu_si128((__m128i*)p, ConvertBgraToAgrb(output.ImageRows_U8[2]));

			p += stride;

			_mm_storeu_si128((__m128i*)p, ConvertBgraToAgrb(output.ImageRows_U8[3]));
		}
	}
}

void CompressKernel(const WorkerItem* begin, const WorkerItem* end, int stride, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim) noexcept
{
	Cell input;

	for (auto it = begin; it != end; it++)
	{
		{
			const uint8_t* p = it->_Cell;

			input.ImageRows_U8[0] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

			p += stride;

			input.ImageRows_U8[1] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

			p += stride;

			input.ImageRows_U8[2] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

			p += stride;

			input.ImageRows_U8[3] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));
		}

		{
			const uint8_t* p = it->_Mask;

			input.MaskRows_S8[0] = _mm_loadu_si128((const __m128i*)p);

			p += stride;

			input.MaskRows_S8[1] = _mm_loadu_si128((const __m128i*)p);

			p += stride;

			input.MaskRows_S8[2] = _mm_loadu_si128((const __m128i*)p);

			p += stride;

			input.MaskRows_S8[3] = _mm_loadu_si128((const __m128i*)p);
		}

		CompressBlock(it->_Output, input);

		pErrorAlpha += input.Error.Alpha;
		pErrorColor += input.Error.Total - input.Error.Alpha;

		pssim.Alpha += input.Quality.Alpha;
		pssim.Color += input.Quality.Color;
	}
}

void ProcessTexture(uint8_t* dst, uint8_t* src_bgra, uint8_t* mask_agrb, int stride, int src_w, int src_h, PBlockKernel blockKernel, size_t block_size, int64_t& pErrorAlpha, int64_t& pErrorColor, BlockSSIM& pssim)
{
	Worker worker;

	WorkerJob* job = new WorkerJob();

	uint8_t* output = dst;

	for (int y = 0; y < src_h; y += 4)
	{
		uint8_t* cell = src_bgra + y * stride;
		uint8_t* mask = mask_agrb + y * stride;

		for (int x = 0; x < src_w; x += 4)
		{
			if (job->Add(WorkerItem(output, cell, mask)))
			{
				worker.Add(job);

				job = new WorkerJob();
			}

			output += block_size;
			cell += 16;
			mask += 16;
		}
	}

	worker.Add(job);

	worker.Run(blockKernel, stride, pErrorAlpha, pErrorColor, pssim);
}

void ShowBadBlocks(const uint8_t* src_bgra, const uint8_t* dst_bgra, uint8_t* mask_agrb, int stride, int src_w, int src_h) noexcept
{
	Cell input;
	Cell output;

	for (int y = 0; y < src_h; y += 4)
	{
		const uint8_t* src = src_bgra + y * stride;
		const uint8_t* dst = dst_bgra + y * stride;
		uint8_t* mask = mask_agrb + y * stride;

		for (int x = 0; x < src_w; x += 4)
		{
			{
				const uint8_t* p = src;

				input.ImageRows_U8[0] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

				p += stride;

				input.ImageRows_U8[1] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

				p += stride;

				input.ImageRows_U8[2] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

				p += stride;

				input.ImageRows_U8[3] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));
			}

			{
				const uint8_t* p = mask;

				input.MaskRows_S8[0] = _mm_loadu_si128((const __m128i*)p);

				p += stride;

				input.MaskRows_S8[1] = _mm_loadu_si128((const __m128i*)p);

				p += stride;

				input.MaskRows_S8[2] = _mm_loadu_si128((const __m128i*)p);

				p += stride;

				input.MaskRows_S8[3] = _mm_loadu_si128((const __m128i*)p);
			}

			{
				const uint8_t* p = dst;

				output.ImageRows_U8[0] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

				p += stride;

				output.ImageRows_U8[1] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

				p += stride;

				output.ImageRows_U8[2] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));

				p += stride;

				output.ImageRows_U8[3] = ConvertBgraToAgrb(_mm_loadu_si128((const __m128i*)p));
			}

			bool bad = DetectGlitches(input, output);
			if (bad)
			{
				const uint8_t* r = src;
				uint8_t* p = mask;

				_mm_storeu_si128((__m128i*)p, _mm_loadu_si128((const __m128i*)r));

				r += stride;
				p += stride;

				_mm_storeu_si128((__m128i*)p, _mm_loadu_si128((const __m128i*)r));

				r += stride;
				p += stride;

				_mm_storeu_si128((__m128i*)p, _mm_loadu_si128((const __m128i*)r));

				r += stride;
				p += stride;

				_mm_storeu_si128((__m128i*)p, _mm_loadu_si128((const __m128i*)r));
			}
			else
			{
				__m128i mc = _mm_setzero_si128();

				uint8_t* p = mask;

				_mm_storeu_si128((__m128i*)p, mc);

				p += stride;

				_mm_storeu_si128((__m128i*)p, mc);

				p += stride;

				_mm_storeu_si128((__m128i*)p, mc);

				p += stride;

				_mm_storeu_si128((__m128i*)p, mc);
			}

			src += 16;
			dst += 16;
			mask += 16;
		}
	}
}
