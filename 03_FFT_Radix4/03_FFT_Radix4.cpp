#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

using sig_t = float;
using fft_t = complex<float>;
using sigd = vector<sig_t>;
using fftd = vector<fft_t>;

const float PI = 3.1415926f;

unsigned int reverse_bits(unsigned int b)
{
	const auto a1 = 0b01010101010101010101010101010101u;
	const auto a2 = 0b00110011001100110011001100110011u;
	const auto a3 = 0b00001111000011110000111100001111u;
	const auto a4 = 0b00000000000000001111111100000000u;
	b = ((b & a1) << 1) | ((b >> 1) & a1);
	b = ((b & a2) << 2) | ((b >> 2) & a2);
	b = ((b & a3) << 4) | ((b >> 4) & a3);
	b = (b << 24) | ((b & a4) << 8) | ((b >> 8) & a4) | (b >> 24);
	return b;
}

int main(int argc, char** argv)
{
	const int N = 64;
	const int N_BITS = 6;

	// Generate data
	sigd data;
	data.resize(N);
	for (auto i = 0u; i < data.size(); i++)
	{
		data[i] = sinf(2.0f * PI * float(i) / 4);
	}

	// Real number to complex number
	fftd fd;
	fd.resize(data.size());
	copy(data.cbegin(), data.cend(), fd.begin());

	// Fourier transform (Cooley-Tukey, Radix 4, Decimation in frequency)
	float the = 1.0f;
	for (auto m = data.size(); (m / 4) >= 1; m /= 4)
	{
		for (auto i = 0u; i < (m / 4); i++)
		{
#if 0
			float re1 = cosf(1 * -2.0f * PI * i * the / N);
			float im1 = sinf(1 * -2.0f * PI * i * the / N);
			float re2 = cosf(2 * -2.0f * PI * i * the / N);
			float im2 = sinf(2 * -2.0f * PI * i * the / N);
			float re3 = cosf(3 * -2.0f * PI * i * the / N);
			float im3 = sinf(3 * -2.0f * PI * i * the / N);
#else
			auto w1 = exp(fft_t(0, 1 * -2.0f * PI * i * the / N));
			auto w2 = exp(fft_t(0, 2 * -2.0f * PI * i * the / N));
			auto w3 = exp(fft_t(0, 3 * -2.0f * PI * i * the / N));
#endif
			for (auto j = i; j < data.size(); j += m)
			{
				auto k0 = j;
				auto k1 = k0 + (m / 4);
				auto k2 = k1 + (m / 4);
				auto k3 = k2 + (m / 4);
				fft_t x0 = fd[k0] + fd[k2];
				fft_t x1 = fd[k0] - fd[k2];
				fft_t x2 = fd[k3] + fd[k1];
				fft_t x3 = fd[k3] - fd[k1];
				fft_t y0 = x0 + x2;
				fft_t y1 = x0 - x2;
				fft_t y2 = fft_t(x1.real() - x3.real(), x1.imag() + x3.imag());
				fft_t y3 = fft_t(x1.real() + x3.real(), x1.imag() - x3.imag());
#if 0
				fd[k0] = y0;
				fd[k1] = fft_t(re2 * y1.real() - im2 * y1.imag(), re2 * y1.imag() + im2 * y1.real());
				fd[k2] = fft_t(re1 * y2.real() - im1 * y2.imag(), re1 * y2.imag() + im1 * y2.real());
				fd[k3] = fft_t(re3 * y3.real() - im3 * y3.imag(), re3 * y3.imag() + im3 * y3.real());
#else
				fd[k0] = y0;
				fd[k1] = fft_t(w2.real() * y1.real() - w2.imag() * y1.imag(), w2.real() * y1.imag() + w2.imag() * y1.real());
				fd[k2] = fft_t(w1.real() * y2.real() - w1.imag() * y2.imag(), w1.real() * y2.imag() + w1.imag() * y2.real());
				fd[k3] = fft_t(w3.real() * y3.real() - w3.imag() * y3.imag(), w3.real() * y3.imag() + w3.imag() * y3.real());
#endif
			}
		}
		the *= 4.0f;
	}

	// Bit reversal
	fftd temp;
	temp.resize(data.size());
	copy(fd.cbegin(), fd.cend(), temp.begin());
	for (auto i = 0u; i < data.size(); i++)
	{
		unsigned int b = reverse_bits(i) >> (32 - N_BITS);
		fd[b] = temp[i];
	}

	// Calcurate spectrum
	sigd ampl;
	ampl.resize(N);
	transform(fd.cbegin(), fd.cend(), ampl.begin(), [N](auto& c) {
		return sqrtf(norm(c)) / N;
	});

	// Estimate frequency
	int hz = 0;
	sig_t p = 0;
	for (auto i = 1u; i < data.size() / 2; i++)
	{
		if (ampl[i] > p)
		{
			p = ampl[i];
			hz = i;
		}
	}
	cout << "Freq = " << hz << "[Hz]" << endl;

	return 0;
}
