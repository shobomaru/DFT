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
	const int N = 32;
	const int N_BITS = 5;

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

	// Fourier transform (Cooley-Tukey, Decimation in frequency)
	float the = 1.0f;
	for (auto m = data.size(); (m / 2) >= 1; m /= 2)
	{
		for (auto i = 0u; i < (m / 2); i++)
		{
#if 0
			float re = cosf(-2.0f * PI * i * the / N);
			float im = sinf(-2.0f * PI * i * the / N);
#else
			auto w = exp(fft_t(0, -2.0f * PI * i * the / N));
#endif
			for (auto j = i; j < data.size(); j += m)
			{
				auto k = j + (m / 2);
				fft_t x = fd[j] - fd[k];
				fd[j] += fd[k];
#if 0
				fd[k] = fft_t(x.real() * re - x.imag() * im, x.imag() * re + x.real() * im);
#else
				fd[k] = fft_t(x.real() * w.real() - x.imag() * w.imag(), x.imag() * w.real() + x.real() * w.imag());
#endif
			}
		}
		the *= 2.0f;
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
	transform(fd.cbegin(), fd.cend(), ampl.begin(), [](auto& c) {
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
