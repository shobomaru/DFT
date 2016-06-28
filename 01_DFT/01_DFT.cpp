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

int main(int argc, char** argv)
{
	const int N = 32;

	// Generate data
	sigd data;
	data.resize(N);
	for (auto i = 0u; i < data.size(); i++)
	{
		data[i] = sinf(2.0f * PI * float(i) / 4); // Sine wave
		//data[i] = (i % 4) >= 2 ? 0.0f : 1.0f; // Square wave (duty ratio 50%)
		//data[i] = (i % 4) / 3.0f; // Saw wave
	}

	// Fourier transform
	fftd fd;
	fd.resize(data.size());
	for (auto i = 0u; i < data.size(); i++)
	{
		auto& c = fd[i];
		for (auto j = 0u; j < data.size(); j++)
		{
#if 0
			float r = data[j] * cosf(2.0f * PI * j * i / N);
			float m = -data[j] * sinf(2.0f * PI * j * i / N);
			c += fft_t(r, m);
#else
			c += data[j] * exp(fft_t(0, -2.0f * PI * j * i / N));
#endif
		}
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
