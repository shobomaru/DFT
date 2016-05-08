#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

using sig_t = float;
using fft_t = complex<float>;
using sigd = vector<sig_t>;
using fftd = vector<fft_t>;

const float PI = 3.1415926f;

int main(int argc, char** argv)
{
	const int X = 32;
	const int Y = 32;

	// Generate data
	sigd data;
	data.resize(X * Y);
	for (auto y = 0u; y < Y; y++)
	{
		for (auto x = 0u; x < X; x++)
		{
			auto i = y * X + x;
			//data[i] = sinf(2.0f * PI * (float)(x + y) / 4);
			data[i] = (x == 0 && y == 0) ? (float)(X * Y) : 0.0f;
		}
	}

	// Fourier transform @ X-axis
	fftd fd;
	fd.resize(data.size());
	for (auto y = 0u; y < Y; y++)
	{
		for (auto x = 0u; x < X; x++)
		{
			auto i = y * X + x;
			auto& c = fd[i];
			for (auto j = 0u; j < X; j++)
			{
				auto h = y * X + j;
#if 0
				float r = data[j] * cosf(2.0f * PI * j * x / X);
				float m = -data[j] * sinf(2.0f * PI * j * x / X);
				c += fft_t(r, m);
#else
				c += data[h] * exp(fft_t(0, -2.0f * PI * j * x / X));
#endif
			}
		}
	}

	// Fourier transform @ Y-axis
	for (auto x = 0u; x < X; x++)
	{
		for (auto y = 0u; y < Y; y++)
		{
			auto i = y * X + x;
			auto& c = fd[i];
			for (auto j = 0u; j < Y; j++)
			{
				auto h = j * X + x;
#if 0
				float r = data[j] * cosf(2.0f * PI * j * y / Y);
				float m = -data[j] * sinf(2.0f * PI * j * y / Y);
				c += fft_t(r, m);
#else
				c += data[h] * exp(fft_t(0, -2.0f * PI * j * y / Y));
#endif
			}
		}
	}

	// Calcurate spectrum
	sigd ampl;
	ampl.resize(data.size());
	transform(fd.cbegin(), fd.cend(), ampl.begin(), [X, Y](auto& c) {
		return sqrtf(norm(c)) / sqrtf((float)X * (float)Y);
	});

	// Output CSV
	ofstream ofs("_original.csv", ios::trunc);
	if (ofs)
	{
		for (auto y = 0u; y < Y; y++)
		{
			for (auto x = 0u; x < X; x++)
			{
				int i = y * X + x;
				ofs << data[i] << ",";
			}
			ofs << endl;
		}
	}
	ofs = ofstream("_out.csv", ios::trunc);
	if (ofs)
	{
		for (auto y = 0u; y < Y; y++)
		{
			for (auto x = 0u; x < X; x++)
			{
				int i = y * X + x;
				ofs << ampl[i] << ",";
			}
			ofs << endl;
		}
	}

	return 0;
}
