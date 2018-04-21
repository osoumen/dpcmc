//
//  main.c
//  dpcmc
//
//  Created by osoumen on 2018/04/03.
//  Copyright © 2018年 test. All rights reserved.
//

#include "FileRead.h"
#include "FileWrite.h"
#include "Stk.h"
#include <cassert>
#include <algorithm>
#include <list>
#include <random>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cstdint>

extern "C" {
	void rdft(int, int, double *, int *, double *);
}


class NoiseShape {
public:
	NoiseShape(const double *coeff, int taps) {
		mCoeff = coeff;
		mTaps = taps;
		mErrPtr = 0;
		mErr = new double[taps];
		for (int i=0; i<taps; ++i) {
			mErr[i] = .0;
		}
	}
	~NoiseShape() {
		delete [] mErr;
	}
	double getOutput() {
		double err = .0;
		for (int i=0; i < mTaps; ++i) {
			err += mCoeff[i] * mErr[(mErrPtr + i) % mTaps];
		}
		return err;
	}
	void inputSample(double sample) {
		mErrPtr--;
		if (mErrPtr < 0) {
			mErrPtr = mTaps - 1;
		}
		mErr[mErrPtr] = sample;
	}
	
private:
	NoiseShape();
	
	const double* mCoeff;
	int			mTaps;
	double		*mErr;
	int			mErrPtr;
};

class Convolution {
public:
	Convolution(int filterLength, const double *filter)
	{
		m_filter_len = filterLength;
		m_fft_len = filterLength * 2;
		m_fft_n = 1;
		while ((2 << m_fft_n) < m_fft_len) {
			m_fft_n++;
		}
		m_ip_len = 2 + (2 << ((m_fft_n - 1)/2));
		m_w_len = 2 << (m_fft_n - 1);
		m_ip = new int[m_ip_len];
		m_work = new double[m_w_len];
		m_filter = new double[m_fft_len];
		m_fft_buff = new double[m_fft_len];
		m_filter_result = new double[m_fft_len];
		
		m_ip[0] = 0;
		
		for (int i=0; i<m_fft_len; ++i) {
			m_filter_result[i] = .0;
		}
		
		for (int i=0; i<m_filter_len; ++i) {
			m_filter[i] = filter[i] * (1.0 / m_filter_len);
		}
		for (int i=m_filter_len; i<m_fft_len; ++i) {
			m_filter[i] = .0;
		}
		rdft(m_fft_len, 1, m_filter, m_ip, m_work);
	}
	~Convolution()
	{
		delete [] m_ip;
		delete [] m_work;
		delete [] m_filter;
		delete [] m_fft_buff;
		delete [] m_filter_result;
	}
	
	void convolve(double *in, int length)
	{
		assert(length <= m_filter_len);
		
		for (int i=0; i<length; ++i) {
			m_fft_buff[i] = in[i];
		}
		for (int i=length; i<m_fft_len; ++i) {
			m_fft_buff[i] = .0;
		}
		rdft(m_fft_len, 1, m_fft_buff, m_ip, m_work);
		m_fft_buff[0] *= m_filter[0];
		m_fft_buff[1] *= m_filter[1];
		for (int i=1; i<m_filter_len; ++i) {
			double realp = m_fft_buff[2*i] * m_filter[2*i] - m_fft_buff[2*i+1] * m_filter[2*i+1];
			double imagp = m_fft_buff[2*i] * m_filter[2*i+1] + m_fft_buff[2*i+1] * m_filter[2*i];
			m_fft_buff[2*i] = realp;
			m_fft_buff[2*i+1] = imagp;
		}
		rdft(m_fft_len, -1, m_fft_buff, m_ip, m_work);
		for (int i=0; i<m_fft_len; ++i) {
			m_filter_result[i] += m_fft_buff[i];
		}
		for (int i=0; i<length; ++i) {
			in[i] = m_filter_result[i];
		}
		for (int i=0; i<m_filter_len; ++i) {
			m_filter_result[i] = m_filter_result[i+m_filter_len];
			m_filter_result[i+m_filter_len] = .0;
		}
	}
	
private:
	Convolution();
	
	int m_filter_len;
	int m_fft_n;
	int m_fft_len;
	double *m_filter;
	double *m_fft_buff;
	double *m_filter_result;
	
	int m_ip_len;
	int m_w_len;
	int *m_ip;
	double *m_work;
	
};

#if _WIN32
#define DIRECTORY_SEPARATOR_CHAR '\\'
#else
#define DIRECTORY_SEPARATOR_CHAR '/'
#endif

int resample(const stk::StkFrames &src,
			 double srcRate,
			 stk::StkFrames &dst,
			 int dstSamples,
			 double dstRate,
			 double stopband,
			 int window_len);
int firLowpassFilter(const stk::StkFrames &src,
					 stk::StkFrames &dst,
					 double stopband,
					 double initial_value);
int preprocessInputBuffer(const std::string &in_file_path,
						  stk::StkFrames &preprocess_buff, double dst_rate,
						  double shifter_weight, bool in_verbose_mode, int centerBiasLevel,
						  bool use_linearity_correction);
void getFileNameDeletingDirectory(const std::string &path, std::string &out);
void getFileNameDeletingExt(const std::string &path, std::string &out);
int outputDpcmToFile(const std::string &in_file_path,
					 const unsigned char *dpcm, int dpcmFrames, int dpcm_initial_volume);
int outputBufferToFile(const std::string &in_file_path, stk::StkFrames &buff,
					   const std::string &fnameSuffix, double sampleRate);
int encodeDpcm(stk::StkFrames &preprocess_buff, const int dpcm_initial_volume,
			   unsigned char *dpcm_buff, const int outBufSize,
			   int ditherMode, int noiseshapeMode, bool use_linearity_correction);
int decodeDpcm(const std::string &in_file_path, const unsigned char *dpcm,
			   int dpcm_frames, int dpcm_initial_volume, int bias_level, bool outputToFile,
			   double dst_rate);

inline unsigned char dpcmBitMask(int x)
{
	return (0x01 << (x % 8));
}

inline double linearToDacCurve(double input)
{
	return (159.79 / ((22638.0 / input) + 100)) / 0.5771522095;
}

double sinc(double p_x1)
{
	if ( p_x1==.0 ) return 1.0;
	else {
		double x1 = p_x1*M_PI;
		return (std::sin(x1)/x1);
	}
}

typedef struct {
	const double	*coeffs;
	int				taps;
	bool			is_iir;
} NoiseShapeFilter;


double wave_gain = 1.0;
int force_initial_volume = -1;
double shifter_weight = 1.0;
int sample_rate_ind = 15;
bool use_linearity_correction = true;
bool help_mode = false;
int dither_mode = 3;
int noise_shape_mode = 3;
int center_bias_level = 64;
bool no_resampling = false;

bool output_envelope = false;
bool output_preprocessed = false;
bool output_encoded_wav = false;

const double s_ns_coeffs_flat[] = {
	.0
};
const double s_ns_coeffs_lp[] = {
	-7.0/16.0
};
const double s_ns_coeffs_hp[] = {
	7.0/16.0
};
const double s_ns_coeffs_notch[] = {
	0.3, -0.2, -0.1, -0.2, -0.1, 0.04
//		0.2, -0.15, -0.1, -0.05, -0.02
//		0.2, -0.1, -0.1, -0.1, -0.05
};

NoiseShapeFilter	s_ns_filter[] = {
	{
		s_ns_coeffs_flat, 1, true
	},
	{
		s_ns_coeffs_lp, 1, true
	},
	{
		s_ns_coeffs_hp, 1, true
	},
	{
		s_ns_coeffs_notch, 6, true
	},
};

const int num_ns_filters = sizeof(s_ns_filter) / sizeof(NoiseShapeFilter);

const int	resample_window_len = 256;
const int	filter_window_len = 16384;

const double system_clock_rate_ntsc = 1789772.5;
const int sample_rate_cycles[16] = {
	428,
	380,
	340,
	320,
	286,
	254,
	226,
	214,
	190,
	160,
	142,
	128,
	106,
	85,
	72,
	54,
};

inline double dpcmSampleRateValue(int value) {
	return system_clock_rate_ntsc / sample_rate_cycles[value];
}

void getFileNameDeletingExt(const std::string &path, std::string &out)
{
	// 拡張子、パス除去処理
	size_t	len = path.length();
	size_t	extPos = len;
	size_t	bcPos = 0;
	
	extPos = path.find_last_of('.');
	
	bcPos = path.find_last_of(DIRECTORY_SEPARATOR_CHAR) + 1;
	
	if (bcPos > extPos) {
		extPos = bcPos;
	}
	
	out = path.substr(0, extPos);
}

void getFileNameDeletingDirectory(const std::string &path, std::string &out)
{
	// パス除去処理
	auto bcPos = path.find_last_of(DIRECTORY_SEPARATOR_CHAR) + 1;
	if (bcPos == std::string::npos) {
		out = path;
	}
	else {
		out = path.substr(bcPos, std::string::npos);
	}
}

int resample(const stk::StkFrames &src,
			 double srcRate,
			 stk::StkFrames &dst,
			 int dstSamples,
			 double dstRate,
			 double stopband,
			 int window_len)
{
	int				half_window_len = window_len / 2;
	double			srcSteps = srcRate / dstRate;
	double			cutoffRate = (1.0/srcSteps) < stopband?(1.0/srcSteps)*stopband:stopband;
	int				dstSize = dstSamples;
	int				actualDstSize = static_cast<int>(src.size() / srcSteps);
	
	if ( actualDstSize < dstSize ) {
		dstSize = actualDstSize;
	}
	for (int i=0; i<dstSize; i++) {
		double	src_pos = i*srcSteps;
		double	dstSum = .0;
		for (int j=-half_window_len; j<half_window_len; j++) {
			int	src_index = static_cast<int>(std::floor(src_pos+0.5)) + j;
			if ((src_index >= 0) && (src_index < (int)src.size())) {
				double	x = src_index - src_pos;
				double	window_x = x/half_window_len + 1.0;
				if (window_x < 0.0)	window_x = 0.0;
				if (window_x > 2.0)	window_x = 2.0;
				double	window = 0.54 - 0.46 * std::cos(M_PI * window_x);
				double	value = src[src_index] * sinc(x*cutoffRate) * window;
				dstSum += value * cutoffRate;
			}
		}
		dst[i] = dstSum;
	}
	return dstSize;
}

int firLowpassFilter(const stk::StkFrames &src,
					 stk::StkFrames &dst,
					 double stopband,
					 double initial_value)
{
	int		half_window_len = filter_window_len / 2;
	double	coeff[filter_window_len];
	double	input[filter_window_len];
	
	// フィルタ作成
	for (int i=0; i<filter_window_len; ++i) {
		double	x = i - half_window_len;
		double	window_x = x/half_window_len + 1.0;
		double	window = 0.54 - 0.46 * std::cos(M_PI * window_x);
		coeff[i] = sinc(x * stopband) * window * stopband;
	}
	
	Convolution conv(filter_window_len, coeff);

	// initial_valueの反映
	for (int i=0; i<filter_window_len; ++i) {
		input[i] = initial_value;
	}
	conv.convolve(input, filter_window_len);
	
	int src_index = 0;
	int dst_index = -half_window_len;
	
	while (dst_index < static_cast<int>(src.size())) {
		int input_len = static_cast<int>(src.size()) - src_index;
		if (input_len > filter_window_len) {
			input_len = filter_window_len;
		}
		if (input_len < 0) {
			input_len = 0;
		}
		for (int i=0; i<input_len; ++i) {
			input[i] = src[src_index + i];
		}
		for (int i=input_len; i<filter_window_len; ++i) {
			input[i] = 0;
		}
		conv.convolve(input, filter_window_len);
		src_index += input_len;
		
		for (int i=0; i<filter_window_len; ++i) {
			if ((dst_index >= 0) && (dst_index < static_cast<int>(src.size()))) {
				dst[dst_index] = input[i];
			}
			dst_index++;
		}
	}

	return (int)src.size();
}

int encodeDpcm(stk::StkFrames &preprocess_buff, const int dpcm_initial_volume,
			   unsigned char *dpcm_buff, const int outBufSize,
			   int ditherMode, int noiseshapeMode, bool use_linearity_correction)
{
	uint32_t src_frames = static_cast<uint32_t>(preprocess_buff.size());
	uint32_t dpcm_frames = ((src_frames + 127) / 128) * 128;
	if (dpcm_frames > (uint32_t)outBufSize*8)  {
		dpcm_frames = (uint32_t)outBufSize*8;
	}

	if (noise_shape_mode > num_ns_filters) {
		noise_shape_mode = num_ns_filters;
	}
	
	const double *ns_coeff = s_ns_filter[noiseshapeMode].coeffs;
	int ns_taps = s_ns_filter[noiseshapeMode].taps;
	bool is_iir = s_ns_filter[noiseshapeMode].is_iir;
	NoiseShape noiseShape(ns_coeff, ns_taps);
	
	const double noisefilter_coeff[] = {0.5};
	NoiseShape noiseFilter(noisefilter_coeff, 1);
	std::mt19937 rand_gen;
	int now_sample = dpcm_initial_volume;
	for (size_t i=0; i<dpcm_frames; ++i) {
		dpcm_buff[i / 8] &= ~dpcmBitMask((int)i);
		double original_sample = (i < src_frames)?
				(preprocess_buff[i]*128.0):center_bias_level;
		double target_sample = original_sample;
		
		// クリッピング
		if (original_sample > 126) original_sample = 126.;
		if (original_sample < 0) original_sample = .0;
		
		// ディザ
		if (ditherMode > 0) {
			if (ditherMode == 1) {
				// ディザ(RPDF)
				double noise = rand_gen();
				noise /= (double)0x80000000;
				noise -= 1.0;
				noise *= 0.5;
				target_sample += noise;
			}
			else {
				// ディザ(TPDF)
				double noise = .0;
				noise += (rand_gen() / (double)0x80000000);
				noise += (rand_gen() / (double)0x80000000);
				noise -= 2.0;
				noise *= 0.25;
				if (ditherMode == 3) {
					noise -= noiseFilter.getOutput();
					noiseFilter.inputSample(noise);
				}
				target_sample += noise;
			}
		}
		// ノイズシェーピング
		double dt = noiseShape.getOutput();
		if (is_iir) {
			target_sample -= dt;
			dt = .0;
		}
		else {
			target_sample += dt;
		}
		
		// 出力ビットと量子化誤差を求める
		if (target_sample >= now_sample) {
			dpcm_buff[i / 8] |= dpcmBitMask((int)i);
			if (now_sample <= 125) {
				now_sample += 2;
			}
		}
		else {
			if (now_sample >= 2) {
				now_sample -= 2;
			}
		}
		noiseShape.inputSample(now_sample - original_sample - dt);
	}
	
	return dpcm_frames;
}

int decodeDpcm(const std::string &in_file_path, const unsigned char *dpcm,
			   int dpcm_frames, int dpcm_initial_volume, int bias_level, bool outputToFile,
			   double dst_rate)
{
	int volume_max = 0;
	
	stk::StkFrames decode_sample(dpcm_frames, 1);
	int now_sample = dpcm_initial_volume;
	for (int i=0; i<dpcm_frames; ++i) {
		if (dpcm[i / 8] & dpcmBitMask(i)) {
			if (now_sample <= 125) {
				now_sample += 2;
			}
		}
		else {
			if (now_sample >= 2) {
				now_sample -= 2;
			}
		}
		
		if (volume_max < now_sample) {
			volume_max = now_sample;
		}
		
		// 実機非線形特性
		decode_sample[i] = linearToDacCurve(now_sample) - linearToDacCurve(bias_level);
	}
	if (outputToFile) {
		outputBufferToFile(in_file_path, decode_sample, "_dpcm", dst_rate);
	}
	
	return volume_max;
}

int outputDpcmToFile(const std::string &in_file_path, const unsigned char *dpcm, int dpcmFrames, int dpcm_initial_volume)
{
	std::string outfname;
	getFileNameDeletingExt(in_file_path, outfname);
	outfname += "_i";
	outfname += std::to_string(dpcm_initial_volume);
	outfname += ".dmc";
	std::cout << "output dmc file: " << outfname << std::endl;
	
	std::ofstream ofs(outfname, std::ios::out|std::ios::binary|std::ios::trunc);
	if (!ofs.bad()) {
		ofs.write((char*)dpcm, dpcmFrames/8);
	}
	else {
		std::cout << "can not open file: " << outfname << std::endl;
		return 1;
	}
	return 0;
}

int outputBufferToFile(const std::string &in_file_path, stk::StkFrames &buff, const std::string &fnameSuffix, double sampleRate)
{
	std::string out_fname;
	try {
		getFileNameDeletingExt(in_file_path, out_fname);
		out_fname += fnameSuffix;
		
		buff.setDataRate(sampleRate);
		stk::FileWrite outfile(out_fname);
		outfile.write(buff);
	}
	catch(stk::StkError err) {
		std::cout << "can not open file: " << out_fname << std::endl;
		return 1;
	}
	return 0;
}

int preprocessInputBuffer(const std::string &in_file_path,
						  stk::StkFrames &preprocess_buff, double dst_rate,
						  double shifter_weight,
						  int centerBiasLevel, bool use_linearity_correction)
{
	int	src_frames = (int)preprocess_buff.size();
	
	// thresholdの補正
	double centerBiasLevelf = centerBiasLevel;
	if (use_linearity_correction) {
		centerBiasLevelf = 128.0 * linearToDacCurve(centerBiasLevelf);
	}
	
	if (centerBiasLevelf > 63) {
		centerBiasLevelf -= 128;
	}
	double threshold_level = -centerBiasLevelf / 64.0;
	
	// 波形を下半分のみ残して半波整流する
	stk::StkFrames rectify_buff(src_frames, 1);
	for (int i=0; i<src_frames; ++i) {
		if (centerBiasLevelf >= 0) {
			if (preprocess_buff[i] >= threshold_level) {
				rectify_buff[i] = 0;
			}
			else {
				rectify_buff[i] = threshold_level - preprocess_buff[i];
			}
		}
		else {
			if (preprocess_buff[i] >= threshold_level) {
				rectify_buff[i] = preprocess_buff[i] - threshold_level;
			}
			else {
				rectify_buff[i] = 0;
			}
		}
	}
	
	// 先頭2048サンプルの平均値をとって初期値とする
	double start_sample = .0;
	{
		int window_len = (src_frames > 2048)?2048 : src_frames;
		for (int i=0; i<window_len; ++i) {
			start_sample += rectify_buff[i] / static_cast<double>(window_len);
		}
	}
	
	// DPCM変換後に近い波形になるように値が2/128ずつだけ変化するようにする
	{
		int now_sample = 2*static_cast<int>((force_initial_volume >= 0)?(force_initial_volume/2):(start_sample * 64.0));
		for (int i=0; i<src_frames; ++i) {
			double target = rectify_buff[i] * 128;
			if (threshold_level != 1.0) {
				target /= 1.0 - std::fabs(threshold_level);
			}
			if (target > now_sample) {
				now_sample += 2;
			}
			else {
				now_sample -= 2;
			}
			if (now_sample < 0) now_sample = 0;
			if (now_sample > 127) now_sample = 127;
			rectify_buff[i] = now_sample / 128.0;
		}
	}
	
	// ピーク&ホールド
	{
		double env = .0;
		double coeffRR = std::exp(-250.0 / dst_rate);
		for (int i=0; i<src_frames; ++i) {
			if (rectify_buff[i] >= env) {
				env = rectify_buff[i];
			}
			else {
				env *= coeffRR;
				env += (1.0 - coeffRR) * rectify_buff[i];
				rectify_buff[i] = env;
			}
		}
	}
	
//	{
//		stk::StkFrames	progress_out_buff((uint32_t)rectify_buff.size(), 1);
//		for (size_t i=0; i<rectify_buff.size(); ++i) {
//			progress_out_buff[i] = rectify_buff[i];
//			if (progress_out_buff[i] > 1.0) progress_out_buff[i] = 1.0;
//			if (progress_out_buff[i] < -1.0) progress_out_buff[i] = -1.0;
//		}
//		outputBufferToFile(in_file_path, progress_out_buff, "_rectify", dst_rate);
//	}
	
	// 整流波形にLPFを掛けて可聴成分を除去
	stk::StkFrames	filtered_shifter_buff(src_frames, 1);
	firLowpassFilter(rectify_buff,
					 filtered_shifter_buff,
					 20.0 / (dst_rate/2),
					 start_sample);
	// ゲイン調整
	for (int i=0; i<src_frames; ++i) {
		filtered_shifter_buff[i] *= 2.0;
	}
	
	if (output_envelope) {
		// エンベロープ波形
		stk::StkFrames	progress_out_buff((uint32_t)filtered_shifter_buff.size(), 1);
		for (size_t i=0; i<filtered_shifter_buff.size(); ++i) {
			progress_out_buff[i] = filtered_shifter_buff[i] * 0.5;
			if (progress_out_buff[i] > 1.0) progress_out_buff[i] = 1.0;
			if (progress_out_buff[i] < -1.0) progress_out_buff[i] = -1.0;
		}
		outputBufferToFile(in_file_path, progress_out_buff, "_envelope", dst_rate);
	}
	
	// 平滑された整流波形を加算してゲイン調整
	for (int i=0; i<src_frames; ++i) {
		filtered_shifter_buff[i] *= shifter_weight;
		if (centerBiasLevelf >= 0) {
			preprocess_buff[i] += filtered_shifter_buff[i] - threshold_level;
		}
		else {
			preprocess_buff[i] += 2.0 - threshold_level - filtered_shifter_buff[i];
		}
		preprocess_buff[i] *= 0.5;		// シフトによって倍になったピークを1.0に落とす
	}
	
	// 非線形DAC特性の補正
	if (use_linearity_correction) {
		for (int i=0; i<src_frames; ++i) {
			preprocess_buff[i] = 22368.0 / (128 * ((159.79 / (0.5771522095 * preprocess_buff[i])) - 100));
		}
	}
	
	return 0;
}

int processDmcInputFile(const std::string &in_file_path)
{
	std::ifstream ifs(in_file_path, std::ios::in | std::ios::binary);
	
	if (ifs.bad()) {
		std::cout << "can not open file: " << in_file_path << std::endl;
		exit(1);
	}
	
	ifs.seekg(0, std::ios::end);
	int dpcm_bytes = static_cast<int>(ifs.tellg());
	ifs.seekg(0, std::ios::beg);
	unsigned char *dpcm = new unsigned char[dpcm_bytes];
	ifs.read((char *)dpcm, dpcm_bytes);
	
	// dmcファイル出力
	std::string in_file_name;
	getFileNameDeletingDirectory(in_file_path, in_file_name);

	stk::Stk::setSampleRate(dpcmSampleRateValue(sample_rate_ind));
	
	decodeDpcm(in_file_name, dpcm,
			   (int)dpcm_bytes*8, force_initial_volume, center_bias_level, true,
			   dpcmSampleRateValue(sample_rate_ind));
	
	delete [] dpcm;

	return 0;
}

int processInputFile(const std::string &in_file_path)
{
	std::string in_file_name;
	getFileNameDeletingDirectory(in_file_path, in_file_name);
	
	std::cout << "----------------------------------" << std::endl;
	std::cout << "input file: " << in_file_path << std::endl;
	
	// 入力ファイルを読み込む
	stk::FileRead in_file;
	try {
		in_file.open(in_file_path);
	}
	catch(stk::StkError err) {
		std::cout << "can not open file: " << in_file_path << std::endl;
		exit(1);
	}
	uint32_t		in_frames = static_cast<uint32_t>(in_file.fileSize());
	stk::StkFloat	in_samplerate = in_file.fileRate();
	stk::StkFrames	filebuff(in_frames, in_file.channels());
	in_file.read(filebuff, 0, false);
	
	// モノラル、float型フォーマットに変換する
	double normalize_gain = 1.0;
	if (in_file.format() == stk::FileRead::STK_SINT8) {
		normalize_gain = 128.0;
	}
	if (in_file.format() == stk::FileRead::STK_SINT16) {
		normalize_gain = 32768.0;
	}
	if (in_file.format() == stk::FileRead::STK_SINT24) {
		normalize_gain = 8388608.0;
	}
	if (in_file.format() == stk::FileRead::STK_SINT32) {
		normalize_gain = 2147483648.0;
	}
	normalize_gain *= in_file.channels();

	stk::StkFrames	floatbuff(in_frames, 1);
	for (uint32_t i=0; i<in_frames; ++i) {
		floatbuff[i] = 0;
		for (unsigned int ch=0; ch<in_file.channels(); ++ch) {
			floatbuff[i] += (filebuff[i*in_file.channels() + ch] / normalize_gain) * wave_gain;
		}
	}
	if (in_file.channels() > 1) {
		std::cout << "convert " << in_file.channels() << " channels to mono." << std::endl;
	}
	
	// 目的のDPCMのサンプリングレートに変換
	double dst_rate = dpcmSampleRateValue(sample_rate_ind);
	int	src_frames = static_cast<int>(in_frames * dst_rate / in_samplerate + 1);
	stk::StkFrames preprocess_buff(src_frames, 1);
	
	if (no_resampling) {
		std::cout << "no resampling." << std::endl;
		src_frames = in_frames;
		preprocess_buff.resize(src_frames);
		for (uint32_t i=0; i<in_frames; ++i) {
			preprocess_buff[i] = floatbuff[i];
		}
		src_frames = in_frames;
	}
	else {
		std::cout << "resample " << in_samplerate << "Hz to " << dst_rate << "Hz" << std::endl;
		src_frames = resample(floatbuff, in_samplerate, preprocess_buff,
							  src_frames, dst_rate, 0.995, resample_window_len);
		preprocess_buff.resize(src_frames);
	}
	
	
	// stkの動作サンプリングレートを変更する
	stk::Stk::setSampleRate(dst_rate);
	
	// 波形の前処理を行う
	preprocessInputBuffer(in_file_name, preprocess_buff, dst_rate,
						  shifter_weight, center_bias_level,
						  use_linearity_correction);
	
	if (output_preprocessed) {
		// 中心補正済み波形を出力
		stk::StkFrames	progress_out_buff((uint32_t)preprocess_buff.size(), 1);
		for (size_t i=0; i<preprocess_buff.size(); ++i) {
			progress_out_buff[i] = preprocess_buff[i];
			if (progress_out_buff[i] > 1.0) progress_out_buff[i] = 1.0;
			if (progress_out_buff[i] < -1.0) progress_out_buff[i] = -1.0;
		}
		outputBufferToFile(in_file_name, progress_out_buff, "_processed", dst_rate);
	}
	
	// DPCM変換の開始
	// サイズは16byte=128sample境界に合うように切り上げ
	int dpcm_frames = ((src_frames + 127) / 128) * 128;
	unsigned char *dpcm = new unsigned char[dpcm_frames / 8];
	int dpcm_initial_volume = 2*static_cast<int>((force_initial_volume >= 0)?(force_initial_volume/2):(preprocess_buff[0] * 64.0));
	if (dpcm_initial_volume < 0) {
		dpcm_initial_volume = 0;
	}
	dpcm_frames = encodeDpcm(preprocess_buff, dpcm_initial_volume,
							 dpcm, dpcm_frames/8, dither_mode, noise_shape_mode,
							 use_linearity_correction);
	
	// 波形最大値を求めるために復号する
	int volume_max = decodeDpcm(in_file_name, dpcm, dpcm_frames, dpcm_initial_volume, center_bias_level, output_encoded_wav, dst_rate);
	
	// dmcファイル出力
	if (outputDpcmToFile(in_file_name, dpcm, dpcm_frames, dpcm_initial_volume)) {
		exit(1);
	}
	
	delete [] dpcm;
	
	std::cout << dpcm_frames/8 << " bytes";
	std::cout << "(" << src_frames << "samples+" << (dpcm_frames - src_frames) << "samples padding)" << std::endl;
	std::cout << "initial volume: " << dpcm_initial_volume << std::endl;
	std::cout << "max volume: " << volume_max << std::endl;
	
	return 0;
}

int main(int argc, char * argv[])
{
	std::list<std::string> in_file_list;

	bool dmc2wav_mode = false;
	
	// コマンドラインオプションを処理する
	for (int i=1; i<argc; ++i) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'g':
					wave_gain = std::stod(&argv[i][2]);
					break;
					
				case 'i':
					force_initial_volume = std::stoi(&argv[i][2]);
					if (force_initial_volume > 127) help_mode = true;
					if (force_initial_volume < 0) help_mode = true;
					break;
					
				case 'c':
					center_bias_level = std::stoi(&argv[i][2]);
					if (center_bias_level > 127) help_mode = true;
					if (center_bias_level < 0) help_mode = true;
					break;

				case 'w':
					shifter_weight = std::stod(&argv[i][2]);
					break;
					
				case 'r':
					sample_rate_ind = std::stoi(&argv[i][2]);
					if (sample_rate_ind > 15 || sample_rate_ind < 0) {
						help_mode = true;
					}
					break;
					
				case 'n':
					no_resampling = true;
					break;
					
				case 'd':
					dither_mode = std::stoi(&argv[i][2]);
					if (dither_mode < 0 || dither_mode > 3) {
						help_mode = true;
					}
					break;
					
				case 's':
					noise_shape_mode = std::stoi(&argv[i][2]);
					if (noise_shape_mode < 0 || noise_shape_mode > 3) {
						help_mode = true;
					}
					break;

				case 'l':
					use_linearity_correction = false;
					break;
					
				case 'e':
					output_envelope = true;
					break;
					
				case 'p':
					output_preprocessed = true;
					break;

				case 'o':
					output_encoded_wav = true;
					break;

				case 'h':
					help_mode = true;
					break;
					
				case 'f':
					dmc2wav_mode = true;
					break;
					
				default:
					std::cout << "Invalid option: " << argv[i] << std::endl;
					help_mode = true;
					break;
			}
		}
		else {
			std::string path(argv[i]);
			in_file_list.push_back(path);
		}
	}
	
	if (help_mode || (in_file_list.size() == 0)) {
		std::cout << "Usage: dpcmc [options] <input> .." << std::endl;
		std::cout << "<input> supports uncompressed WAV, AIFF/AIFC, SND (AU), MAT-file (Matlab)" << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "  -g[gain]   Input volume(default:gain=1.0)" << std::endl;
		std::cout << "  -i[0..127] Initial volume(default:auto)" << std::endl;
		std::cout << "  -c[0..127] Center bias level(default:64)" << std::endl;
		std::cout << "  -w[weight] Slope weight(default:weight=1.0)" << std::endl;
		std::cout << "  -r[rate]   Output sampling rate(default:15) rate:" << std::endl;
		for (int i=0; i<16; ++i) {
			if ((i % 4) == 0) {
				std::cout << "                 ";
			}
			if (i < 10) {
				std::cout << " ";
			}
			std::cout << i << ":" << std::showpoint << dpcmSampleRateValue(i) << "Hz ";
			if ((i % 4) == 3) {
				std::cout << std::endl;
			}
		}
		std::cout << "  -n         No resampling" << std::endl;
		std::cout << "  -d[mode]   Dither mode:" << std::endl;
		std::cout << "                 0 : off" << std::endl;
		std::cout << "                 1 : Rectangle" << std::endl;
		std::cout << "                 2 : Triangle" << std::endl;
		std::cout << "                 3 : Highpassed Triangle(default)" << std::endl;
		std::cout << "  -s[mode]   Noise shaping mode:" << std::endl;
		std::cout << "                 0 : off" << std::endl;
		std::cout << "                 1 : Lowpass" << std::endl;
		std::cout << "                 2 : Highpass" << std::endl;
		std::cout << "                 3 : equal-loudness(default)" << std::endl;
		std::cout << "  -l         Disable correction of linearity(default:enable)" << std::endl;
		std::cout << "  -e         Output envelope to wav file(default:off)" << std::endl;
		std::cout << "  -p         Output preprocessed waveform to wav file(default:off)" << std::endl;
		std::cout << "  -o         Output waveform after encode to wav file(default:off)" << std::endl;
		std::cout << "  -f         Convert dmc to wav file" << std::endl;
		std::cout << "  -h         Show this help" << std::endl;
		exit(0);
	}
	
	if (dmc2wav_mode) {
		if (force_initial_volume < 0) {
			force_initial_volume = 64;
		}
		
		std::for_each(in_file_list.cbegin(), in_file_list.cend(), [](std::string in_file_path) {
			processDmcInputFile(in_file_path);
		});
	}
	else {
		std::cout << "gain: " << std::showpoint << wave_gain << std::endl;
		if (force_initial_volume < 0) {
			std::cout << "initial volume: auto" << std::endl;
		}
		std::cout << "center bias level: " << center_bias_level << std::endl;
		std::cout << "slope weight: " << std::showpoint << shifter_weight << std::endl;
		switch (dither_mode) {
			case 0:
				std::cout << "dither: off" << std::endl;
				break;
			case 1:
				std::cout << "dither: Rectangle" << std::endl;
				break;
			case 2:
				std::cout << "dither: Triangle" << std::endl;
				break;
			case 3:
				std::cout << "dither: Highpass Triangle" << std::endl;
				break;
			default:
				break;
		}
		switch (noise_shape_mode) {
			case 0:
				std::cout << "noise shape: off" << std::endl;
				break;
			case 1:
				std::cout << "noise shape: Lowpass" << std::endl;
				break;
			case 2:
				std::cout << "noise shape: Highpass" << std::endl;
				break;
			case 3:
				std::cout << "noise shape: equal-loudness" << std::endl;
				break;
			default:
				break;
		}
		if (use_linearity_correction) {
			std::cout << "linearity correction: on" << std::endl;
		}
		else {
			std::cout << "linearity correction: off" << std::endl;
		}

		std::for_each(in_file_list.cbegin(), in_file_list.cend(), [](std::string in_file_path) {
			processInputFile(in_file_path);
		});
	}
	
	return 0;
}
