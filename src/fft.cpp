/**********************************************************************

  fft.cpp

  
  This class is a C++ wrapper for original code 
  written by Dominic Mazzoni in September 2000

  This file contains a few FFT routines, including a real-FFT
  routine that is almost twice as fast as a normal complex FFT,
  and a power spectrum routine which is more convenient when
  you know you don't care about phase information.  It now also
  contains a few basic windowing functions.

  Some of this code was based on a free implementation of an FFT
  by Don Cross, available on the web at:

    http://www.intersrv.com/~dcross/fft.html

  The basic algorithm for his code was based on Numerical Recipes
  in Fortran.  I optimized his code further by reducing array
  accesses, caching the bit reversal table, and eliminating
  float-to-double conversions, and I added the routines to
  calculate a real FFT and a real power spectrum.

  Note: all of these routines use single-precision floats.
  I have found that in practice, floats work well until you
  get above 8192 samples.  If you need to do a larger FFT,
  you need to use doubles.

**********************************************************************/

#include "fft.h"	


/**
 * コンストラクタ
 */
fft::fft() : _pi(3.141592654), _maxFastBits(16){

	this->_fftBitTable = NULL;

}

/**
 * デストラクタ
 */
fft::~fft(){
}


/**
 * ある値が2のべき乗かどうかをチェックする
 */
int fft::IsPowerOfTwo(int val){

	if(val < 2){
		return false;
	}

	if(val & (val - 1)){
		return false;
	}

	return true;
}

/**
 * 値が2の何乗かを返す
 */
int fft::NumberOfBitsNeeded(int powerOfTwo){
	
	if(powerOfTwo < 2){
		std::cerr << "Error: FFT called with size " << powerOfTwo << std::endl;
		exit(1);
	}

	for(int i = 0;; i++){

		if(powerOfTwo & (1 << i)){
			return i;
		}

	}
	
}

/**
 * ビット演算を使ってデータソートの順番を決定する
 */
int fft::ReverseBits(int index, int numBits){

	int i, rev;

	for(i = rev = 0; i < numBits; i++){
	
		rev = (rev << 1) | (index & 1);
		index >>= 1;
	
	}

	return rev;

}

/**
 * データソートの順番を取得する
 */
inline int fft::FastReverseBits(int i, int numBits){

	if(numBits <= this->_maxFastBits){

		return this->_fftBitTable[numBits - 1][i];
	
	}else{
	
		return this->ReverseBits(i, numBits);

	}

}


/**
 * 入力データを整形するための窓関数を実行
 */
void fft::WindowFunc(int windowFuncId, int numSamples, float *in){

	if(windowFuncId == 1){
	
		// Bartlett(triangular) window
		for(int i = 0; i < numSamples / 2; i ++){
			in[i]                    *= (i / (float)(numSamples / 2));
			in[i + (numSamples / 2)] *= (1.0 - (i / (float)(numSamples / 2)));
		}

	}

	if(windowFuncId == 2){
	
		// Hamming
		for(int i = 0; i < numSamples; i++){
			in[i] *= 0.54 - 0.46 * cos(2.0 * this->_pi * (float)i / (float)(numSamples - 1));
		}
	
	}

	if(windowFuncId == 3){
	
		// Hanning
		for(int i = 0; i < numSamples; i++){
			in[i] *= 0.5 - 0.5 * cos(2.0 * this->_pi * (float)i / (float)(numSamples - 1));
		}
	
	}


}



/**
 * シャッフリングを行うため、入力データのソート順を計算する
 */
void fft::InitFFT(){

	this->_fftBitTable = new int *[this->_maxFastBits];

	int len = 2;
	for(int b = 1; b <= this->_maxFastBits; b++){
	
		this->_fftBitTable[b - 1] = new int[len];

		for(int i = 0; i < len; i++){
			this->_fftBitTable[b - 1][i] = this->ReverseBits(i, b); 
		}

		len <<= 1;
	
	}

}


/**
 * FFTを実行する
 */
void fft::ExecFFT(int numSamples, bool inverseTransform,
	float *inReal, float *inImag, float *outReal, float *outImag){

	int numBits;   // インデックスを格納するのに必要なビット数
	int blockSize, blockEnd;

	double angleNumerator = 2.0 * this->_pi; // 角度の分母
	float tr, ti;                           // temp Real, temp Imaginary

	if(!IsPowerOfTwo(numSamples)){

		std::cerr << numSamples << "is not a power of two " << std::endl;
		exit(1);
	
	}

	// バタフライ演算を行うための、入力データの並び順を計算
	if(!this->_fftBitTable){
		this->InitFFT();
	}


	// 逆フリーエ変換の場合は角度の分母をマイナスにする
	if(inverseTransform){
		angleNumerator = -angleNumerator;
	}

	numBits = this->NumberOfBitsNeeded(numSamples);

	// 入力データの並び変え(シャッフリング)
	for(int i = 0; i < numSamples; i++){

		int j      = this->FastReverseBits(i, numBits);
		outReal[j] = inReal[i];
		outImag[j] = (inImag == NULL) ? 0.0 : inImag[i];
	}

	// FFTを実行する
	blockEnd = 1;


	//////// バタフライ演算 ///////

	// バタフライ演算のN階層(2^Nの"N"が階層数となる)分のループを回す
	// blockSizeはバタフライ演算の対象となるブロックの要素数
	// ※ 1階層目 ⇒ バタフライ演算は2個の要素間で行う
	//    2階層目 ⇒ バタフライ演算は4個の要素間で行う
	//    3階層目 ⇒ バタフライ演算は8個の要素間で行う
	for(blockSize = 2; blockSize <= numSamples; blockSize <<= 1){
		

		// 波長(deltaAngle)を短くしていく(π, π/2, π/4, …)
		double deltaAngle = angleNumerator / (double) blockSize;

		float sm2 = sin(-2.0 * deltaAngle);
		float sm1 = sin(-deltaAngle);
		float cm2 = cos(-2.0 * deltaAngle);
		float cm1 = cos(-deltaAngle);
		float w   = 2.0 * cm1;              // 回転因子(波長を変えていくイメージ)
		
		float ar0, ar1, ar2, ai0, ai1, ai2;

		// ブロック数分のループ
		for(int i = 0; i < numSamples; i += blockSize){
		
			ar2 = cm2;
			ar1 = cm1;

			ai2 = sm2;
			ai1 = sm1;

			// ブロック内の要素数分のループ
			for(int j = i, n = 0; n < blockEnd; j++, n++){
			

				// 実部(cos)の計算
				ar0 = w * ar1 - ar2;
				ar2 = ar1;
				ar1 = ar0;

				// 虚部(sin)の計算
				ai0 = w * ai1 - ai2;
				ai2 = ai1;
				ai1 = ai0;

				// 時間間引きアルゴリズム
				// j番目とk番目の要素でバタフライ演算
				int k = j + blockEnd;

				tr = ar0 * outReal[k] - ai0 * outImag[k];
				ti = ar0 * outImag[k] + ai0 * outReal[k];
			
				outReal[k] = outReal[j] - tr;
				outImag[k] = outImag[j] - ti;

				outReal[j] += tr;
				outImag[j] += ti;

			}
		
		}

		blockEnd = blockSize;
	
	}

	// 逆フーリエ変換の場合、正規化する
	if(inverseTransform){
	
		float denom = (float)numSamples;

		for(int i = 0; i < numSamples; i++){
			outReal[i] /= denom;
			outImag[i] /= denom;
		}

	}


}


/**
 * 実数専用のFFTを実行する
 */
void fft::RealFFT(int numSamples, float *inReal, float *outReal, float *outImag){

	int   half  = numSamples / 2;
	float theta = this->_pi / half;

	float *tmpReal = new float[half];
	float *tmpImag = new float[half];

	// 計算量を減らすため、奇数番目の実数データを虚数配列に保存する
	for(int i = 0; i < half; i++){
	
		tmpReal[i] = inReal[2 * i];
		tmpImag[i] = inReal[2 * i + 1];

	}

	// FFTを実行する.
	// 入力データは半分ずつに分かれてtmpRealとtmpImagに保存されているので、
	// データ数はhalf(numSamples/2)となる.
	// 出力結果のoutReal, outImagはFFT後の周波数領域(横軸が周波数)のデータになっているが、
	// それぞれhalf(numSamples/2)番目の要素までしか値が埋まっていない.
	this->ExecFFT(half, 0, tmpReal, tmpImag, outReal, outImag);

	float wTmp = float(sin(0.5 * theta));  // 回転因子

	float wpr = -2.0 * wTmp * wTmp;
	float wpi = float(sin(theta));
	float wr  = 1.0 + wpr;
	float wi  = wpi;

	int i3;
	float h1r, h1i, h2r, h2i;

	// バタフライ演算
	for(int i = 1; i < half / 2; i++){
	
		i3 = half - i;

		// バタフライ演算の入力データを修正.
		// ポイントは、outImagのデータも本質的に実数値だというとこ.
		// 極座標の虚数部と実数部を入れかれるには、90°回転させる.
		// そのため、h2rにはoutImag、h2iにはoutRealの値を代入する.
		// ※ 0.5を掛けているのは、2要素分を足しているため.
		h1r =  0.5 * (outReal[i] + outReal[i3]);
		h1i =  0.5 * (outImag[i] - outImag[i3]);
		h2r =  0.5 * (outImag[i] + outImag[i3]);
		h2i = -0.5 * (outReal[i] - outReal[i3]);

		outReal[i]  =  h1r + wr * h2r - wi * h2i;
		outImag[i]  =  h1i + wr * h2i + wi * h2r;
		outReal[i3] =  h1r - wr * h2r + wi * h2i;
		outImag[i3] = -h1r + wr * h2i + wi * h2r;

		wr = (wTmp = wr) * wpr - wi * wpi + wr;
		wi = wi * wpr + wTmp * wpi + wi;

	}

	// 0はナイキスト周波数になるため、ループ内で計算できない.
	// ※ そもそもN/2が存在しない.
	outReal[0] = (h1r = outReal[0]) + outImag[0];
	outImag[0] = h1r - outImag[0];

	delete[] tmpReal;
	delete[] tmpImag;


}


/**
 * パワースペクトルを計算する
 * 
 * @param[in] start 入力サンプルの開始位置
 * @param[in] half 使用する入力サンプル数の半分
 * @param[in] data 入力サンプル
 * @param[in] windowSize 入力サンプル数(周波数分解数)
 * @param[in,out] magnitude 振幅
 * @param[in,out] phase 位相
 * @param[in,out] power パワースペクトル(振幅の2乗[dB])
 * @param[in,out] ave_power 全周波数のパワーの平均
 *
 */
void fft::PowerSpectrum(int start, int half, float *data, int windowSize,
	float *magnitude, float *phase, float *power, float *avgPower){


	int windowFuncId = 1;    // WindowFunc ID
	float totalPower = 0.0f;

	float *inReal  = new float[windowSize];  // 入力の実数部
	float *inImag  = new float[windowSize];  // 入力の虚数部
	float *outReal = new float[windowSize];  // 出力の実数部
	float *outImag = new float[windowSize];  // 出力の虚数部

	for(int i = 0; i < windowSize; i++){
		inReal[i] = data[start + i];
	}

	// 入力データに窓関数を適用する
	this->WindowFunc(windowFuncId, windowSize, inReal);

	// 実数専用のFFTを実行
	this->RealFFT(windowSize, inReal, outReal, outImag);

	for(int i = 0; i < half; i++){

		power[i]    = outReal[i] * outReal[i] + outImag[i] * outImag[i]; // パワースペクトル
		totalPower += power[i];

		//magnitude[i] = 2.0 * sqrt(power[i]);
		magnitude[i] = sqrt(power[i]);                    // 振幅
		phase[i]     = atan2(outImag[i], outReal[i]);     // 位相

	}

	// パワースペクトルの平均値を計算する
	*(avgPower) = totalPower / (float)half;

	delete[] inReal;
	delete[] inImag;
	delete[] outReal;
	delete[] outImag;

}


/**
 * 逆フーリエ変換を実行する
 */
void fft::InversePoserSpectrum(int start, int half, int windowSize, 
	float *finalOut,float *magnitude,float *phase){

	int windowFuncId = 1;    // WindowFunc ID

	float *inReal  = new float[windowSize];  // 入力の実数部
	float *inImag  = new float[windowSize];  // 入力の虚数部
	float *outReal = new float[windowSize];  // 出力の実数部
	float *outImag = new float[windowSize];  // 出力の虚数部

	for(int i = 0; i < half; i++){
	
		inReal[i] = magnitude[i] * cos(phase[i]);
		inImag[i] = magnitude[i] * sin(phase[i]);
	
	}

	for(int i = half; i < windowSize; i++){
	
		inReal[i] = 0.0;
		inImag[i] = 0.0;
	
	}

	// 逆フーリエ変換を実行
	this->ExecFFT(windowSize, 1, inReal, inImag, outReal, outImag);
	this->WindowFunc(windowFuncId, windowSize, outReal);

	for(int i = 0; i < windowSize; i++){
		finalOut[start + i] += outReal[i];
	}

	delete[] inReal;
	delete[] inImag;
	delete[] outReal;
	delete[] outImag;

}

