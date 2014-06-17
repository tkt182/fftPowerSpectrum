#ifndef __FFT_H_
#define __FFT_H_

#include <iostream>
#include <math.h>


class fft{

private:

	const float _pi;
	const int   _maxFastBits;

	int **_fftBitTable;


	/**
	 * ある値が2のべき乗かどうかをチェックする
	 */
	int IsPowerOfTwo(int val);


	/**
	 * 値が2の何乗かを返す
	 */
	int NumberOfBitsNeeded(int powerOfTwo);

	/**
     * ビット演算を使ってデータソートの順番を決定する
	 */
	int ReverseBits(int index, int numBits);

	/**
     * データソートの順番を取得する
	 */
	inline int FastReverseBits(int i, int numBits);


	/**
     * 入力データを整形するための窓関数を実行
	 */
	void WindowFunc(int windowFuncId, int numSamples, float *in);


	/**
	 * シャッフリングを行うため、入力データのソート順を計算する
	 */
	void InitFFT();



public:

	fft();
	~fft();


    /**
	 * FFTを実行する
	 */
	void ExecFFT(int numSamples, bool inverseTransform,
		float *inReal, float *inImag, float *outReal, float *outImag);

	/**
	 * 実数専用のFFTを実行する
	 */
	void RealFFT(int numSamples, float *inReal, float *outReal, float *outImag);

	/**
	 * パワースペクトルを計算する
	 */
	void PowerSpectrum(int start, int half, float *data, int windowSize,
		float *magnitude, float *phase, float *power, float *avgPower);


	/**
	 * 逆フーリエ変換を実行する
	 */
	void InversePoserSpectrum(int start, int half, int windowSize, 
		float *finalOut,float *magnitude,float *phase); 

};




#endif /* __FFT_H_ */