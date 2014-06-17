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
	 * ����l��2�ׂ̂��悩�ǂ������`�F�b�N����
	 */
	int IsPowerOfTwo(int val);


	/**
	 * �l��2�̉��悩��Ԃ�
	 */
	int NumberOfBitsNeeded(int powerOfTwo);

	/**
     * �r�b�g���Z���g���ăf�[�^�\�[�g�̏��Ԃ����肷��
	 */
	int ReverseBits(int index, int numBits);

	/**
     * �f�[�^�\�[�g�̏��Ԃ��擾����
	 */
	inline int FastReverseBits(int i, int numBits);


	/**
     * ���̓f�[�^�𐮌`���邽�߂̑��֐������s
	 */
	void WindowFunc(int windowFuncId, int numSamples, float *in);


	/**
	 * �V���b�t�����O���s�����߁A���̓f�[�^�̃\�[�g�����v�Z����
	 */
	void InitFFT();



public:

	fft();
	~fft();


    /**
	 * FFT�����s����
	 */
	void ExecFFT(int numSamples, bool inverseTransform,
		float *inReal, float *inImag, float *outReal, float *outImag);

	/**
	 * ������p��FFT�����s����
	 */
	void RealFFT(int numSamples, float *inReal, float *outReal, float *outImag);

	/**
	 * �p���[�X�y�N�g�����v�Z����
	 */
	void PowerSpectrum(int start, int half, float *data, int windowSize,
		float *magnitude, float *phase, float *power, float *avgPower);


	/**
	 * �t�t�[���G�ϊ������s����
	 */
	void InversePoserSpectrum(int start, int half, int windowSize, 
		float *finalOut,float *magnitude,float *phase); 

};




#endif /* __FFT_H_ */