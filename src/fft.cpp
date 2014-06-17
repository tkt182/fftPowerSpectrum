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
 * �R���X�g���N�^
 */
fft::fft() : _pi(3.141592654), _maxFastBits(16){

	this->_fftBitTable = NULL;

}

/**
 * �f�X�g���N�^
 */
fft::~fft(){
}


/**
 * ����l��2�ׂ̂��悩�ǂ������`�F�b�N����
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
 * �l��2�̉��悩��Ԃ�
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
 * �r�b�g���Z���g���ăf�[�^�\�[�g�̏��Ԃ����肷��
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
 * �f�[�^�\�[�g�̏��Ԃ��擾����
 */
inline int fft::FastReverseBits(int i, int numBits){

	if(numBits <= this->_maxFastBits){

		return this->_fftBitTable[numBits - 1][i];
	
	}else{
	
		return this->ReverseBits(i, numBits);

	}

}


/**
 * ���̓f�[�^�𐮌`���邽�߂̑��֐������s
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
 * �V���b�t�����O���s�����߁A���̓f�[�^�̃\�[�g�����v�Z����
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
 * FFT�����s����
 */
void fft::ExecFFT(int numSamples, bool inverseTransform,
	float *inReal, float *inImag, float *outReal, float *outImag){

	int numBits;   // �C���f�b�N�X���i�[����̂ɕK�v�ȃr�b�g��
	int blockSize, blockEnd;

	double angleNumerator = 2.0 * this->_pi; // �p�x�̕���
	float tr, ti;                           // temp Real, temp Imaginary

	if(!IsPowerOfTwo(numSamples)){

		std::cerr << numSamples << "is not a power of two " << std::endl;
		exit(1);
	
	}

	// �o�^�t���C���Z���s�����߂́A���̓f�[�^�̕��я����v�Z
	if(!this->_fftBitTable){
		this->InitFFT();
	}


	// �t�t���[�G�ϊ��̏ꍇ�͊p�x�̕�����}�C�i�X�ɂ���
	if(inverseTransform){
		angleNumerator = -angleNumerator;
	}

	numBits = this->NumberOfBitsNeeded(numSamples);

	// ���̓f�[�^�̕��ѕς�(�V���b�t�����O)
	for(int i = 0; i < numSamples; i++){

		int j      = this->FastReverseBits(i, numBits);
		outReal[j] = inReal[i];
		outImag[j] = (inImag == NULL) ? 0.0 : inImag[i];
	}

	// FFT�����s����
	blockEnd = 1;


	//////// �o�^�t���C���Z ///////

	// �o�^�t���C���Z��N�K�w(2^N��"N"���K�w���ƂȂ�)���̃��[�v����
	// blockSize�̓o�^�t���C���Z�̑ΏۂƂȂ�u���b�N�̗v�f��
	// �� 1�K�w�� �� �o�^�t���C���Z��2�̗v�f�Ԃōs��
	//    2�K�w�� �� �o�^�t���C���Z��4�̗v�f�Ԃōs��
	//    3�K�w�� �� �o�^�t���C���Z��8�̗v�f�Ԃōs��
	for(blockSize = 2; blockSize <= numSamples; blockSize <<= 1){
		

		// �g��(deltaAngle)��Z�����Ă���(��, ��/2, ��/4, �c)
		double deltaAngle = angleNumerator / (double) blockSize;

		float sm2 = sin(-2.0 * deltaAngle);
		float sm1 = sin(-deltaAngle);
		float cm2 = cos(-2.0 * deltaAngle);
		float cm1 = cos(-deltaAngle);
		float w   = 2.0 * cm1;              // ��]���q(�g����ς��Ă����C���[�W)
		
		float ar0, ar1, ar2, ai0, ai1, ai2;

		// �u���b�N�����̃��[�v
		for(int i = 0; i < numSamples; i += blockSize){
		
			ar2 = cm2;
			ar1 = cm1;

			ai2 = sm2;
			ai1 = sm1;

			// �u���b�N���̗v�f�����̃��[�v
			for(int j = i, n = 0; n < blockEnd; j++, n++){
			

				// ����(cos)�̌v�Z
				ar0 = w * ar1 - ar2;
				ar2 = ar1;
				ar1 = ar0;

				// ����(sin)�̌v�Z
				ai0 = w * ai1 - ai2;
				ai2 = ai1;
				ai1 = ai0;

				// ���ԊԈ����A���S���Y��
				// j�Ԗڂ�k�Ԗڂ̗v�f�Ńo�^�t���C���Z
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

	// �t�t�[���G�ϊ��̏ꍇ�A���K������
	if(inverseTransform){
	
		float denom = (float)numSamples;

		for(int i = 0; i < numSamples; i++){
			outReal[i] /= denom;
			outImag[i] /= denom;
		}

	}


}


/**
 * ������p��FFT�����s����
 */
void fft::RealFFT(int numSamples, float *inReal, float *outReal, float *outImag){

	int   half  = numSamples / 2;
	float theta = this->_pi / half;

	float *tmpReal = new float[half];
	float *tmpImag = new float[half];

	// �v�Z�ʂ����炷���߁A��Ԗڂ̎����f�[�^�������z��ɕۑ�����
	for(int i = 0; i < half; i++){
	
		tmpReal[i] = inReal[2 * i];
		tmpImag[i] = inReal[2 * i + 1];

	}

	// FFT�����s����.
	// ���̓f�[�^�͔������ɕ������tmpReal��tmpImag�ɕۑ�����Ă���̂ŁA
	// �f�[�^����half(numSamples/2)�ƂȂ�.
	// �o�͌��ʂ�outReal, outImag��FFT��̎��g���̈�(���������g��)�̃f�[�^�ɂȂ��Ă��邪�A
	// ���ꂼ��half(numSamples/2)�Ԗڂ̗v�f�܂ł����l�����܂��Ă��Ȃ�.
	this->ExecFFT(half, 0, tmpReal, tmpImag, outReal, outImag);

	float wTmp = float(sin(0.5 * theta));  // ��]���q

	float wpr = -2.0 * wTmp * wTmp;
	float wpi = float(sin(theta));
	float wr  = 1.0 + wpr;
	float wi  = wpi;

	int i3;
	float h1r, h1i, h2r, h2i;

	// �o�^�t���C���Z
	for(int i = 1; i < half / 2; i++){
	
		i3 = half - i;

		// �o�^�t���C���Z�̓��̓f�[�^���C��.
		// �|�C���g�́AoutImag�̃f�[�^���{���I�Ɏ����l���Ƃ����Ƃ�.
		// �ɍ��W�̋������Ǝ���������ꂩ���ɂ́A90����]������.
		// ���̂��߁Ah2r�ɂ�outImag�Ah2i�ɂ�outReal�̒l��������.
		// �� 0.5���|���Ă���̂́A2�v�f���𑫂��Ă��邽��.
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

	// 0�̓i�C�L�X�g���g���ɂȂ邽�߁A���[�v���Ōv�Z�ł��Ȃ�.
	// �� ��������N/2�����݂��Ȃ�.
	outReal[0] = (h1r = outReal[0]) + outImag[0];
	outImag[0] = h1r - outImag[0];

	delete[] tmpReal;
	delete[] tmpImag;


}


/**
 * �p���[�X�y�N�g�����v�Z����
 * 
 * @param[in] start ���̓T���v���̊J�n�ʒu
 * @param[in] half �g�p������̓T���v�����̔���
 * @param[in] data ���̓T���v��
 * @param[in] windowSize ���̓T���v����(���g������)
 * @param[in,out] magnitude �U��
 * @param[in,out] phase �ʑ�
 * @param[in,out] power �p���[�X�y�N�g��(�U����2��[dB])
 * @param[in,out] ave_power �S���g���̃p���[�̕���
 *
 */
void fft::PowerSpectrum(int start, int half, float *data, int windowSize,
	float *magnitude, float *phase, float *power, float *avgPower){


	int windowFuncId = 1;    // WindowFunc ID
	float totalPower = 0.0f;

	float *inReal  = new float[windowSize];  // ���͂̎�����
	float *inImag  = new float[windowSize];  // ���͂̋�����
	float *outReal = new float[windowSize];  // �o�͂̎�����
	float *outImag = new float[windowSize];  // �o�͂̋�����

	for(int i = 0; i < windowSize; i++){
		inReal[i] = data[start + i];
	}

	// ���̓f�[�^�ɑ��֐���K�p����
	this->WindowFunc(windowFuncId, windowSize, inReal);

	// ������p��FFT�����s
	this->RealFFT(windowSize, inReal, outReal, outImag);

	for(int i = 0; i < half; i++){

		power[i]    = outReal[i] * outReal[i] + outImag[i] * outImag[i]; // �p���[�X�y�N�g��
		totalPower += power[i];

		//magnitude[i] = 2.0 * sqrt(power[i]);
		magnitude[i] = sqrt(power[i]);                    // �U��
		phase[i]     = atan2(outImag[i], outReal[i]);     // �ʑ�

	}

	// �p���[�X�y�N�g���̕��ϒl���v�Z����
	*(avgPower) = totalPower / (float)half;

	delete[] inReal;
	delete[] inImag;
	delete[] outReal;
	delete[] outImag;

}


/**
 * �t�t�[���G�ϊ������s����
 */
void fft::InversePoserSpectrum(int start, int half, int windowSize, 
	float *finalOut,float *magnitude,float *phase){

	int windowFuncId = 1;    // WindowFunc ID

	float *inReal  = new float[windowSize];  // ���͂̎�����
	float *inImag  = new float[windowSize];  // ���͂̋�����
	float *outReal = new float[windowSize];  // �o�͂̎�����
	float *outImag = new float[windowSize];  // �o�͂̋�����

	for(int i = 0; i < half; i++){
	
		inReal[i] = magnitude[i] * cos(phase[i]);
		inImag[i] = magnitude[i] * sin(phase[i]);
	
	}

	for(int i = half; i < windowSize; i++){
	
		inReal[i] = 0.0;
		inImag[i] = 0.0;
	
	}

	// �t�t�[���G�ϊ������s
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

