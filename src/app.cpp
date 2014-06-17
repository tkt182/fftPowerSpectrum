#include "app.h"

//--------------------------------------------------------------
void App::setup(){

	ofSetVerticalSync(true);
	ofBackground(0, 0, 0);

	
	this->soundStream.listDevices();
	//this->soundStream.setDeviceID(1);   // デバイス番号の設定は、環境によって異なる

	this->soundStream.setup(this, 0, 2, 44100, BUFFER_SIZE, 4);


	this->_left  = new float[BUFFER_SIZE];
	this->_right = new float[BUFFER_SIZE];
	

	for(int i = 0; i < NUM_WINDOWS; i++){
	
		for(int j = 0; j < BUFFER_SIZE / 2; j++){
		
			this->_freq[i][j] = 0.0f;
		
		}

	}


}

//--------------------------------------------------------------
void App::update(){

	

}

//--------------------------------------------------------------
void App::draw(){


	float avg_power = 0.0f;

	// 左チャンネルのデータに対してのみFFTを実行
	this->_fft.PowerSpectrum(0, (int)BUFFER_SIZE/2, 
		this->_left, BUFFER_SIZE, &this->_magnitude[0], &this->_phase[0], &this->_power[0], &avg_power);


	// 振幅を描画(グラフィックスイコライザ)
	// _magnitude[0]は直流成分(時間変化のない=振動しない成分)のため、無視する
	float width = (float)ofGetWidth() / (float)(BUFFER_SIZE/2 - 1);
	for(int i = 1; i < (int)(BUFFER_SIZE/2); i++){
		ofRect((i - 1)*width, ofGetHeight(), width, -(this->_magnitude[i] * 200.0));
	}


}


//--------------------------------------------------------------
void App::exit(){

	this->soundStream.close();
	ofBaseApp::exit();

}

//--------------------------------------------------------------
void App::audioIn(float *input, int bufferSize, int nChannels){


	// ユーザが指定したバッファと、実際にドライバが確保するバッファのサイズを
	// 比較し、小さい方を採用する.もし大きい値を使うと、SoundStreamを閉じる際に
	// エラーが発生する
	int minBufferSize = min(BUFFER_SIZE, bufferSize);

	for(int i = 0; i < minBufferSize; i++){

		this->_left[i]  = input[i * 2];
		this->_right[i] = input[i * 2 + 1];

	}
	

}


//--------------------------------------------------------------
void App::keyPressed(int key){

	// ESC
	if(key == 27){

		delete this->_left;
		delete this->_right;

		this->exit();

	}

}

//--------------------------------------------------------------
void App::keyReleased(int key){

}

//--------------------------------------------------------------
void App::mouseMoved(int x, int y){

}

//--------------------------------------------------------------
void App::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void App::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void App::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void App::windowResized(int w, int h){

}

//--------------------------------------------------------------
void App::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void App::dragEvent(ofDragInfo dragInfo){ 

}