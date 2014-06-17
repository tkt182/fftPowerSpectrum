#ifndef __APP_H_
#define __APP_H_

#include "ofMain.h"
#include "fft.h"


#define BUFFER_SIZE 256
#define NUM_WINDOWS 80


class App : public ofBaseApp{


	public:
		void setup();
		void update();
		void draw();
		void exit();
		
		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y);
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);


		// FFTŠÖ˜A
		void audioIn(float * input, int bufferSize, int nChannels);
	
		float *_left;
		float *_right;

		fft    _fft;

		float  _magnitude[BUFFER_SIZE];
		float  _phase[BUFFER_SIZE];
		float  _power[BUFFER_SIZE];
		float  _freq[NUM_WINDOWS][BUFFER_SIZE/2];


		ofSoundStream soundStream;



};

#endif /* __APP_H_ */
