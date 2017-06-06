
#include <stdlib.h>
#include <math.h>
#include <ti/dsplib/dsplib.h>

/* Global definitions */
//#define _LITTLE_ENDIAN
/* ����� ������ ��� ���� ��������� �������� ���*/
#define N 512

unsigned char brev[64] = {
    0x0, 0x20, 0x10, 0x30, 0x8, 0x28, 0x18, 0x38,
    0x4, 0x24, 0x14, 0x34, 0xc, 0x2c, 0x1c, 0x3c,
    0x2, 0x22, 0x12, 0x32, 0xa, 0x2a, 0x1a, 0x3a,
    0x6, 0x26, 0x16, 0x36, 0xe, 0x2e, 0x1e, 0x3e,
    0x1, 0x21, 0x11, 0x31, 0x9, 0x29, 0x19, 0x39,
    0x5, 0x25, 0x15, 0x35, 0xd, 0x2d, 0x1d, 0x3d,
    0x3, 0x23, 0x13, 0x33, 0xb, 0x2b, 0x1b, 0x3b,
    0x7, 0x27, 0x17, 0x37, 0xf, 0x2f, 0x1f, 0x3f
};

/* ����������� ����� �� ���� ������� �����  8 */
#define AL 8

/* ����� ����� �������� ������� �� ��������� */
#pragma DATA_ALIGN(in_signal,AL) // �����������
float   in_signal [2*N]; // ��� ���������� ��������
// 2*N - ��� REAL i IMAGIN ���������  ������������ �������

/* ����� ����� �������� ������� �� ��������� ��� ��� */
#pragma DATA_ALIGN(in_sig_fft,AL) // �����������
float   in_sig_fft [2*N]; // ��� ���������� ��������
// ���� ��� ��� ����� ��������

#pragma DATA_ALIGN(out_signal,AL) // �����������
float   out_signal [2*N]; // ��� ������� �������

#pragma DATA_ALIGN(w_sp, AL)
float   w_sp [N]; // ��� ������� cos | sin - Twiddle factor

// ������ ��� ����� � ����� ��������� ��������� �������
float y_real_sp [N];
float y_imag_sp [N];
// ����� ��� ������ ����� � ����� ���������
float a4h[N];
float a4h_out[N/2];

// ��������� �-��
void generateInput();
void seperateRealImg (float *real, float * img, float * cmplx, int size_cmplx);
void tw_gen (float *w, int n);
void magnitude(float *in_real, float * in_img, float * out_arr, int size);


int main(void) {
	int i;
	/* ����������� ������� */
	generateInput();

	/* ������� ������ ����� (�������) �������� �� �����*/
	tw_gen(w_sp, N);

	// ��������� ���
	DSPF_sp_fftSPxSP(N, in_sig_fft, w_sp, out_signal, brev, 4, 0, N);

	/* ��������� ������������ ������� �� real � imaginary ��� */
	seperateRealImg (y_real_sp, y_imag_sp, out_signal, N);

	// ������
	magnitude(y_real_sp, y_imag_sp, a4h, N);

	// ���������� ������� �������
	for (i = 1; i < N/2; ++i){
		a4h_out[i] = a4h[i] * 2 / N;
	}
	a4h_out[0] /= N;

	i=0;// for breakpoint

	return 0;
}



/*
    ������� ������ ������ ����� ������� � ������ ����� � �����
*/
void generateInput () {
    int   i;
    float FreqSignal1, sinWaveMag1 ;
    float FreqSignal2, sinWaveMag2 ;
    float real_s, img_s ; // ��� ������� �� ����� ��������

    float FreqSample;

    /* �������  ������������� Hz*/
    FreqSample = 48000;

    /* �������� ������� */
    sinWaveMag1 = 5;
    sinWaveMag2 = 10;

    /* �������  ������� Hz*/
    FreqSignal1 = 2000;
    FreqSignal2 = 5550;

    /* ����������� ������� ������ �� ������ ��������� ������� */
	for (i = 0; i < N; i++) {
		real_s  = (float) (
				sinWaveMag1 * cos(2*3.14*i * FreqSignal1 /FreqSample)+
				sinWaveMag2 * cos(2*3.14*i* FreqSignal2 /FreqSample)
		);

		img_s = (float)0.0;

		in_signal[2*i] = real_s;
		in_signal[2*i + 1]  = img_s;
	}

	// ���� ������ � in_sig_fft
	for (i = 0; i < N; i++) {
		in_sig_fft[2*i] = in_signal[2*i];
		in_sig_fft[2*i + 1]  = in_signal[2*i+1];
	}

}

/*
    ������� ������ ����� (�������) �������� �� �����  ��� ���
*/

void tw_gen (float *w, int n)
{
    int i, j, k;
    const double PI = 3.141592654;

    for (j = 1, k = 0; j <= n >> 2; j = j << 2)
    {
        for (i = 0; i < n >> 2; i += j)
        {
#ifdef _LITTLE_ENDIAN
            w[k]     = (float) sin (2 * PI * i / n);
            w[k + 1] = (float) cos (2 * PI * i / n);
            w[k + 2] = (float) sin (4 * PI * i / n);
            w[k + 3] = (float) cos (4 * PI * i / n);
            w[k + 4] = (float) sin (6 * PI * i / n);
            w[k + 5] = (float) cos (6 * PI * i / n);
#else
            w[k]     = (float)  cos (2 * PI * i / n);
            w[k + 1] = (float) -sin (2 * PI * i / n);
            w[k + 2] = (float)  cos (4 * PI * i / n);
            w[k + 3] = (float) -sin (4 * PI * i / n);
            w[k + 4] = (float)  cos (6 * PI * i / n);
            w[k + 5] = (float) -sin (6 * PI * i / n);
#endif
            k += 6;
        }
    }
}



/*
    ������� seperateRealImg
    ��� ��������� real and imaginary ����� ���� ���
    ���� ������� ��� ����, ��� ���������� �������� ���,
    �������������� CCS graph
*/

void seperateRealImg (float *real, float * img, float * cmplx, int size_cmplx) {
    int i, j;

    for (i = 0, j = 0; j < size_cmplx; i+=2, j++) {
    	real[j] = cmplx[i];
    	img[j] = cmplx[i + 1];
    }
}

/*
 �������  ���������� ������ ������������ �����
 */
void magnitude(float *in_real, float * in_img, float * out_arr, int size){
	int i;
	 for (i = 0; i < size; i++) {
		 out_arr[i] = (float)sqrt(in_real[i]* in_real[i] + in_img[i]*in_img[i]);
	    }
}

// from some paper
/*
 void tw_gen (float *w, int n)
{
    int i, j, k;
    double x_t, y_t, theta1, theta2, theta3;
    const double PI = 3.141592654;

    for (j = 1, k = 0; j <= n >> 2; j = j << 2)
    {
        for (i = 0; i < n >> 2; i += j)
        {
            theta1 = 2 * PI * i / n;
            x_t = cos (theta1);
            y_t = sin (theta1);
            w[k] = (float) x_t;
            w[k + 1] = (float) y_t;

            theta2 = 4 * PI * i / n;
            x_t = cos (theta2);
            y_t = sin (theta2);
            w[k + 2] = (float) x_t;
            w[k + 3] = (float) y_t;

            theta3 = 6 * PI * i / n;
            x_t = cos (theta3);
            y_t = sin (theta3);
            w[k + 4] = (float) x_t;
            w[k + 5] = (float) y_t;
            k += 6;
        }
    }
}
*/
