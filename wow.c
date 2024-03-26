#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <immintrin.h>
#include <time.h>

unsigned char p[32], q[32], r[32], inv_r[32], inv_p[32], kkk[32] = {0},out[32];

// x^7+c
static const uint8_t s_box[256] = {198, 199, 70, 57, 92, 213, 163, 101, 243, 107, 133, 248, 150, 155, 108, 181, 3, 216, 45, 145, 16, 53, 244, 135, 83, 192, 137, 230, 128, 240, 99, 222, 167, 22, 95, 59, 77, 30, 254, 73, 228, 161, 249, 242, 174, 20, 11, 206, 50, 55, 235, 34, 74, 219, 109, 1, 166, 134, 152, 239, 65, 232, 114, 189, 160, 97, 201, 149, 104, 229, 184, 153, 171, 113, 165, 121, 217, 82, 157, 33, 118, 39, 141, 88, 116, 31, 131, 93, 76, 215, 210, 86, 203, 119, 170, 56, 84, 245, 226, 117, 183, 48, 140, 12, 6, 209, 196, 190, 0, 236, 188, 204, 32, 127, 139, 251, 18, 29, 129, 144, 241, 9, 44, 96, 25, 172, 15, 63, 13, 225, 90, 194, 7, 125, 200, 247, 182, 193, 246, 110, 185, 69, 146, 176, 250, 75, 130, 94, 187, 234, 238, 227, 223, 36, 178, 191, 164, 162, 23, 214, 47, 173, 58, 103, 124, 2, 197, 85, 52, 64, 37, 21, 61, 168, 115, 147, 43, 154, 158, 233, 40, 19, 132, 72, 28, 138, 175, 100, 122, 123, 35, 78, 159, 156, 46, 186, 91, 10, 180, 67, 120, 98, 79, 89, 252, 169, 102, 14, 17, 5, 179, 41, 221, 237, 148, 27, 60, 224, 26, 211, 143, 106, 177, 112, 151, 42, 195, 66, 81, 212, 111, 231, 255, 54, 62, 80, 38, 202, 126, 4, 24, 220, 208, 49, 205, 71, 218, 68, 8, 105, 87, 136, 253, 207, 142, 51};
static const uint8_t inv_s_box[256] = {108, 55, 165, 16, 239, 209, 104, 132, 248, 121, 197, 46, 103, 128, 207, 126, 20, 208, 116, 181, 45, 171, 33, 158, 240, 124, 218, 215, 184, 117, 37, 85, 112, 79, 51, 190, 153, 170, 236, 81, 180, 211, 225, 176, 122, 18, 194, 160, 101, 243, 48, 255, 168, 21, 233, 49, 95, 3, 162, 35, 216, 172, 234, 127, 169, 60, 227, 199, 247, 141, 2, 245, 183, 39, 52, 145, 88, 36, 191, 202, 235, 228, 77, 24, 96, 167, 91, 250, 83, 203, 130, 196, 4, 87, 147, 34, 123, 65, 201, 30, 187, 7, 206, 163, 68, 249, 221, 9, 14, 54, 139, 230, 223, 73, 62, 174, 84, 99, 80, 93, 200, 75, 188, 189, 164, 133, 238, 113, 28, 118, 146, 86, 182, 10, 57, 23, 251, 26, 185, 114, 102, 82, 254, 220, 119, 19, 142, 175, 214, 67, 12, 224, 58, 71, 177, 13, 193, 78, 178, 192, 64, 41, 157, 6, 156, 74, 56, 32, 173, 205, 94, 72, 125, 161, 44, 186, 143, 222, 154, 210, 198, 15, 136, 100, 70, 140, 195, 148, 110, 63, 107, 155, 25, 137, 131, 226, 106, 166, 0, 1, 134, 66, 237, 92, 111, 244, 47, 253, 242, 105, 90, 219, 229, 5, 159, 89, 17, 76, 246, 53, 241, 212, 31, 152, 217, 129, 98, 151, 40, 69, 27, 231, 61, 179, 149, 50, 109, 213, 150, 59, 29, 120, 43, 8, 22, 97, 138, 135, 11, 42, 144, 115, 204, 252, 38, 232};

int Round;

uint64_t seki(uint64_t a, uint64_t b)
{
	uint64_t c = 0;

	while (a != 0)
	{
		// printf("b %b %b %b\n",a,b,c);
		if ((a & 1) == 1)
			c ^= b;

		b <<= 1;

		a >>= 1;
	}

	return c;
}

uint16_t inv(uint16_t a, uint16_t n)
{
	uint16_t d = n;
	uint16_t x = 0;
	uint16_t s = 1, q = 0, r = 0, t = 0;
	while (a != 0)
	{
		q = d / a;
		r = d % a;
		d = a;
		a = r;
		t = (n + x - q * s);
		x = s;
		s = t;
	}
	// gcd = d  # $\gcd(a, n)$

	return ((x + n) % (n / d)) % n;
}

uint16_t cheki(uint16_t a, uint16_t n)
{
	int i;

	if (a == 0)
		return 0;
	for (i = 0; i < n; i++)
		if (a * i % n == 1)
			return i;
	printf("baka\n");
	exit(1);
}

/*
 * Multiplication in GF(2^8)
 * http://en.wikipedia.org/wiki/Finite_field_arithmetic
 * Irreducible polynomial m(x) = x8 + x4 + x3 + x + 1
 *
 * NOTE: This function can be easily replaced with a look up table for a speed
 *       boost, at the expense of an increase in memory size (around 65 KB). See
 *       the aes.h header file to find the macro definition.
 *
 */
uint8_t gmult(uint8_t a, uint8_t b)
{

	uint8_t p = 0, i = 0, hbs = 0;

	for (i = 0; i < 8; i++)
	{
		if (b & 1)
		{
			p ^= a;
		}

		hbs = a & 0x80;
		a <<= 1;
		if (hbs)
			a ^= 0x1b; // 0000 0001 0001 1011
		b >>= 1;
	}

	return (uint8_t)p;
}

void ha(__uint128_t a, __uint128_t b)
{
	__uint128_t p = {0}, hbs = {0}, out = 1;

	for (int i = 0; i < 8; i++)
	{
		if (b & out)
		{
			p ^= a;
		}

		hbs = a & 0x80;
		a <<= 1;
		if (hbs)
			a ^= 0xf3; // 0000 0000 1111 0011
		b >>= 1;
	}
}

uint8_t rotl(uint8_t x, uint8_t r)
{
	if (r == 0)
		return x;
	if (r < 0)
		return (x >> r) | (x << (8 - r));

	return (x << r) | (x >> (8 - r));
}

void rounder()
{
	for (int i = 0; i < 32; i++)
		q[i] = p[r[inv_p[i]]];
	memcpy(r, q, 32);
	for (int i = 0; i < 32; i++)
		inv_r[r[i]] = i;
}

void reverse()
{
	for (int i = 0; i < 32; i++)
		q[i] = inv_p[r[p[i]]];
	memcpy(r, q, 32);
	for (int i = 0; i < 32; i++)
		inv_r[r[i]] = i;
}

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(short))
#define SWAP(type, a, b) \
	{                    \
		type work = a;   \
		a = b;           \
		b = work;        \
	}

/*
	Fisher-Yates shuffle による方法
	配列の要素をランダムシャッフルする
*/
void random_shuffle(unsigned char *array, size_t size)
{
	for (size_t i = size; i > 1; --i)
	{
		size_t a = i - 1;
		size_t b = rand() % i;
		SWAP(int, array[a], array[b]);
	}
}

int mlt(int x, int y)
{

	if (x == 0 || y == 0)
		return 0;

	return ((x + y - 2) % (256 - 1)) + 1;
}

uint8_t der[4][4] ={

	{2, 3, 4, 5},
	{4, 5, 16, 17},
	{8, 15, 64, 85},
	{16, 17, 29, 28}
	};
//{1,2,3,4},
//{1,4,5,16},
//{1,8,15,64},
//{1,16,17,29}

uint8_t snoot[4][4] = {
	{104, 115, 192, 96},
	{210, 233, 192, 64},
	{26, 80, 192, 48},
	{58, 139, 192, 203}};
	

/*
 * Multiplication of 4 byte words
 * m(x) = x4+1
 */
void coef_mult(uint8_t *a, uint8_t *b, uint8_t *d) {

	d[0] = gmult(a[0],b[0])^gmult(a[0],b[1])^gmult(a[0],b[2])^gmult(a[0],b[3]);
	d[1] = gmult(a[1],b[0])^gmult(a[1],b[1])^gmult(a[1],b[2])^gmult(a[1],b[3]);
	d[2] = gmult(a[2],b[0])^gmult(a[2],b[1])^gmult(a[2],b[2])^gmult(a[2],b[3]);
	d[3] = gmult(a[3],b[0])^gmult(a[3],b[1])^gmult(a[3],b[2])^gmult(a[3],b[3]);
}


/*
 * Transformation in the Cipher that processes the State by cyclically
 * shifting the last three rows of the State by different offsets.
 */
void shift_rows(uint8_t *state)
{

	uint8_t i, k, s, tmp;

	for (i = 1; i < 4; i++)
	{
		// shift(1,4)=1; shift(2,4)=2; shift(3,4)=3
		// shift(r, 4) = r;
		s = 0;
		while (s < i)
		{
			tmp = state[8 * i + 0];

			for (k = 1; k < 8; k++)
			{
				state[8 * i + k - 1] = state[8 * i + k];
			}

			state[8 * i + 8 - 1] = tmp;
			s++;
		}
	}
}

/*
 * Transformation in the Inverse Cipher that is the inverse of
 * ShiftRows().
 */
void inv_shift_rows(uint8_t *state)
{

	uint8_t i, k, s, tmp;

	for (i = 1; i < 4; i++)
	{
		s = 0;
		while (s < i)
		{
			tmp = state[8 * i + 8 - 1];

			for (k = 8 - 1; k > 0; k--)
			{
				state[8 * i + k] = state[8 * i + k - 1];
			}

			state[8 * i + 0] = tmp;
			s++;
		}
	}
}

void add(uint8_t *m, uint8_t *k)
{
	int i, n;
	for (i = 0; i < 32; i++)
	{
		n = (m[i] + k[r[i]]) % 256;
		// m[i]=s_box[((n % 16) + (n >> 4) * 16)];
		m[i] = n;
	}
}

void sub(uint8_t *c, uint8_t *k)
{
	int i;
	for (i = 0; i < 32; i++)
	{
		c[i] = c[i];
		// c[i]=inv_s_box[((n % 16) + (n >> 4) * 16)];
		c[i] = (256 + c[i] - k[r[i]]) % 256;
	}
}

void perm(uint8_t *m, uint8_t *r)
{
	uint8_t u[32];
	for (int i = 0; i < 32; i++)
		u[i] = m[r[i]];
	memcpy(m, u, 32);
}

unsigned char mkbox(unsigned char x)
{
	// GF(256)
	// return gmult(gmult(gmult(x, x), gmult(x, x)), gmult(gmult(x, x), x)) ^ 198;
	// GF(257)
	return x * x * x % 257;
}

uint16_t invb(uint16_t i)
{
	return inv(i * i * i % 257, 257);
}

uint64_t pd(uint64_t a, uint64_t b, uint64_t d)
{
	uint64_t c = 0, hbs = 0, ll = (1 << (8 - 1)), l = d ^ (1 << 8);

	while (a != 0)
	{
		// printf("b %b %b %b\n",a,b,c);
		if ((a & 1) == 1)
			c ^= b;
		hbs = b & (ll);
		b <<= 1;
		if (hbs)
			b ^= l;
		a >>= 1;
	}

	return c & 0xff; // mask
}

void gen_t_box(unsigned *gf, unsigned *fg)
{
	uint32_t i;

	gf[0] = 0;
	for (i = 0; i < 256; i++)
	{
		gf[i] = pd(seki(seki(seki(i, i), seki(i, i)), seki(i, i)), i, 0b100011011) ^ 0b11000110; // pmod(pmod(i,i,normal[E]),i,normal[E]);
	}

	// exit(1);
	printf("static const unsigned short gf[256]={\n");
	for (i = 0; i < 256; i++)
	{
		fg[gf[i]] = i;
		printf("%3d,", gf[i]);
	}
	printf("};\n");
	printf("\n");
	printf("static const unsigned short fg[256]={\n");
	for (i = 0; i < 256; i++)
	{
		printf("%3d,", fg[i]);
	}
	printf("};\n");
	printf("\n");
}

uint8_t box(unsigned char x)
{
	int i;
	//	uint8_t a;
	uint8_t A[8] = {
		0b10001111,
		0b11000111,
		0b11100011,
		0b11110001,
		0b11111000,
		0b01111100,
		0b00111110,
		0b00011111};
	uint8_t c = 0b11000110;
	uint8_t y = 0;

	for (i = 0; i < 8; i++)
	{
		y <<= 1;
		y ^= __builtin_parity((gmult(x, x) & A[i]) ^ c);
	}

	return y;
}

unsigned short vb[6][6] = {0};

unsigned short gf[256] = {0, 1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35, 70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 175, 67, 134, 17, 34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197, 151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115, 230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219, 171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25, 50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9, 18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11, 22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142};
unsigned short fg[256] = {0, 1, 2, 26, 3, 51, 27, 199, 4, 224, 52, 239, 28, 105, 200, 76, 5, 101, 225, 15, 53, 142, 240, 130, 29, 194, 106, 249, 201, 9, 77, 114, 6, 139, 102, 48, 226, 37, 16, 34, 54, 148, 143, 219, 241, 19, 131, 70, 30, 182, 195, 126, 107, 40, 250, 186, 202, 155, 10, 121, 78, 229, 115, 167, 7, 192, 140, 99, 103, 222, 49, 254, 227, 153, 38, 180, 17, 146, 35, 137, 55, 209, 149, 207, 144, 151, 220, 190, 242, 211, 20, 93, 132, 57, 71, 65, 31, 67, 183, 164, 196, 73, 127, 111, 108, 59, 41, 85, 251, 134, 187, 62, 203, 95, 156, 160, 11, 22, 122, 44, 79, 213, 230, 173, 116, 244, 168, 88, 8, 113, 193, 248, 141, 129, 100, 14, 104, 75, 223, 238, 50, 198, 255, 25, 228, 166, 154, 120, 39, 185, 181, 125, 18, 69, 147, 218, 36, 33, 138, 47, 56, 64, 210, 92, 150, 189, 208, 206, 145, 136, 152, 179, 221, 253, 191, 98, 243, 87, 212, 172, 21, 43, 94, 159, 133, 61, 58, 84, 72, 110, 66, 163, 32, 46, 68, 217, 184, 124, 165, 119, 197, 24, 74, 237, 128, 13, 112, 247, 109, 162, 60, 83, 42, 158, 86, 171, 252, 97, 135, 178, 188, 205, 63, 91, 204, 90, 96, 177, 157, 170, 161, 82, 12, 246, 23, 236, 123, 118, 45, 216, 80, 175, 214, 234, 231, 232, 174, 233, 117, 215, 245, 235, 169, 81, 89, 176};

int mltn(int n, int x)
{
	int ret = 1;
	while (n > 0)
	{
		if (n & 1)
			ret = mlt(ret, x); // n の最下位bitが 1 ならば x^(2^i) をかける
		x = mlt(x, x);
		n >>= 1; // n を1bit 左にずらす
	}
	return ret;
}

void mds(int kk)
{
	int i, j;

	printf("van der\n");

	for (i = 2; i < 5; i++)
		vb[0][i] = 1;
	// #pragma omp parallel for private(i, j)
	for (i = 1; i < 6; i++)
	{
		for (j = 2; j < 6; j++)
		{
			vb[i][j] = gf[mltn(i, fg[j])];
			printf("%d,", vb[i][j]);
		}
		printf("\n");
	}
}

void line(uint8_t *a,uint8_t *m,uint8_t *c){
int i,j,k;

for (i = 0; i < 4; i++){
    for (k = 0; k < 4; k++){
        for (j = 0; j < 8; j++){
            c[i*8+j]^=gf[mlt(fg[a[i*8+k]],fg[m[k*8+j]])];
        }
    }
}
}

void matmax(uint8_t g[4][4], uint8_t *h, uint8_t *c)
{
	int i, j, k;
	// GH=0であることの確認。
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 8; j++)
		{
			c[i*8+j] = 0;
			for (k = 0; k < 4; k++)
				c[i*8+j] ^= gf[mlt(fg[g[i][k]], fg[h[k*8+j]])];
			// printf("c%d,", c[i][j]);
		}
		// printf("\n");
	}
	// printf("\n");
}


int oinv(unsigned short b)
{

	if (b == 0)
		return 0;

	return (256 - fg[b]) % (256 - 1) + 1;
}

#define MATRIX_SIZE 4
// 行列の逆行列を計算する関数
void inverseMatrix(uint8_t A[MATRIX_SIZE][MATRIX_SIZE], uint8_t A_inv[MATRIX_SIZE][MATRIX_SIZE])
{
	int i, j, k;
	uint8_t temp;

	// 単位行列を初期化
	for (i = 0; i < MATRIX_SIZE; i++)
	{
		for (j = 0; j < MATRIX_SIZE; j++)
		{
			A_inv[i][j] = (i == j) ? 1 : 0;
		}
	}

	// ガウス・ジョルダン法による逆行列の計算
	for (k = 0; k < MATRIX_SIZE; k++)
	{
		temp = gf[oinv(A[k][k])];
		for (j = 0; j < MATRIX_SIZE; j++)
		{
			A[k][j] = gf[mlt(fg[A[k][j]], fg[temp])];
			A_inv[k][j] = gf[mlt(fg[A_inv[k][j]], fg[(temp)])];
		}
		for (i = 0; i < MATRIX_SIZE; i++)
		{
			if (i != k)
			{
				temp = A[i][k];
				for (j = 0; j < MATRIX_SIZE; j++)
				{
					A[i][j] ^= gf[mlt(fg[A[k][j]], fg[temp])];
					A_inv[i][j] ^= gf[mlt(fg[A_inv[k][j]], fg[temp])];
				}
			}
		}
	}

	for (i = 0; i < MATRIX_SIZE; i++)
	{
		for (j = 0; j < MATRIX_SIZE; j++)
			printf("%d,", A_inv[i][j]);
		printf("\n");
	}
	printf("\n");
}


void milk(uint8_t vc[4][4])
{
	// uint16_t cx[8][8]={0};
	uint8_t i, j, k = 0, l = 0;

	for (i = 3; i < 10; i = i + 2)
	{
		l=0;
		for (j = 2; j < 9; j = j + 2)
		{
			vc[k][l] = oinv((256 + i - j) % 256);
			printf("%d,", vc[k][l]);
			l++;
		}
		l = 0;
		k++;
		printf("\n");
	}
	printf("\n");
	// exit(1);
}

uint8_t table[16][32];
void dec(uint8_t *m, uint8_t *k, uint8_t *inv_ss)
{
	int i, l;
	uint8_t mm[4][8], con[32];

	//milk(der);
	//inverseMatrix(der,snoot);
	memcpy(r, out, 32);
	for (i = 0; i < 16; i++)
		rounder();

	// aes_inv_cipher(out, m, w);
	for (i = 0; i < 16; i++)
	{
		for (int l = 0; l < 32; l++)
		{
			// printf("%d ",k[l]);
			k[l] = table[15 - i][l];
		}
		// printf("\n");
		perm(m, inv_r);
		matmax(snoot, m, con);
		//inv_shift_rows(con);
		for (int l = 0; l < 32; l++)
		{
			if (l % 2 == 0)
				con[l] = inv_ss[con[l]]; // inv_s_box[m[l]];
		}
		sub(con, k);
		memcpy(m,con,32);
		reverse();
	}
	
	printf("Original message (after inv cipher):\n");
	for (i = 0; i < 32; i++)
	{
		printf("%02x ", inv_s_box[s_box[m[i]]]);
	}
	printf("\n");
	
}

int enc(uint8_t *k, uint8_t *m, uint8_t *ss)
{

	int i, j;
	uint8_t mm[4][8], con[32];

	memcpy(r, out, 32);
	//milk(der);
	// exit(1);
	//milk(der);
	for (i = 0; i < 16; i++)
	{
		for (int l = 0; l < 32; l++)
		{
			k[l] = table[i][l];
		}
		// printf("\n");
		rounder();
		add(m, k);
		for (int l = 0; l < 32; l++)
		{
			if (l % 2 == 0)
				m[l] = ss[m[l]]; // s_box[m[l]];
		}
		//shift_rows(m);
		matmax(der, m, con);
		perm(con, r);
		
		memcpy(m,con,32);
	}
	
	printf("Ciphered message:\n");
	for (i = 0; i < 32; i++)
	{
		printf("%02x ", m[i]);
	}
	printf("\n");
	
}

void main()
{
	unsigned int i;
	unsigned char m[32];
	unsigned char k[32] = {0};
	unsigned char s[32], nonce[32] = {103, 198, 105, 115, 81, 255, 74, 236, 41, 205, 186, 171, 242, 251, 227, 70, 124, 194, 84, 248, 27, 232, 231, 141, 118, 90, 46, 99, 51, 159, 201, 154};

	for (i = 0; i < 32; i++)
	{
		m[i] = 0; // 255-i;
		//printf("%d,", nonce[i]);
		k[i] = 0; // random() % 256;
		k[i] ^= nonce[i];
		p[i] = i;
		r[i] = i;
	}
	//uint8_t out[32]; //, mo[256], inv_mo[256];
	uint8_t mm[4][8] = {0}, con[8][8] = {0};

	// mds(4);
	// exit(1);
	random_shuffle(p, 32);
	random_shuffle(r, 32);
	memcpy(s, r, 32);
	for (i = 0; i < 32; i++)
	{
		inv_p[p[i]] = i;
		inv_r[r[i]] = i;
	}
	memcpy(out, r, 32);
	printf("\n");
	printf("Plaintext message:\n");
	for (i = 0; i < 32; i++)
	{
		printf("%02x ", m[i]);
	}
	printf("\n");
	
	uint8_t ss[256], inv_ss[256];
	for (i = 0; i < 256; i++)
		ss[i] = mkbox(i);
	for (i = 0; i < 256; i++)
		inv_ss[ss[i]] = i;
	//memcpy(r,out,32);
	for (i = 0; i < 16; i++)
	{
		for (int j = 0; j < 32; j++)
			k[j] ^= rotl(k[r[j]], 3);
		rounder();
		for (int j = 0; j < 32; j++)
			table[i][j] = k[j];
	}
	i=0;
	


	while(i<1000000){
	enc(k, m, ss);
	dec(m, k, inv_ss);
	i++;
	}
}