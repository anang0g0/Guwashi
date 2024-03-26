/*

  AES Versoin 1.4 pjc.co.jp

  Copyright (C) 2009  PJC.CO.JP

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

/*
 
　本ソースは、AESを学ぶために作成したソースです。
　当ホームページより配信しているPDFファイルを参照して組込しておりますので、
　PDFファイルもご参考ください。なおPDFファイルは英語表記となっております。
  http://free.pjc.co.jp/AES/fips-197.pdf

*/

/*

  ソースコメント内容について補足

  P.** は、当ホームページ配信の仕様書PDFファイルと照らし合わせてご覧ください。
  http://free.pjc.co.jp/AES/fips-197.pdf

  なおこの P.** 表記は、PDFのページではなく、仕様書の実ページと対応しております。
  例) P.1 とあった場合、
  Announcing the ADVANCED ENCRYPTION STANDARD (AES) のページではなく
  PDFでの P.5 Specification for the ADVANCED ENCRYPTION STANDARD (AES) (目次のページ)が 
  P.1 となります。

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NB 4                          /* 128bit 固定として規格されている(データの長さ) */
#define NBb 16

/************************************************************/
unsigned char key[32];
int w[60];                            /* FIPS 197 P.19 5.2 Key Expansion */
int data[NB];
int nk;                               /* 4,6,8(128,192,256 bit) 鍵の長さ */
int nr;                               /* 10,12,14 ラウンド数 */

void SubBytes(int data[]);            /* FIPS 197  P.16 Figure  6 */
void ShiftRows(int data[]);           /* FIPS 197  P.17 Figure  8 */
void MixColumns(int *);               /* FIPS 197  P.18 Figure  9 */
void AddRoundKey(int *,int);          /* FIPS 197  P.19 Figure 10 */
int SubWord(int in);                  /* FIPS 197  P.20 Figure 11 */ /* FIPS 197  P.19  5.2 */
int RotWord(int in);                  /* FIPS 197  P.20 Figure 11 */ /* FIPS 197  P.19  5.2 */
void KeyExpansion(void *);            /* FIPS 197  P.20 Figure 11 */
void invShiftRows(int data[]);        /* FIPS 197  P.22 Figure 13 */
void invSubBytes(int data[]);         /* FIPS 197  P.22 5.3.2 */
void invMixColumns(int data[],int n); /* FIPS 197  P.23 5.3.3 */

int Cipher(int *);                    /* FIPS 197  P.15 Figure  5 */
int invCipher(int *);                 /* FIPS 197  P.21 Figure 12 */

/************************************************************/
void datadump(char c[],void *dt,int len)
{
  int i;
  unsigned char *cdt = (unsigned char *)dt;
  printf("%s",c);
  for(i=0;i<len*4;i++)
  {
    printf("%02x",cdt[i]);
  }
  printf("\n");
}

/************************************************************/
int main()
{
  unsigned char keys[]={0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
                        0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
                        0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
                        0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f};
  unsigned char init[]={0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,
                        0x88,0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xff};

  /* FIPS 197  P.35 Appendix C.1 AES-128 Test */
  /*
  memcpy(key,keys,16);
  nk = 4;               //鍵の長さ 4,6,8(128,192,256 bit)
  nr = nk + 6;          //ラウンド数 10,12,14

  KeyExpansion(key);    //暗号化するための鍵の準備
  memcpy(data,init,NBb); //NBにて 4ワード 16バイトと定義している

  printf("  <FIPS 197  P.35 Appendix C.1 AES-128 TEST>\n\n");
  datadump("PLAINTEXT: ",data,4);
  datadump("KEY:       ",key,4);
  Cipher(data);
  datadump("暗号化:    ",data,4);
  invCipher(data);
  datadump("復号化:    ",data,4);
  printf("\n");
*/
  /* FIPS 197  P.38 Appendix C.2 AES-192 Test */
  /*
  memcpy(key,keys,24);
  nk = 6;               //鍵の長さ 4,6,8(128,192,256 bit)
  nr = nk + 6;          //ラウンド数 10,12,14

  KeyExpansion(key);    //暗号化するための鍵の準備
  memcpy(data,init,NBb); //NBにて 4ワード 16バイトと定義している

  printf("  <FIPS 197  P.38 Appendix C.2 AES-192 TEST>\n\n");
  datadump("PLAINTEXT: ",data,4);
  datadump("KEY:       ",key,6);
  Cipher(data);
  datadump("暗号化:    ",data,4);
  invCipher(data);
  datadump("復号化:    ",data,4);
  printf("\n");
*/
  /* FIPS 197  P.42 Appendix C.3 AES-256 Test */
  memcpy(key,keys,32);
  nk = 8;               //鍵の長さ 4,6,8(128,192,256 bit)
  nr = nk + 6;          //ラウンド数 10,12,14

  KeyExpansion(key);    //暗号化するための鍵の準備
  memcpy(data,init,NBb); //NBにて 4ワード 16バイトと定義している

  printf("  <FIPS 197  P.42 Appendix C.3 AES-256 TEST>\n\n");
  datadump("PLAINTEXT: ",data,4);
  datadump("KEY:       ",key,8);
  int i=0;
  while(i<1000000){
  Cipher(data);
  //datadump("暗号化:    ",data,4);
  invCipher(data);
  i++;
  }
  datadump("復号化:    ",data,4);
  return 0;
}

/************************************************************/
/* FIPS 197  P.16 Figure 7 */
int Sbox[256]={
  0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
  0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
  0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
  0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
  0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
  0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
  0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
  0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
  0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
  0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
  0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
  0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
  0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
  0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
  0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
  0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16
};

/************************************************************/
/* FIPS 197  P.22 Figure 14 */
int invSbox[256]={
  0x52,0x09,0x6a,0xd5,0x30,0x36,0xa5,0x38,0xbf,0x40,0xa3,0x9e,0x81,0xf3,0xd7,0xfb,
  0x7c,0xe3,0x39,0x82,0x9b,0x2f,0xff,0x87,0x34,0x8e,0x43,0x44,0xc4,0xde,0xe9,0xcb,
  0x54,0x7b,0x94,0x32,0xa6,0xc2,0x23,0x3d,0xee,0x4c,0x95,0x0b,0x42,0xfa,0xc3,0x4e,
  0x08,0x2e,0xa1,0x66,0x28,0xd9,0x24,0xb2,0x76,0x5b,0xa2,0x49,0x6d,0x8b,0xd1,0x25,
  0x72,0xf8,0xf6,0x64,0x86,0x68,0x98,0x16,0xd4,0xa4,0x5c,0xcc,0x5d,0x65,0xb6,0x92,
  0x6c,0x70,0x48,0x50,0xfd,0xed,0xb9,0xda,0x5e,0x15,0x46,0x57,0xa7,0x8d,0x9d,0x84,
  0x90,0xd8,0xab,0x00,0x8c,0xbc,0xd3,0x0a,0xf7,0xe4,0x58,0x05,0xb8,0xb3,0x45,0x06,
  0xd0,0x2c,0x1e,0x8f,0xca,0x3f,0x0f,0x02,0xc1,0xaf,0xbd,0x03,0x01,0x13,0x8a,0x6b,
  0x3a,0x91,0x11,0x41,0x4f,0x67,0xdc,0xea,0x97,0xf2,0xcf,0xce,0xf0,0xb4,0xe6,0x73,
  0x96,0xac,0x74,0x22,0xe7,0xad,0x35,0x85,0xe2,0xf9,0x37,0xe8,0x1c,0x75,0xdf,0x6e,
  0x47,0xf1,0x1a,0x71,0x1d,0x29,0xc5,0x89,0x6f,0xb7,0x62,0x0e,0xaa,0x18,0xbe,0x1b,
  0xfc,0x56,0x3e,0x4b,0xc6,0xd2,0x79,0x20,0x9a,0xdb,0xc0,0xfe,0x78,0xcd,0x5a,0xf4,
  0x1f,0xdd,0xa8,0x33,0x88,0x07,0xc7,0x31,0xb1,0x12,0x10,0x59,0x27,0x80,0xec,0x5f,
  0x60,0x51,0x7f,0xa9,0x19,0xb5,0x4a,0x0d,0x2d,0xe5,0x7a,0x9f,0x93,0xc9,0x9c,0xef,
  0xa0,0xe0,0x3b,0x4d,0xae,0x2a,0xf5,0xb0,0xc8,0xeb,0xbb,0x3c,0x83,0x53,0x99,0x61,
  0x17,0x2b,0x04,0x7e,0xba,0x77,0xd6,0x26,0xe1,0x69,0x14,0x63,0x55,0x21,0x0c,0x7d
};

/************************************************************/
/* FIPS 197  P.15 Figure 5 */ //暗号化
int Cipher(int data[])
{
  int i;

  AddRoundKey(data,0);

  for(i=1;i<nr;i++)
  {
    SubBytes(data);
    ShiftRows(data);
    MixColumns(data);
    AddRoundKey(data,i);
  }

  SubBytes(data);
  ShiftRows(data);
  AddRoundKey(data,i);
  return(i);
}

/************************************************************/
/* FIPS 197  P.21 Figure 12 */ //復号化
int invCipher(int data[])
{
  int i;

  AddRoundKey(data,nr);

  for(i=nr-1;i>0;i--)
  {
    invShiftRows(data);
    invSubBytes(data);
    AddRoundKey(data,i);
    invMixColumns(data,i);
  }

  invShiftRows(data);
  invSubBytes(data);
  AddRoundKey(data,0);
  return(nr);
}

/************************************************************/
/* FIPS 197  P.16 Figure 6 */
void SubBytes(int data[])
{
  int i,j;
  unsigned char *cb=(unsigned char*)data;
  for(i=0;i<NBb;i+=4)//理論的な意味から二重ループにしているが意味は無い
  {
    for(j=0;j<4;j++)
    {
      cb[i+j] = Sbox[cb[i+j]];
    }
  }
}

/************************************************************/
/* FIPS 197  P.22 5.3.2 */
void invSubBytes(int data[])
{
  int i,j;
  unsigned char *cb=(unsigned char*)data;
  for(i=0;i<NBb;i+=4)//理論的な意味から二重ループにしているが意味は無い
  {
    for(j=0;j<4;j++)
    {
      cb[i+j] = invSbox[cb[i+j]];
    }
  }
}

/************************************************************/
/* FIPS 197  P.17 Figure 8 */
void ShiftRows(int data[])
{
  int i,j,i4;
  unsigned char *cb=(unsigned char*)data;
  unsigned char cw[NBb];
  memcpy(cw,cb,sizeof(cw));
  for(i=0;i<NB;i+=4)
  {
    i4 = i*4;
    for(j=1;j<4;j++)
    {
      cw[i4+j+0*4] = cb[i4+j+((j+0)&3)*4];
      cw[i4+j+1*4] = cb[i4+j+((j+1)&3)*4];
      cw[i4+j+2*4] = cb[i4+j+((j+2)&3)*4];
      cw[i4+j+3*4] = cb[i4+j+((j+3)&3)*4];
    }
  }
  memcpy(cb,cw,sizeof(cw));
}

/************************************************************/
/* FIPS 197  P.22 Figure 13 */
void invShiftRows(int data[])
{
  int i,j,i4;
  unsigned char *cb=(unsigned char*)data;
  unsigned char cw[NBb];
  memcpy(cw,cb,sizeof(cw));
  for(i=0;i<NB;i+=4)
  {
    i4 = i*4;
    for(j=1;j<4;j++)
    {
      cw[i4+j+((j+0)&3)*4] = cb[i4+j+0*4];
      cw[i4+j+((j+1)&3)*4] = cb[i4+j+1*4];
      cw[i4+j+((j+2)&3)*4] = cb[i4+j+2*4];
      cw[i4+j+((j+3)&3)*4] = cb[i4+j+3*4];
    }
  }
  memcpy(cb,cw,sizeof(cw));
}

/************************************************************/
/* FIPS 197 P.10 4.2 乗算 (n倍) */
int mul(int dt,int n)
{
  int i,x=0;
  for(i=8;i>0;i>>=1)
  {
    x <<= 1;
    if(x&0x100)
      x = (x ^ 0x1b) & 0xff;
    if((n & i))
      x ^= dt;
  }
  return(x);
}

/************************************************************/
int dataget(void* data,int n)
{
  return(((unsigned char*)data)[n]);
}

/************************************************************/
/* FIPS 197  P.18 Figure 9 */
void MixColumns(int data[])
{
  int i,i4,x;
  for(i=0;i<NB;i++)
  {
    i4 = i*4;
    x  =  mul(dataget(data,i4+0),2) ^
          mul(dataget(data,i4+1),3) ^
          mul(dataget(data,i4+2),1) ^
          mul(dataget(data,i4+3),1);
    x |= (mul(dataget(data,i4+1),2) ^
          mul(dataget(data,i4+2),3) ^
          mul(dataget(data,i4+3),1) ^
          mul(dataget(data,i4+0),1)) << 8;
    x |= (mul(dataget(data,i4+2),2) ^
          mul(dataget(data,i4+3),3) ^
          mul(dataget(data,i4+0),1) ^
          mul(dataget(data,i4+1),1)) << 16;
    x |= (mul(dataget(data,i4+3),2) ^
          mul(dataget(data,i4+0),3) ^
          mul(dataget(data,i4+1),1) ^
          mul(dataget(data,i4+2),1)) << 24;
    data[i] = x;
  } 
}

/************************************************************/
/* FIPS 197  P.23 5.3.3 */
void invMixColumns(int data[],int n)
{
  int i,i4,x;
  for(i=0;i<NB;i++)
  {
    i4 = i*4;
    x  =  mul(dataget(data,i4+0),14) ^
          mul(dataget(data,i4+1),11) ^
          mul(dataget(data,i4+2),13) ^
          mul(dataget(data,i4+3), 9);
    x |= (mul(dataget(data,i4+1),14) ^
          mul(dataget(data,i4+2),11) ^
          mul(dataget(data,i4+3),13) ^
          mul(dataget(data,i4+0), 9)) << 8;
    x |= (mul(dataget(data,i4+2),14) ^
          mul(dataget(data,i4+3),11) ^
          mul(dataget(data,i4+0),13) ^
          mul(dataget(data,i4+1), 9)) << 16;
    x |= (mul(dataget(data,i4+3),14) ^
          mul(dataget(data,i4+0),11) ^
          mul(dataget(data,i4+1),13) ^
          mul(dataget(data,i4+2), 9)) << 24;
    data[i] = x;
  } 
}

/************************************************************/
/* FIPS 197  P.19 Figure 10 */
void AddRoundKey(int data[],int n)
{
  int i;
  for(i=0;i<NB;i++)
  {
    data[i] ^= w[i+NB*n];
  }
}

/************************************************************/
/* FIPS 197  P.20 Figure 11 */ /* FIPS 197  P.19  5.2 */
int SubWord(int in)
{
  int inw=in;
  unsigned char *cin=(unsigned char*)&inw;
  cin[0] = Sbox[cin[0]];
  cin[1] = Sbox[cin[1]];
  cin[2] = Sbox[cin[2]];
  cin[3] = Sbox[cin[3]];
  return(inw);
}

/************************************************************/
/* FIPS 197  P.20 Figure 11 */ /* FIPS 197  P.19  5.2 */
int RotWord(int in)
{
  int inw=in,inw2=0;
  unsigned char *cin=(unsigned char*)&inw;
  unsigned char *cin2=(unsigned char*)&inw2;
  cin2[0] = cin[1];
  cin2[1] = cin[2];
  cin2[2] = cin[3];
  cin2[3] = cin[0];
  return(inw2);
}

/************************************************************/
/* FIPS 197  P.20 Figure 11 */
void KeyExpansion(void *key)
{
  /* FIPS 197  P.27 Appendix A.1 Rcon[i/Nk] */ //又は mulを使用する
  int Rcon[10]={0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1b,0x36};
  int i,temp;

  memcpy(w,key,nk*4);
  for(i=nk;i<NB*(nr+1);i++)
  {
    temp = w[i-1];
    if((i%nk) == 0)
      temp = SubWord(RotWord(temp)) ^ Rcon[(i/nk)-1];
    else if(nk > 6 && (i%nk) == 4)
      temp = SubWord(temp);
    w[i] = w[i-nk] ^ temp;
  }
}

