//  synthwave.h
//  
//  Copyright (C) 2017 Robin Southern   https://github.com/betajaen/synthwave
//  
//  This software may be modified and distributed under the terms
//  of the MIT license.  See the LICENSE file for details.

#include "synthwave.h"
#include <math.h>
#include <stdio.h>

f32 Vector2_Length(Vector* v)
{
  return sqrtf(Vector2_Length2(v));
}

f32 Vector_Length(Vector* v)
{
  return sqrtf(Vector_Length2(v));
}

f32 Vector4_Length(Vector* v)
{
  return sqrtf(Vector4_Length2(v));
}

void Matrix_Multiply(Matrix* out, Matrix* a, Matrix* b)
{
 f32 a00 = a->m[0],  a01 = a->m[1],  a02 = a->m[2],  a03 = a->m[3],
     a10 = a->m[4],  a11 = a->m[5],  a12 = a->m[6],  a13 = a->m[7],
     a20 = a->m[8],  a21 = a->m[9],  a22 = a->m[10], a23 = a->m[11],
     a30 = a->m[12], a31 = a->m[13], a32 = a->m[14], a33 = a->m[15],
     b00 = b->m[0],  b01 = b->m[1],  b02 = b->m[2],  b03 = b->m[3],
     b10 = b->m[4],  b11 = b->m[5],  b12 = b->m[6],  b13 = b->m[7],
     b20 = b->m[8],  b21 = b->m[9],  b22 = b->m[10], b23 = b->m[11],
     b30 = b->m[12], b31 = b->m[13], b32 = b->m[14], b33 = b->m[15];

  out->m[0]  = b00 * a00 + b01 * a10 + b02 * a20 + b03 * a30;
  out->m[1]  = b00 * a01 + b01 * a11 + b02 * a21 + b03 * a31;
  out->m[2]  = b00 * a02 + b01 * a12 + b02 * a22 + b03 * a32;
  out->m[3]  = b00 * a03 + b01 * a13 + b02 * a23 + b03 * a33;
  out->m[4]  = b10 * a00 + b11 * a10 + b12 * a20 + b13 * a30;
  out->m[5]  = b10 * a01 + b11 * a11 + b12 * a21 + b13 * a31;
  out->m[6]  = b10 * a02 + b11 * a12 + b12 * a22 + b13 * a32;
  out->m[7]  = b10 * a03 + b11 * a13 + b12 * a23 + b13 * a33;
  out->m[8]  = b20 * a00 + b21 * a10 + b22 * a20 + b23 * a30;
  out->m[9]  = b20 * a01 + b21 * a11 + b22 * a21 + b23 * a31;
  out->m[10] = b20 * a02 + b21 * a12 + b22 * a22 + b23 * a32;
  out->m[11] = b20 * a03 + b21 * a13 + b22 * a23 + b23 * a33;
  out->m[12] = b30 * a00 + b31 * a10 + b32 * a20 + b33 * a30;
  out->m[13] = b30 * a01 + b31 * a11 + b32 * a21 + b33 * a31;
  out->m[14] = b30 * a02 + b31 * a12 + b32 * a22 + b33 * a32;
  out->m[15] = b30 * a03 + b31 * a13 + b32 * a23 + b33 * a33;
}

void Matrix_Inverse(Matrix* out, Matrix* m)
{
  f32 a00 = m->m[0],  a01 = m->m[1],  a02 = m->m[2],  a03 = m->m[3],
      a10 = m->m[4],  a11 = m->m[5],  a12 = m->m[6],  a13 = m->m[7],
      a20 = m->m[8],  a21 = m->m[9],  a22 = m->m[10], a23 = m->m[11],
      a30 = m->m[12], a31 = m->m[13], a32 = m->m[14], a33 = m->m[15],
      b00 = a00 * a11 - a01 * a10,
      b01 = a00 * a12 - a02 * a10,
      b02 = a00 * a13 - a03 * a10,
      b03 = a01 * a12 - a02 * a11,
      b04 = a01 * a13 - a03 * a11,
      b05 = a02 * a13 - a03 * a12,
      b06 = a20 * a31 - a21 * a30,
      b07 = a20 * a32 - a22 * a30,
      b08 = a20 * a33 - a23 * a30,
      b09 = a21 * a32 - a22 * a31,
      b10 = a21 * a33 - a23 * a31,
      b11 = a22 * a33 - a23 * a32,
      d = (b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06),
      invDet;

  if (!d)
  {
    printf("Bad determinant");
    Matrix_Identity(m); 
    return;
  }

  invDet = 1.0f / d;

  m->m[0]  = (a11 * b11 - a12 * b10 + a13 * b09) * invDet;
  m->m[1]  = (-a01 * b11 + a02 * b10 - a03 * b09) * invDet;
  m->m[2]  = (a31 * b05 - a32 * b04 + a33 * b03) * invDet;
  m->m[3]  = (-a21 * b05 + a22 * b04 - a23 * b03) * invDet;
  m->m[4]  = (-a10 * b11 + a12 * b08 - a13 * b07) * invDet;
  m->m[5]  = (a00 * b11 - a02 * b08 + a03 * b07) * invDet;
  m->m[6]  = (-a30 * b05 + a32 * b02 - a33 * b01) * invDet;
  m->m[7]  = (a20 * b05 - a22 * b02 + a23 * b01) * invDet;
  m->m[8]  = (a10 * b10 - a11 * b08 + a13 * b06) * invDet;
  m->m[9]  = (-a00 * b10 + a01 * b08 - a03 * b06) * invDet;
  m->m[10] = (a30 * b04 - a31 * b02 + a33 * b00) * invDet;
  m->m[11] = (-a20 * b04 + a21 * b02 - a23 * b00) * invDet;
  m->m[12] = (-a10 * b09 + a11 * b07 - a12 * b06) * invDet;
  m->m[13] = (a00 * b09 - a01 * b07 + a02 * b06) * invDet;
  m->m[14] = (-a30 * b03 + a31 * b01 - a32 * b00) * invDet;
  m->m[15] = (a20 * b03 - a21 * b01 + a22 * b00) * invDet;
}

void Matrix_RotX(Matrix* m, f32 x)
{
  f32 c = cosf(x), s = sinf(x);

  m->M[0][0] = 1.0f;  m->M[0][1] = 0.0f;  m->M[0][2] = 0.0f;  m->M[0][3] = 0.0f;
  m->M[1][0] = 0.0f;  m->M[1][1] = c;     m->M[1][2] = -s;    m->M[1][3] = 0.0f;
  m->M[2][0] = 0.0f;  m->M[2][1] = s;     m->M[2][2] = c;     m->M[2][3] = 0.0f;
  m->M[3][0] = 0.0f;  m->M[3][1] = 0.0f;  m->M[3][2] = 0.0f;  m->M[3][3] = 1.0f;
}

void Matrix_RotY(Matrix* m, f32 y)
{
  f32 c = cosf(y), s = sinf(y);

  m->M[0][0] = c;     m->M[0][1] = 0.0f;  m->M[0][2] = s;     m->M[0][3] = 0.0f;
  m->M[1][0] = 0.0f;  m->M[1][1] = 1.0f;  m->M[1][2] = 0.0f;  m->M[1][3] = 0.0f;
  m->M[2][0] = -s;    m->M[2][1] = 0.0f;  m->M[2][2] = c;     m->M[2][3] = 0.0f;
  m->M[3][0] = 0.0f;  m->M[3][1] = 0.0f;  m->M[3][2] = 0.0f;  m->M[3][3] = 1.0f;
}

void Matrix_RotZ(Matrix* m, f32 z)
{
  f32 c = cosf(z), s = sinf(z);

  m->M[0][0] = c;     m->M[0][1] = -s;    m->M[0][2] = 0.0f;  m->M[0][3] = 0.0f;
  m->M[1][0] = s;     m->M[1][1] = c;     m->M[1][2] = 0.0f;  m->M[1][3] = 0.0f;
  m->M[2][0] = 0.0f;  m->M[2][1] = 0.0f;  m->M[2][2] = 1.0f;  m->M[2][3] = 0.0f;
  m->M[3][0] = 0.0f;  m->M[3][1] = 0.0f;  m->M[3][2] = 0.0f;  m->M[3][3] = 1.0f;
}

int main(int argc, char** argv)
{

}
