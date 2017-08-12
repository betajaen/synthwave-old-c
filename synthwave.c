//  synthwave.h
//  
//  Copyright (C) 2017 Robin Southern   https://github.com/betajaen/synthwave
//  
//  This software may be modified and distributed under the terms
//  of the MIT license.  See the LICENSE file for details.

#include "synthwave.h"
#include <math.h>

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

int main(int argc, char** argv)
{
}
