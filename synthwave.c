//  synthwave.h
//  
//  Copyright (C) 2017 Robin Southern   https://github.com/betajaen/synthwave
//  
//  This software may be modified and distributed under the terms
//  of the MIT license.  See the LICENSE file for details.

#include "synthwave.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <limits.h>
#include <assert.h>

#include <SDL.h>
#include <SDL_main.h>

#define $Assert(X, T) assert(X && T)
#define $MaybeFunction(FN, ...)  do { if (FN) { FN(__VA_ARGS__); } } while(0)

typedef struct
{
  Synthwave_Description desc;

  struct
  {
    SDL_Window*   window;
    SDL_Renderer* renderer;
    u32           flags;
  } Screen;

  struct
  {
    u8     num, padding[3];
    Colour colour[256];
  } Palette;
  
  struct
  {
    f32 frameSec, fixedSec;
    Timer frameTimer, fixedLimiter, fpsTimer;
    u32 accumulator;
  } Time;

} Synthwave_Internal;

Synthwave_Internal $$;
Synthwave_Interface $;

f32 $NormaliseDegrees(f32 deg)
{
  deg = fmodf(deg + 180.0f, 360.0f);
  if (deg < 0.0f)
      deg += 360;
  return deg - 180;
}

f32 $NormaliseRadians(f32 rad)
{
  rad = fmodf(rad + $PI, $TWO_PI);
  if (rad < 0.0f)
      rad += $TWO_PI;
  return rad - $PI;
}

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

void Vector_MultiplyMatrix(Vector* v, Matrix* m)
{
  v->x = v->x * m->row[0].x + v->y * m->row[1].x + v->z * m->row[2].x;
  v->y = v->x * m->row[0].y + v->y * m->row[1].y + v->z * m->row[2].y;
  v->z = v->x * m->row[0].z + v->y * m->row[1].z + v->z * m->row[2].z;
  v->w = 1.0f;
}

void Vector4_MultiplyMatrix(Vector* v, Matrix* m)
{
  v->x = v->x * m->row[0].x + v->y * m->row[1].x + v->z * m->row[2].x + v->w * m->row[3].x;
  v->y = v->x * m->row[0].y + v->y * m->row[1].y + v->z * m->row[2].y + v->w * m->row[3].y;
  v->z = v->x * m->row[0].z + v->y * m->row[1].z + v->z * m->row[2].z + v->w * m->row[3].z;
  v->w = v->x * m->row[0].w + v->y * m->row[1].w + v->z * m->row[2].w + v->w * m->row[3].w;
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

void Matrix_Transpose(Matrix* out, Matrix* m)
{
  out->M[0][0] = m->M[0][0];
  out->M[0][1] = m->M[1][0];
  out->M[0][2] = m->M[2][0];
  out->M[0][3] = m->M[3][0];
  out->M[1][0] = m->M[0][1];
  out->M[1][1] = m->M[1][1];
  out->M[1][2] = m->M[2][1];
  out->M[1][3] = m->M[3][1];
  out->M[2][0] = m->M[0][2];
  out->M[2][1] = m->M[1][2];
  out->M[2][2] = m->M[2][2];
  out->M[2][3] = m->M[3][2];
  out->M[3][0] = m->M[0][3];
  out->M[3][1] = m->M[1][3];
  out->M[3][2] = m->M[2][3];
  out->M[3][3] = m->M[3][3];
}

void Matrix_Rotate_Roll(Matrix* m, f32 roll)
{
  f32 cr = cosf(roll), sr = sinf(roll);

  m->M[0][0] = 1.0f;  m->M[0][1] = 0.0f;  m->M[0][2] = 0.0f;  m->M[0][3] = 0.0f;
  m->M[1][0] = 0.0f;  m->M[1][1] = cr;    m->M[1][2] = -sr;   m->M[1][3] = 0.0f;
  m->M[2][0] = 0.0f;  m->M[2][1] = sr;    m->M[2][2] = cr;    m->M[2][3] = 0.0f;
  m->M[3][0] = 0.0f;  m->M[3][1] = 0.0f;  m->M[3][2] = 0.0f;  m->M[3][3] = 1.0f;
}

void Matrix_Rotate_Pitch(Matrix* m, f32 pitch)
{
  f32 cp = cosf(pitch), sp = sinf(pitch);

  m->M[0][0] = cp;    m->M[0][1] = 0.0f;   m->M[0][2] = sp;    m->M[0][3] = 0.0f;
  m->M[1][0] = 0.0f;  m->M[1][1] = 1.0f;   m->M[1][2] = 0.0f;  m->M[1][3] = 0.0f;
  m->M[2][0] = -sp;   m->M[2][1] = 0.0f;   m->M[2][2] = cp;    m->M[2][3] = 0.0f;
  m->M[3][0] = 0.0f;  m->M[3][1] = 0.0f;   m->M[3][2] = 0.0f;  m->M[3][3] = 1.0f;
}

void Matrix_Rotate_Yaw(Matrix* m, f32 yaw)
{
  f32 cy = cosf(yaw), sy = sinf(yaw);

  m->M[0][0] = cy;    m->M[0][1] = -sy;    m->M[0][2] = 0.0f;  m->M[0][3] = 0.0f;
  m->M[1][0] = sy;    m->M[1][1] = cy;     m->M[1][2] = 0.0f;  m->M[1][3] = 0.0f;
  m->M[2][0] = 0.0f;  m->M[2][1] = 0.0f;   m->M[2][2] = 1.0f;  m->M[2][3] = 0.0f;
  m->M[3][0] = 0.0f;  m->M[3][1] = 0.0f;   m->M[3][2] = 0.0f;  m->M[3][3] = 1.0f;
}

void Matrix_LookAt(Matrix* m, Vector* position, Vector* target, Vector* up)
{
  Vector z;
  Vector_Sub(&z, target, position);
  Vector_Normalise(&z);

  Vector x;
  Vector_Cross(&x, up, &z);
  Vector_Normalise(&x);

  Vector y;
  Vector_Cross(&y, &z, &x);
  
  m->M[0][0] = x.x;   m->M[0][1] = y.x;   m->M[0][2] = z.x;   m->M[0][3] = 0.0f;
  m->M[1][0] = x.y;   m->M[1][1] = y.y;   m->M[1][2] = z.y;   m->M[1][3] = 0.0f;
  m->M[2][0] = x.z;   m->M[2][1] = y.z;   m->M[2][2] = z.z;   m->M[2][3] = 0.0f;

  m->M[3][0] = -Vector_Dot(&x, position);
  m->M[3][1] = -Vector_Dot(&z, position);
  m->M[3][2] = -Vector_Dot(&x, position);
  m->M[3][3] = 1.0f;
}

void Rotation_Normalize(Rotation* rotation)
{
  rotation->roll  = $NormaliseDegrees(rotation->roll);
  rotation->pitch = $NormaliseDegrees(rotation->pitch);
  rotation->yaw   = $NormaliseDegrees(rotation->yaw);
}

u32 Synthwave_Time_GetTicks()
{
  return SDL_GetTicks();
}

void Synthwave_Palette_Add(Colour* colour)
{
  $Assert($$.Palette.num <= 255, "Maximum number of colours exceeded");
  $$.Palette.colour[$$.Palette.num++] = *colour;
}

void Synthwave_Palette_AddRgb(u8 r, u8 g, u8 b)
{
  Colour col = {.r = r, .g = g, .b = b};
  Synthwave_Palette_Add(&col);
}

void Synthwave_Palette_AddU32(u32 colour)
{
  Colour col;
  col.b = colour & 0xFF;
  colour >>= 8;
  col.g = colour & 0xFF;
  colour >>= 8;
  col.r = colour & 0xFF;
  colour >>= 8;
  
  Synthwave_Palette_Add(&col);
}

u8 Synthwave_Palette_GetCount()
{
  return $$.Palette.num;
}

void Synthwave_Palette_Reset()
{
  $$.Palette.num = 0;
  memset($$.Palette.colour, 0, sizeof($$.Palette.colour));
}

void Synthwave_Frame()
{
  $.Time.frameMs = Timer_GetTimeAndReset(&$$.Time.frameTimer);
  $.Time.frame   = $.Time.frameMs / 1000.0f;

  SDL_Event event;

  while (SDL_PollEvent(&event))
  {
    switch(event.type)
    {
      case SDL_QUIT:
      {
        $.quit = true;
      }
    }
  }

  u32 frameTime = $.Time.frameMs;
  if (frameTime > 250)
    frameTime = 250;
  
  if ($$.desc.Events.OnFixedFrame != NULL)
  {
    $$.Time.accumulator += frameTime;

    while($$.Time.accumulator >= $$.desc.Time.fixedMs)
    {
      $.Time.fixedMs = $$.desc.Time.fixedMs;
      $.Time.fixed   = $$.desc.Time.fixedMs / 1000.0f;
      
      $$.desc.Events.OnFixedFrame();
      
      $$.Time.accumulator += $.Time.fixedMs;
    }
  }

  $MaybeFunction($$.desc.Events.OnFrame, /**/);
  
  SDL_RenderPresent($$.Screen.renderer);
  $.Time.numFrames++;
}

#if $IsWindows == 1
static void Synthwave_Win32_SetupConsoleHandler();
#endif

int main(int argc, char** argv)
{
  memset(&$$, 0, sizeof($$));
  memset(&$,  0, sizeof($));
  
  $$.desc.Screen.title = NULL;
  $$.desc.Screen.x     = INT32_MAX;
  $$.desc.Screen.y     = INT32_MAX;
  $$.desc.Screen.w     = 320;
  $$.desc.Screen.h     = 200;
  $$.desc.Screen.scale = 2;
  $$.desc.Time.frameMs = 20;
  $$.desc.Time.fixedMs = 20;
  
  $Setup(&$$.desc);
  
  $.Palette.Add = Synthwave_Palette_Add;
  $.Palette.AddRgb = Synthwave_Palette_AddRgb;
  $.Palette.AddU32 = Synthwave_Palette_AddU32;
  $.Palette.GetCount = Synthwave_Palette_GetCount;
  $.Palette.Reset    = Synthwave_Palette_Reset;
  $.Time.GetTicks = Synthwave_Time_GetTicks;
  $$.desc.type = $Clamp($$.desc.type, ST_Windowed, ST_Windowed);

  $$.desc.Time.fixedMs = $Max(4, $$.desc.Time.fixedMs);
  $$.desc.Time.frameMs = $Max(4, $$.desc.Time.frameMs);
  $$.Time.fixedSec = $$.desc.Time.fixedMs / 1000.0f;
  $$.Time.frameSec = $$.desc.Time.frameMs / 1000.0f;


  $$.Screen.flags = 0;

  switch($$.desc.type)
  {
    case ST_Windowed:
      $$.Screen.flags = SDL_INIT_VIDEO | SDL_INIT_AUDIO | SDL_INIT_EVENTS | SDL_INIT_TIMER | SDL_INIT_NOPARACHUTE;
    break;
    case ST_Console:
      $$.Screen.flags = SDL_INIT_TIMER | SDL_INIT_EVENTS | SDL_INIT_NOPARACHUTE;
    break;
    case ST_RunOnce:
      $$.Screen.flags = 0;
    break;
  }
  
  if ($$.Screen.flags != 0)
  {
    SDL_Init($$.Screen.flags);
  }
  
   if ($$.Screen.flags != 0)
  {

    if (($$.Screen.flags & SDL_INIT_VIDEO) == 0)
    {
#if $IsWindows == 1
      Synthwave_Win32_SetupConsoleHandler();
#endif
    }
    
    if (($$.Screen.flags & SDL_INIT_VIDEO) != 0)
    {
      i32 x = $$.desc.Screen.x == INT_MAX ? SDL_WINDOWPOS_CENTERED : $$.desc.Screen.x,
          y = $$.desc.Screen.y == INT_MAX ? SDL_WINDOWPOS_CENTERED : $$.desc.Screen.y;

      u32 w = $$.desc.Screen.w * $$.desc.Screen.scale, 
          h = $$.desc.Screen.h * $$.desc.Screen.scale;

      const char* title = $$.desc.Screen.title == NULL ? "Synthwave" : $$.desc.Screen.title;

      $$.Screen.window = SDL_CreateWindow(title, x, y, w, h, SDL_WINDOW_SHOWN);
      
      SDL_GetWindowPosition($$.Screen.window, &$.Screen.x, &$.Screen.y);
      SDL_GetWindowSize($$.Screen.window, (i32*) &$.Screen.width, (i32*) &$.Screen.height);
      
      $$.Screen.renderer = SDL_CreateRenderer(
        $$.Screen.window, 
        -1, 
        SDL_RENDERER_ACCELERATED
      );

    }
    
    Timer_Reset(&$$.Time.fixedLimiter);
    Timer_Reset(&$$.Time.frameTimer);
    Timer_Reset(&$$.Time.fpsTimer);

#if $IsWindows == 1
    while($.quit == false)
    {
      Synthwave_Frame();
      
      u32 frameMs = Timer_GetTime(&$$.Time.frameTimer);

      if (frameMs < $$.desc.Time.frameMs)
      {
        SDL_Delay($$.desc.Time.frameMs - frameMs);
      }
    }
#else
    emscripten_set_main_loop(Synthwave_Frame, 0, 1);
#endif
  }

  $MaybeFunction($$.desc.Events.OnQuit, /**/);
  
  if (($$.Screen.flags & SDL_INIT_VIDEO) != 0 && $$.Screen.window != NULL && $$.Screen.renderer != NULL)
  {
    SDL_DestroyRenderer($$.Screen.renderer);
    $$.Screen.renderer;

    SDL_DestroyWindow($$.Screen.window);
    $$.Screen.window = NULL;
  }

  if ($$.Screen.flags != 0)
  {
    SDL_Quit();
  }

  memset(&$$, 0, sizeof($$));
  memset(&$, 0, sizeof($));

  return 0;
}


#include <windows.h>

#if $IsWindows == 1
static BOOL WINAPI Synthwave_Win32_ConsoleHandler(DWORD signal) {

    if (signal == CTRL_C_EVENT)
    {
      $.quit = true;
    }
    return TRUE;
}

static void Synthwave_Win32_SetupConsoleHandler()
{
  if (!SetConsoleCtrlHandler(Synthwave_Win32_ConsoleHandler, TRUE))
  {
    printf("\nERROR: Could not set control handler"); 
  }
}
#endif
