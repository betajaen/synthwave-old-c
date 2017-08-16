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

#define $Assert(X, T)            assert(X && T)
#define $MaybeFunction(FN, ...)  do { if (FN) { FN(__VA_ARGS__); } } while(0)

typedef enum
{
  ST_DrawMesh = 0,
} Synthwave_Scene_DrawCommand_Type;

typedef struct
{
  Synthwave_Scene_DrawCommand_Type type;
  union
  {
    struct
    {
      Matrix                       worldMatrix;
      Mesh*                        mesh;
      Synthwave_TriangleShader_Fn  triangle;
      Synthwave_FragmentShader_Fn  fragment;
    } drawMesh;
  };
} Synthwave_Scene_DrawCommand;

typedef struct
{
  i16 x, y;
  u16 w, h;

  u8*  palette;
  f32* depth;
  f32* brightness;

  u8*  paletteMem;
  f32* depthMem;
  f32* brightnessMem;

  SDL_Texture* texture;
} Synthwave_FrameBuffer;

typedef struct
{
  Synthwave_FrameBuffer        fb;
  Synthwave_Scene_DrawCommand* drawCmds;
  Matrix                       invView, projection, screen, viewProjectionScreen;
  Vector                       cameraPosition, cameraTarget, lightPosition;
  bool                         vpsOutOfDate;
} Synthwave_Scene;

typedef struct
{
  SDL_Texture* texture;
} Synthwave_Surface;

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
    u8* tempStart, *tempCurrent, *tempEnd;
  } Mem;

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

void Matrix_Projection(Matrix* out, f32 w, f32 h, f32 fovY, f32 nearClip, f32 farClip)
{
  f32 aspect    = (f32) w / (f32) h;
  f32 farPlane  = farClip;
  f32 nearPlane = nearClip;

  f32 width  = tanf(fovY * 0.5f);
  f32 height = width / aspect;
  f32 diff   = farPlane - nearPlane;
  f32 div    = farPlane / diff;

  out->M[0][0] = 1.0f / width;
  out->M[0][1] = 0.0f;
  out->M[0][2] = 0.0f;
  out->M[0][3] = 0.0f;
  out->M[1][0] = 0.0f;
  out->M[1][1] = 1.0f / height;
  out->M[1][2] = 0.0f;
  out->M[1][3] = 0.0f;
  out->M[2][0] = 0.0f;
  out->M[2][1] = 0.0f;
  out->M[2][2] = div;
  out->M[2][3] = 1.0f;
  out->M[3][0] = 0.0f;
  out->M[3][1] = 0.0f;
  out->M[3][2] = -nearPlane * div;
  out->M[3][3] = 0.0f;
}

void Matrix_Screen(Matrix* out, f32 w, f32 h)
{
  f32 halfW = w * 0.5f;
  f32 halfH = h * 0.5f;

  out->M[0][0] = halfW;
  out->M[0][1] = 0.0f;
  out->M[0][2] = 0.0f;
  out->M[0][3] = 0.0f;
  out->M[1][0] = 0.0f;
  out->M[1][1] = -halfH;
  out->M[1][2] = 0.0f;
  out->M[1][3] = 0.0f;
  out->M[2][0] = 0.0f;
  out->M[2][1] = 0.0f;
  out->M[2][2] = 1.0f;
  out->M[2][3] = 0.0f;
  out->M[3][0] = halfW;
  out->M[3][1] = halfH;
  out->M[3][2] = 0.0f;
  out->M[3][3] = 1.0f;
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

void Synthwave_Mem_Temp_Setup(u32 memorySize)
{
  $$.Mem.tempStart   = (u8*) malloc(memorySize);
  $$.Mem.tempCurrent = $$.Mem.tempStart;
  $$.Mem.tempEnd     = $$.Mem.tempStart + memorySize;
}

void Synthwave_Mem_Temp_Shutdown(u32 memorySize)
{
  free($$.Mem.tempStart);
  $$.Mem.tempStart   = NULL;
  $$.Mem.tempCurrent = NULL;
  $$.Mem.tempEnd     = NULL;
}

void Synthwave_Mem_Temp_Reset(u32 memorySize)
{
  $$.Mem.tempCurrent = $$.Mem.tempStart;
}

void* Synthwave_Mem_Temp(u32 numBytes)
{
  void* mem = $$.Mem.tempCurrent;
  $$.Mem.tempCurrent += numBytes;
 
  $Assert($$.Mem.tempCurrent < $$.Mem.tempEnd, "Out of memory in the temporay allocator");
 
  return mem;
}

void* Synthwave_Mem_Heap(void* mem, u32 numBytes)
{
  void* m = NULL;

  if (mem == NULL && numBytes)
  {
    m = malloc(numBytes);
    memset(m, 0, numBytes);
  }
  else if (mem != NULL && numBytes)
    m =  realloc(mem, numBytes);
  else if (mem != NULL && numBytes == 0)
    free(mem);

  return m;
}

void Synthwave_Debug_Log(const char* text)
{
  printf("%s\n", text);
}

void Synthwave_Debug_LogF(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  printf("\n");
  va_end(args);
}
char* Synthwave_Text_Format(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  u32 len = vsnprintf(NULL, 0, fmt, args);

  va_end(args);

  char* text = $.Mem.Temp(len + 1);
  va_start(args, fmt);

#if $IsWindows == 1
  #pragma warning( push )  
  #pragma warning( disable : 4996 ) // This function or variable may be unsafe.
#endif

  vsprintf(text, fmt, args);

#if $IsWindows == 1
  #pragma warning( pop )
#endif

  va_end(args);

  return text;
}

char* Synthwave_Text_FormatHeap(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  u32 len = vsnprintf(NULL, 0, fmt, args);
  va_end(args);

  char* text = $.Mem.Heap(NULL, len + 1);
  va_start(args, fmt);

#if $IsWindows == 1
  #pragma warning( push )  
  #pragma warning( disable : 4996 ) // This function or variable may be unsafe.
#endif

  vsprintf(text, fmt, args);

#if $IsWindows == 1
  #pragma warning( pop )
#endif

  va_end(args);

  return text;
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

Surface Synthwave_Surface_New()
{
  Synthwave_Surface* surface = (Synthwave_Surface*) Synthwave_Mem_Heap(NULL, sizeof(Synthwave_Surface));
  
  surface->texture = SDL_CreateTexture(
    $$.Screen.renderer,
    SDL_PIXELFORMAT_ABGR8888,
    SDL_TEXTUREACCESS_TARGET,
    $.Screen.width, $.Screen.height
  );
  
  SDL_SetTextureBlendMode(surface->texture, SDL_BLENDMODE_BLEND);

  Surface s;
  s.opaque = (u64) surface;
  return s;
}

void Synthwave_Surface_Delete(Surface* surface)
{
  if (surface->opaque == 0)
    return;

  Synthwave_Surface* s = (Synthwave_Surface*) (surface->opaque);
  
  SDL_DestroyTexture(s->texture);
  s->texture = NULL;
  
  Synthwave_Mem_Heap(s, 0);

  surface->opaque = 0;
}

void Synthwave_Surface_Render(Surface* surface)
{
  if (surface->opaque == 0)
    return;

  Synthwave_Surface* s = (Synthwave_Surface*) (surface->opaque);
  
  SDL_RenderCopy($$.Screen.renderer, s->texture, NULL, NULL);
}

Scene Synthwave_Scene_NewXywh(i16 x, i16 y, u16 w, u16 h)
{
  Synthwave_Scene* scene = (Synthwave_Scene*) Synthwave_Mem_Heap(NULL, sizeof(Synthwave_Scene));
  Array_New(scene->drawCmds, 64);
  
  scene->fb.x = x;
  scene->fb.y = y;
  scene->fb.w = w;
  scene->fb.h = h;
  scene->fb.paletteMem     = $.Mem.Heap(NULL, sizeof(u8)  * (w * 2) + (w * h));
  scene->fb.depthMem       = $.Mem.Heap(NULL, sizeof(f32) * (w * 2) + (w * h));
  scene->fb.brightnessMem  = $.Mem.Heap(NULL, sizeof(f32) * (w * 2) + (w * h));
  
  scene->fb.palette    = scene->fb.paletteMem    + scene->fb.w;
  scene->fb.depth      = scene->fb.depthMem      + scene->fb.w;
  scene->fb.brightness = scene->fb.brightnessMem + scene->fb.w;
  scene->fb.texture = SDL_CreateTexture(
    $$.Screen.renderer,
    SDL_PIXELFORMAT_BGR888,
    SDL_TEXTUREACCESS_STREAMING,
    w, h
  );
  
  Matrix_Projection(&scene->projection, scene->fb.w, scene->fb.h, 70.0f, 1.0f, 100.0f);
  Matrix_Screen(&scene->screen, scene->fb.w, scene->fb.h);

  Scene s;
  s.opaque = (u64) scene;
  return s;
}

Scene Synthwave_Scene_New()
{
  return Synthwave_Scene_NewXywh(0,0,  $.Screen.width, $.Screen.height);
}

void Synthwave_Scene_Delete(Scene* scene)
{
  if (scene->opaque == 0)
    return;

  Synthwave_Scene* s = (Synthwave_Scene*) (scene->opaque);

  Array_Delete(s->drawCmds);
  
  Synthwave_Mem_Heap(s->fb.brightnessMem, 0);
  Synthwave_Mem_Heap(s->fb.paletteMem, 0);
  Synthwave_Mem_Heap(s->fb.depthMem, 0);
  Synthwave_Mem_Heap(s, 0);

  scene->opaque = 0;

}

void Synthwave_Scene_SetCamera(Scene* scene, Vector* position, Vector* target)
{
  if (scene->opaque == 0)
    return;

  Synthwave_Scene* s = (Synthwave_Scene*) (scene->opaque);
  
  s->cameraPosition  = *position;
  s->cameraTarget    = *target;

  Vector up;
  Vector_Set(&up, 0, 0, 1);
  Matrix_LookAt(&s->invView, position, target, &up);
  Matrix_Inverse(&s->invView, &s->invView);

  s->vpsOutOfDate    = true;
}

void Synthwave_Scene_SetProjection(Scene* scene, f32 fovY, f32 nearPlane, f32 farPlane)
{
  if (scene->opaque == 0)
    return;

  Synthwave_Scene* s = (Synthwave_Scene*) (scene->opaque);
  
  Matrix_Projection(&s->projection, s->fb.w, s->fb.h, fovY, nearPlane, farPlane);
  s->vpsOutOfDate = true;
}

void Synthwave_Scene_Submit(Scene* scene, Surface* surface)
{
  if (scene->opaque == 0 || surface->opaque == 0)
    return;

  Synthwave_Scene* s = (Synthwave_Scene*) (scene->opaque);
  Synthwave_Surface* g = (Synthwave_Surface*) (surface->opaque);
  
  Synthwave_FrameBuffer* fb = &s->fb;
  u32 fbSize = fb->w * fb->h;

  memset(fb->palette, 0, sizeof(u8) * fbSize);
  
  for(u32 i=0;i < fbSize;i++)
  {
    fb->depth[i] = 0.0f;
    fb->brightness[i] = 1.0f;
  }
  
  if (s->vpsOutOfDate)
  {
    Matrix_Identity(&s->viewProjectionScreen);
    Matrix_Multiply(&s->viewProjectionScreen, &s->invView,    &s->viewProjectionScreen);
    Matrix_Multiply(&s->viewProjectionScreen, &s->projection, &s->viewProjectionScreen);
    Matrix_Multiply(&s->viewProjectionScreen, &s->screen,     &s->viewProjectionScreen);

    s->vpsOutOfDate = false;
  }

  // @TODO Draw Cmds Here.
  
  u8* pixels = NULL;
  i32 pixelPitch = 0;

  SDL_LockTexture(fb->texture, NULL, (void*) &pixels, &pixelPitch);

  Colour* palette = $$.Palette.colour;

#define $$RENDER_TYPE 0

  for(u32 i=0, j=0;i < fbSize;i++,j+=4 /* RGBA for some reason */)
  {
#if $$RENDER_TYPE == 0
    u8 index = fb->palette[i];
    f32 brightness = fb->brightness[i];
    Colour* colour = &palette[index];
    pixels[j + 0] = (u8) (colour->r * brightness); // R
    pixels[j + 1] = (u8) (colour->g * brightness); // G
    pixels[j + 2] = (u8) (colour->b * brightness); // B
#elif $$RENDER_TYPE == 1
    u8 index = (u8) ((1.0f - fb->depth[i]) * 255.0f);
    pixels[j + 0] = index; // R
    pixels[j + 1] = index; // G
    pixels[j + 2] = index; // B
#elif $$RENDER_TYPE == 2
    u8 index = (u8) ((fb->brightness[i]) * 255.0f);
    pixels[j + 0] = index; // R
    pixels[j + 1] = index; // G
    pixels[j + 2] = index; // B
#endif
  }
  
  SDL_UnlockTexture(fb->texture);
  
  SDL_Rect dst;
  dst.x = fb->x;
  dst.y = fb->y;
  dst.w = fb->w;
  dst.h = fb->h;

  SDL_SetRenderTarget($$.Screen.renderer, g->texture);
  SDL_RenderCopy($$.Screen.renderer, fb->texture, NULL, &dst);
  SDL_SetRenderTarget($$.Screen.renderer, NULL);

}

void Synthwave_Scene_RenderMesh(Scene* scene, Matrix* worldMatrix, Mesh* mesh)
{
  
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
      
      $$.Time.accumulator -= $.Time.fixedMs;
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
  
  $.Palette.Add         = Synthwave_Palette_Add;
  $.Palette.AddRgb      = Synthwave_Palette_AddRgb;
  $.Palette.AddU32      = Synthwave_Palette_AddU32;
  $.Palette.GetCount    = Synthwave_Palette_GetCount;
  $.Palette.Reset       = Synthwave_Palette_Reset;
  $.Time.GetTicks       = Synthwave_Time_GetTicks;
  $.Mem.Heap            = Synthwave_Mem_Heap;
  $.Mem.Temp            = Synthwave_Mem_Temp;
  $.Debug.Log           = Synthwave_Debug_Log;
  $.Debug.LogF          = Synthwave_Debug_LogF;
  $.Text.Format         = Synthwave_Text_Format;
  $.Text.FormatHeap     = Synthwave_Text_FormatHeap;
  $.Scene.New           = Synthwave_Scene_New;
  $.Scene.NewXywh       = Synthwave_Scene_NewXywh;
  $.Scene.Delete        = Synthwave_Scene_Delete;
  $.Scene.Submit        = Synthwave_Scene_Submit;
  $.Scene.SetCamera     = Synthwave_Scene_SetCamera;
  $.Scene.SetProjection = Synthwave_Scene_SetProjection;
  $.Scene.RenderMesh    = Synthwave_Scene_RenderMesh;

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
  }

  $MaybeFunction($$.desc.Events.OnStart, /**/);

  Timer_Reset(&$$.Time.fixedLimiter);
  Timer_Reset(&$$.Time.frameTimer);
  Timer_Reset(&$$.Time.fpsTimer);
  
  if ($$.Screen.flags != 0)
  {
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

#if $IsWindows == 1
#include <windows.h>
#endif

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
