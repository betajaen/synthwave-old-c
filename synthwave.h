//  synthwave.h
//  
//  Copyright (C) 2017 Robin Southern   https://github.com/betajaen/synthwave
//  
//  This software may be modified and distributed under the terms
//  of the MIT license.  See the LICENSE file for details.

#ifndef SYNTHWAVE_H
#define SYNTHWAVE_H

#if defined(_MSC_VER)
#define $IsWindows 1
#define $IsBrowser 0
#else
#define $IsWindows 0
#define $IsBrowser 1
#endif

#include <stdint.h>
#include <stdbool.h>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;
typedef float    f32;
typedef double   f64;

#define $PI     3.14159265358979323846264338327950288f
#define $TWO_PI 6.28318530717958647692528676655900576f

#define $Min(A, B)              (((A) < (B)) ? (A) : (B))
#define $Max(A, B)              (((A) > (B)) ? (A) : (B))
#define $Min3(A, B, C)          ($Min(A, $Min(B, C)))
#define $Max3(A, B, C)          ($Max(A, $Max(B, C)))
#define $Min4(A, B, C, D)       ($Min(A, $Min(B, $Min(C, D))))
#define $Max4(A, B, C, D)       ($Max(A, $Max(B, $Max(C, D))))
#define $Clamp(X, MIN, MAX)     (X > MAX ? MAX : X < MIN ? MIN : X)
#define $Squared(X)             ((X) * (X))
#define $Lerp(X, Y, T)          (X + T * (Y - X))

typedef enum {
  ST_Windowed,
  ST_Console,
  ST_RunOnce
} Synthwave_Type;

typedef struct
{
  u32 start, paused;
  u8  state;
} Timer;

typedef union
{
  struct { f32 x, y, z, w;         };
  struct { f32 forward, right, up; };
  f32  m[4];
} Vector;

typedef union
{
  struct { f32 x, y, z;           };
  struct { f32 roll, pitch, yaw;  };
  f32    m[4];
} Rotation;

typedef union
{
  f32    m[16];
  f32    M[4][4];
  Vector row[4];
} Matrix;

typedef struct
{
  bool quit;
  
  struct
  {
    i32 x, y;
    u32 width, height;
  } Screen;

  struct
  {
    u32 frameMs, fixedMs;
    f32 delta, fixed;
    u32 numFrames;

    u32 (*GetTicks)();
  } Time;

} Synthwave_Interface;

extern Synthwave_Interface $;

typedef struct
{
  Synthwave_Type type;
  
  struct
  {
    i32         x, y;
    u32         w, h, scale;
    const char* title;
  } Screen;
  
  struct
  {
    u32  fixedMs, frameMs;
  } Time;
  
  struct
  {
    void (*OnStart)();
    void (*OnQuit)();
    void (*OnFrame)();
    void (*OnFixedFrame)();
  } Events;

} Synthwave_Description;


inline f32 $Deg2Rad(f32 rad) { return ((rad) * 180.0f / $PI); }
inline f32 $Rad2Deg(f32 deg) { return ((deg) * $PI    / 180.0f); }
f32 $NormaliseDegrees(f32 deg);
f32 $NormaliseRadians(f32 rad);

inline void Vector2_Set(Vector* v, f32 x, f32 y)               { v->m[0] = x; v->m[1] = y; v->m[2] = 0; v->m[3] = 1.0f; }
inline void Vector_Set(Vector* v, f32 x, f32 y, f32 z)         { v->m[0] = x; v->m[1] = y; v->m[2] = z; v->m[3] = 1.0f; }
inline void Vector4_Set(Vector* v, f32 x, f32 y, f32 z, f32 w) { v->m[0] = x; v->m[1] = y; v->m[2] = z; v->m[3] = w;    }
inline void Vector2_Add(Vector* v, Vector* a, Vector* b)       { Vector2_Set(v, a->m[0] + b->m[0], a->m[1] + b->m[1]); }
inline void Vector_Add(Vector* v, Vector* a, Vector* b)        { Vector_Set(v, a->m[0] + b->m[0], a->m[1] + b->m[1], a->m[2] + b->m[2]); }
inline void Vector4_Add(Vector* v, Vector* a, Vector* b)       { Vector4_Set(v, a->m[0] + b->m[0], a->m[1] + b->m[1], a->m[2] + b->m[2], a->m[3] + b->m[3]); }
inline void Vector2_Add_s(Vector* v, Vector* a, f32 s)         { Vector2_Set(v, a->m[0] + s, a->m[1] + s); }
inline void Vector_Add_s(Vector* v, Vector* a, f32 s)          { Vector_Set(v, a->m[0] + s, a->m[1] + s, a->m[2] + s); }
inline void Vector4_Add_s(Vector* v, Vector* a, f32 s)         { Vector4_Set(v, a->m[0] + s, a->m[1] + s, a->m[2] + s, a->m[3] + s); } ;
inline void Vector2_Sub(Vector* v, Vector* a, Vector* b)       { Vector2_Set(v, a->m[0] - b->m[0], a->m[1] - b->m[1]); }
inline void Vector_Sub(Vector* v, Vector* a, Vector* b)        { Vector_Set(v, a->m[0] - b->m[0], a->m[1] - b->m[1], a->m[2] - b->m[2]); }
inline void Vector4_Sub(Vector* v, Vector* a, Vector* b)       { Vector4_Set(v, a->m[0] - b->m[0], a->m[1] - b->m[1], a->m[2] - b->m[2], a->m[3] - b->m[3]); }
inline void Vector2_Sub_s(Vector* v, Vector* a, f32 s)         { Vector2_Set(v, a->m[0] - s, a->m[1] - s); }
inline void Vector_Sub_s(Vector* v, Vector* a, f32 s)          { Vector_Set(v, a->m[0] - s, a->m[1] - s, a->m[2] - s); }
inline void Vector4_Sub_s(Vector* v, Vector* a, f32 s)         { Vector4_Set(v, a->m[0] - s, a->m[1] - s, a->m[2] - s, a->m[3] - s); } ;
inline void Vector2_Mul(Vector* v, Vector* a, Vector* b)       { Vector2_Set(v, a->m[0] * b->m[0], a->m[1] * b->m[1]); }
inline void Vector_Mul(Vector* v, Vector* a, Vector* b)        { Vector_Set(v, a->m[0] * b->m[0], a->m[1] * b->m[1], a->m[2] * b->m[2]); }
inline void Vector4_Mul(Vector* v, Vector* a, Vector* b)       { Vector4_Set(v, a->m[0] * b->m[0], a->m[1] * b->m[1], a->m[2] * b->m[2], a->m[3] * b->m[3]); }
inline void Vector2_Mul_s(Vector* v, Vector* a, f32 s)         { Vector2_Set(v, a->m[0] * s, a->m[1] * s); }
inline void Vector_Mul_s(Vector* v, Vector* a, f32 s)          { Vector_Set(v, a->m[0] * s, a->m[1] * s, a->m[2] * s); }
inline void Vector4_Mul_s(Vector* v, Vector* a, f32 s)         { Vector4_Set(v, a->m[0] * s, a->m[1] * s, a->m[2] * s, a->m[3] * s); } ;
inline void Vector2_Div(Vector* v, Vector* a, Vector* b)       { Vector2_Set(v, a->m[0] / b->m[0], a->m[1] / b->m[1]); }
inline void Vector_Div(Vector* v, Vector* a, Vector* b)        { Vector_Set(v, a->m[0] / b->m[0], a->m[1] / b->m[1], a->m[2] / b->m[2]); }
inline void Vector4_Div(Vector* v, Vector* a, Vector* b)       { Vector4_Set(v, a->m[0] / b->m[0], a->m[1] / b->m[1], a->m[2] / b->m[2], a->m[3] / b->m[3]); }
inline void Vector2_Div_s(Vector* v, Vector* a, f32 s)         { Vector2_Set(v, a->m[0] / s, a->m[1] / s); }
inline void Vector_Div_s(Vector* v, Vector* a, f32 s)          { Vector_Set(v, a->m[0] / s, a->m[1] / s, a->m[2] / s); }
inline void Vector4_Div_s(Vector* v, Vector* a, f32 s)         { Vector4_Set(v, a->m[0] / s, a->m[1] / s, a->m[2] / s, a->m[3] / s); } ;

inline f32  Vector2_Dot(Vector* a, Vector* b)                  { return a->m[0] * b->m[0] + a->m[1] * b->m[1]; }
inline f32  Vector_Dot(Vector* a, Vector* b)                   { return a->m[0] * b->m[0] + a->m[1] * b->m[1] + a->m[2] * b->m[2]; }
inline f32  Vector4_Dot(Vector* a, Vector* b)                  { return a->m[0] * b->m[0] + a->m[1] * b->m[1] + a->m[2] * b->m[2] + a->m[3] * b->m[3]; }

inline void Vector_Cross(Vector* v, Vector* a, Vector* b)      { v->m[0] = a->m[1] * b->m[2] - b->m[1] * a->m[2];
                                                                 v->m[1] = a->m[2] * b->m[0] - b->m[2] * a->m[0];
                                                                 v->m[2] = a->m[0] * b->m[1] - b->m[0] * a->m[1];
                                                               }

inline f32 Vector2_Length2(Vector* v)                          { return Vector2_Dot(v, v); }
inline f32 Vector_Length2(Vector* v)                           { return Vector_Dot(v, v); }
inline f32 Vector4_Length2(Vector* v)                          { return Vector4_Dot(v, v); }

f32 Vector2_Length(Vector* v);
f32 Vector_Length(Vector* v);
f32 Vector4_Length(Vector* v);

inline void Vector2_Normalise(Vector* v)                       { Vector2_Mul_s(v, v, (1.0f / Vector_Length(v))); }
inline void Vector_Normalise(Vector* v)                        { Vector_Mul_s(v, v, (1.0f / Vector_Length(v))); }
inline void Vector4_Normalise(Vector* v)                       { Vector4_Mul_s(v, v, (1.0f / Vector_Length(v))); }

void Vector_MultiplyMatrix(Vector* v, Matrix* m);
void Vector4_MultiplyMatrix(Vector* v, Matrix* m);

void Vector_Rotate(Vector* out, Vector* v, Rotation* r);
void Vector_Rotate_Pitch(Vector* out, Vector* v, f32 roll);
void Vector_Rotate_Roll(Vector* out, Vector* v, f32 pitch);
void Vector_Rotate_Yaw(Vector* out, Vector* v, f32 yaw);

inline void Matrix_Identity(Matrix* m)                         { m->m[0] = 1.0f; m->m[1] = 0.0f; m->m[2] = 0.0f; m->m[3] = 0.0f;
                                                                 m->m[4] = 0.0f; m->m[5] = 1.0f; m->m[6] = 0.0f; m->m[7] = 0.0f;
                                                                 m->m[8] = 0.0f; m->m[9] = 0.0f; m->m[10]= 1.0f; m->m[11]= 0.0f;
                                                                 m->m[12]= 0.0f; m->m[13]= 0.0f; m->m[14]= 0.0f; m->m[15]= 1.0f;
                                                               }

void Matrix_Multiply(Matrix* out, Matrix* a, Matrix* b);

void Matrix_Transpose(Matrix* out, Matrix* v);

void Matrix_Inverse(Matrix* out, Matrix* m);

void Matrix_Rotate_Roll(Matrix* out, f32 x);

void Matrix_Rotate_Pitch(Matrix* out, f32 y);

void Matrix_Rotate_Yaw(Matrix* out, f32 z);

void Matrix_LookAt(Matrix* m, Vector* position, Vector* target, Vector* up);

void Rotation_Normalize(Rotation* rotation);

inline void Rotation_Set(Rotation* r, f32 roll, f32 pitch, f32 yaw)         { r->roll = roll; r->pitch = pitch; r->yaw = yaw;
                                                                              Rotation_Normalize(r);
                                                                            }
inline void Rotation_Add(Rotation* r, Rotation* a, Rotation* b)             { Rotation_Set(r, a->roll + b->roll, a->pitch + b->pitch, a->yaw + b->yaw); }
inline void Rotation_Add_rpw(Rotation* o, Rotation* a, f32 r, f32 p, f32 y) { Rotation_Set(o, a->roll + r, a->pitch + p, a->yaw + y); }
inline void Rotation_Sub(Rotation* r, Rotation* a, Rotation* b)             { Rotation_Set(r, a->roll - b->roll, a->pitch - b->pitch, a->yaw - b->yaw); }
inline void Rotation_Sub_rpw(Rotation* o, Rotation* a, f32 r, f32 p, f32 y) { Rotation_Set(o, a->roll - r, a->pitch - p, a->yaw - y); }
inline void Rotation_Mul(Rotation* r, Rotation* a, Rotation* b)             { Rotation_Set(r, a->roll * b->roll, a->pitch * b->pitch, a->yaw * b->yaw); }
inline void Rotation_Mul_rpw(Rotation* o, Rotation* a, f32 r, f32 p, f32 y) { Rotation_Set(o, a->roll * r, a->pitch * p, a->yaw * y); }
inline void Rotation_Mul_s(Rotation* o, Rotation* a, f32 s)                 { Rotation_Set(o, a->roll * s, a->pitch * s, a->yaw * s); }
inline void Rotation_Div(Rotation* r, Rotation* a, Rotation* b)             { Rotation_Set(r, a->roll / b->roll, a->pitch / b->pitch, a->yaw / b->yaw); }
inline void Rotation_Div_rpw(Rotation* o, Rotation* a, f32 r, f32 p, f32 y) { Rotation_Set(o, a->roll / r, a->pitch / p, a->yaw / y); }
inline void Rotation_Div_s(Rotation* o, Rotation* a, f32 s)                 { Rotation_Set(o, a->roll / s, a->pitch / s, a->yaw / s); }

inline void Timer_Reset(Timer* timer)
{
  timer->start  = $.Time.GetTicks();
  timer->paused = 0;
  timer->state  = 1;
}

inline void Timer_Pause(Timer* timer)
{
  if (timer->state == 1)
  {
    timer->state   = 3;
    timer->paused  = $.Time.GetTicks() - timer->start;
    timer->start   = 0;
  }
}

inline void Timer_Unpause(Timer* timer)
{
  if (timer->state == 3)
  {
    timer->state   = 1;
    timer->start   = $.Time.GetTicks() - timer->paused;
    timer->paused  = 0;
  }
}

inline u32  Timer_GetTime(Timer* timer)
{
  if (timer->state == 1)
    return $.Time.GetTicks() - timer->start;
  else if (timer->state == 3)
    return timer->paused;
  return 0;
}

inline bool Timer_IsRunning(Timer* timer)
{
  return timer->state >= 1;
}

inline bool Timer_IsPaused(Timer* timer)
{
  return timer->state == 3;
}

extern void $Setup(Synthwave_Description* d);

#endif
