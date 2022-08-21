#pragma once

/*
 * STRONG_INLINE is a stronger version of the inline, 
 * using __forceinline on MSVC, always_inline on GCC/clang, and otherwise just use inline.
 */
#ifndef STRONG_INLINE
#if defined(_MSC_VER)
#define STRONG_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define STRONG_INLINE __attribute__((always_inline)) inline
#else
#define STRONG_INLINE inline
#endif
#endif