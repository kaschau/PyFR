<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

#include <metal_stdlib>

using namespace metal;

// AoSoA macros
#define SOA_SZ ${soasz}
#define SOA_IX(a, v, nv) ((((a) / SOA_SZ)*(nv) + (v))*SOA_SZ + (a) % SOA_SZ)

// Typedefs
typedef ${pyfr.npdtype_to_ctype(fpdtype)} fpdtype_t;
typedef ${pyfr.npdtype_to_ctype(ixdtype)} ixdtype_t;
typedef bfloat bf16;

inline bf16 f32_to_bf16(float f) { return (bf16) f; }
inline float bf16_to_f32(bf16 b) { return (float) b; }

// Atomic helpers
% for aspace in ['device', 'threadgroup']:
% for op, op_pos, op_neg in [('min', 'min', 'max'), ('max', 'max', 'min')]:
inline void atomic_${op}_fpdtype(${aspace} fpdtype_t* addr, fpdtype_t val)
{
    union { float f; int i; uint u; } u; u.f = val;
    if (!signbit(val))
        atomic_fetch_${op_pos}_explicit((${aspace} atomic_int*) addr,
                                        u.i, memory_order_relaxed);
    else
        atomic_fetch_${op_neg}_explicit((${aspace} atomic_uint*) addr,
                                        u.u, memory_order_relaxed);
}
% endfor
inline void atomic_sum_fpdtype(${aspace} fpdtype_t* addr, fpdtype_t val)
{
    union { float f; uint u; } e, d;
    e.u = atomic_load_explicit((${aspace} atomic_uint*) addr,
                               memory_order_relaxed);
    do {
        d.f = e.f + val;
    } while (!atomic_compare_exchange_weak_explicit(
        (${aspace} atomic_uint*) addr, &e.u, d.u,
        memory_order_relaxed, memory_order_relaxed));
}
% endfor


// Thread/block helpers
#define PYFR_THREAD_ID tid_.x
#define PYFR_BLOCK_ID bid_.x
#define PYFR_BLOCK_ID_Y bid_.y
#define PYFR_SYNC_THREADS() threadgroup_barrier(mem_flags::mem_threadgroup)
#define PYFR_SYNC_GMEM_THREADS() threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device)
#define PYFR_SHARED threadgroup
#define PYFR_GMEM device
#define PYFR_LMEM threadgroup

// Tiled layout indexing for block-diagonal matrices
#define TILED_IX(e, r, c, block_sz, tile_sz) \
    ((e) * (block_sz) * (block_sz) + \
     ((r) / (tile_sz)) * ((block_sz) / (tile_sz)) * (tile_sz) * (tile_sz) + \
     ((c) / (tile_sz)) * (tile_sz) * (tile_sz) + \
     ((r) % (tile_sz)) * (tile_sz) + \
     ((c) % (tile_sz)))

<%def name="argmax_storage(ftype)">
<% nsg = -(-nthreads // sgsize) %>
    PYFR_SHARED ${ftype} argmaxv[${nsg}];
    PYFR_SHARED short argmaxi[${nsg}];
</%def>

<%def name="argmax_reduce(ftype, val, idx, dst)">
<% nsg = -(-nthreads // sgsize) %>
    {
        short amsg = PYFR_THREAD_ID / ${sgsize}, amsl = PYFR_THREAD_ID % ${sgsize};
        ${ftype} amv = simd_max(${val});
        int ami = simd_min((${val} == amv) ? ${idx} : 0x7fffffff);
        if (amsl == 0) { argmaxv[amsg] = amv; argmaxi[amsg] = ami; }
        PYFR_SYNC_THREADS();
        amv = argmaxv[0]; ami = argmaxi[0];
        for (int amw = 1; amw < ${nsg}; amw++)
            if (argmaxv[amw] > amv || (argmaxv[amw] == amv && argmaxi[amw] < ami)) { amv = argmaxv[amw]; ami = argmaxi[amw]; }
        if (PYFR_THREAD_ID == 0) ${dst} = ami;
        PYFR_SYNC_THREADS();
    }
</%def>

// FP-precise block support
#define PYFR_FP_PRECISE_BEGIN _Pragma("clang fp reassociate(off) contract(off)")

// expm1/log1p shims (absent from MSL); branch-free ports of ARM
// optimized-routines math/aarch64/advsimd/v_{expm1f,log1pf}_inline.h
// (MIT OR Apache-2.0 WITH LLVM-exception).  In lieu of the caller-side
// special-case handling upstream requires: expm1 clamps to [-87, 88]
// and log1p requires x > -1.  The pragma protects the ln2 hi/lo split
// from fast-math reassociation.
inline fpdtype_t expm1(fpdtype_t x)
{
    x = clamp(x, fpdtype_t(-87), fpdtype_t(88));

    // Reduce argument: f in [-ln2/2, ln2/2], i is exact
    fpdtype_t j = rint(x*fpdtype_t(0x1.715476p+0));
    int i = int(j);
    fpdtype_t f;
    {
        PYFR_FP_PRECISE_BEGIN
        f = fma(-j, fpdtype_t(0x1.62e4p-1), x);
        f = fma(-j, fpdtype_t(0x1.7f7d1cp-20), f);
    }

    // expm1(f) ~= f + f^2*P(f)
    fpdtype_t f2 = f*f, f4 = f2*f2;
    fpdtype_t p01 = fma(f, fpdtype_t(0x1.5554aep-3), fpdtype_t(0x1.fffffep-2));
    fpdtype_t p23 = fma(f, fpdtype_t(0x1.12287cp-7), fpdtype_t(0x1.555736p-5));
    fpdtype_t p = fma(f2, p23, p01);
    p = fma(f4, fpdtype_t(0x1.6b55a2p-10), p);
    p = fma(f2, p, f);

    // t = 2^i; expm1(x) ~= p*t + (t - 1)
    fpdtype_t t = as_type<fpdtype_t>((i << 23) + 0x3f800000);
    return fma(p, t, t - fpdtype_t(1));
}

inline fpdtype_t log1p(fpdtype_t x)
{
    // x + 1 = t*2^k with t = m + 1, m in [-0.25, 0.5]
    fpdtype_t m = x + fpdtype_t(1);
    uint ku = (as_type<uint>(m) - 0x3f400000u) & 0xff800000u;
    fpdtype_t s = as_type<fpdtype_t>(0x40800000u - ku);
    fpdtype_t ms = as_type<fpdtype_t>(as_type<uint>(x) - ku);
    ms = ms + fma(fpdtype_t(0.25), s, fpdtype_t(-1));

    // log(1 + ms) on [-0.25, 0.5], pairwise Horner
    fpdtype_t q = fma(ms, fpdtype_t(0x1.5555aap-2), fpdtype_t(-0.5));
    fpdtype_t m2 = ms*ms;
    fpdtype_t p67 = fma(ms, fpdtype_t(-0x1.6f0d5ep-5), fpdtype_t(0x1.abcb6p-4));
    fpdtype_t p45 = fma(ms, fpdtype_t(-0x1.0da91p-3), fpdtype_t(0x1.28a1f4p-3));
    fpdtype_t p23 = fma(ms, fpdtype_t(-0x1.54ef78p-3), fpdtype_t(0x1.99675cp-3));
    fpdtype_t p = fma(m2, p67, p45);
    p = fma(m2, p, p23);
    p = fma(ms, p, fpdtype_t(-0x1.000038p-2));
    p = m2*p;
    p = fma(m2, p, ms);
    p = fma(m2, q, p);

    // + k*ln2
    fpdtype_t sb = fpdtype_t(as_type<int>(ku))*fpdtype_t(0x1.0p-23);
    return fma(sb, fpdtype_t(0x1.62e43p-1), p);
}

<%def name="_kdecl(name, bounds)">kernel void ${name}</%def>
<%def name="_karg(intent, t, n)">\
% if intent == 'in':
device const ${t}* ${n}\
% elif intent == 'out':
device ${t}* ${n}\
% else:
constant ${t}& ${n}\
% endif
</%def>
<%def name="_kextra(body)"><%
    extra = []
    if 'PYFR_THREAD_ID' in body:
        extra.append('uint3 tid_ [[thread_position_in_threadgroup]]')
        extra.append('uint3 bid_ [[threadgroup_position_in_grid]]')
    if 'PYFR_GLOBAL_ID' in body:
        extra.append('uint2 ji [[thread_position_in_grid]]')
%>${', '.join(extra)}</%def>

${next.body()}
