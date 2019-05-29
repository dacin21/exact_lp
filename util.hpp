#ifndef UTIL_HPP
#define UTIL_HPP

#include <bits/stdc++.h>
#include "num.hpp"

namespace dacin{ namespace lp{
    using std::vector;
    using std::move;

    #ifdef LOCAL_RUN
    constexpr int seed = 918273741;
    std::mt19937 rng(918273741);
    #else
    std::mt19937 rng(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
    #endif // LOCAL_run
    void reset_seed(){
        #ifdef LOCAL_RUN
        rng = decltype(rng)(seed);
        #endif
    }

    Num scal_affine(vector<Num> const&a, vector<Num> const&b){
        size_t n = std::min(a.size(), b.size());
        assert(a.size() <= n+1 && b.size() <= n+1);
        Num ret(0);
        for(size_t i=0; i<n; ++i){
            ret += a[i]*b[i];
        }
        return ret;
    }
    Num scal(vector<Num> const&a, vector<Num> const&b){
        assert(a.size() == b.size());
        return scal_affine(a, b);
    }
    void reduce_by_gcd(vector<Num> &v){
        Num g(0);
        for(auto &e:v){
            g = Num::gcd(e%g, g);
        }
        g.set_neg(false);
        for(auto &e:v){
            e /= g;
        }
    }
    vector<Num> reduced(vector<Num> v){
        reduce_by_gcd(v);
        return v;
    }

} }
#endif // UTIL_HPP
