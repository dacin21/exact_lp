#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include "fraction.hpp"
#include "lp_result.hpp"
#include "lp_instance.hpp"
#include "util.hpp"

namespace dacin{ namespace lp{
    namespace detail{

        void reduce_tableau(vector<vector<Num> > &T, Num& scale){
            Num g = scale;
            for(auto const&e:T){
                for(auto &f:e){
                    g = Num::gcd(g, f);
                }
            }
            if(g != Num(1)){
                scale/=g;
                for(auto &e:T){
                    for(auto &f:e){
                        f/=g;
                    }
                }
            }
        }
        void pivot(vector<vector<Num> > &T, Num& scale, vector<int> &basic, vector<int> &nonbasic, int const enter, int const leave){
            const int X = T.size(), Y = T[0].size();
            const Num Drs_abs = Num::abs(T[leave][enter]);
            const int Drs_sign = T[leave][enter].sign();
            for(int i=0;i<X;++i) if(i != leave){
                for(int j=0;j<Y;++j) if(j != enter){
                    T[i][j] = T[i][j] * Drs_abs - T[leave][j] * T[i][enter] * Drs_sign;
                }
            }
            for(int j=0;j<Y;++j) if(j != enter){
                T[leave][j] *= scale * Drs_sign;
            }
            for(int i=0;i<X;++i) if(i != leave){
                T[i][enter] *= scale * -Drs_sign;
            }
            T[leave][enter] = scale * scale * Drs_sign;
            scale*= Drs_abs;
            std::swap(basic[leave], nonbasic[enter]);
            reduce_tableau(T, scale);
        }
        int run_phase(vector<vector<Num> > &T, Num&scale, vector<int> &basic, vector<int> &nonbasic, const int phase){
            const int n = T.size()-2, d = T[0].size()-2;
            const int x = phase==1 ? n+1 : n;
            for(;;){
                int enter = -1;
                // primal steepest edge with lexicographical tie breaking
                Fraction slope = Fraction::inf();
                for(int j=0;j<=d;++j){
                    if(phase==2 && nonbasic[j] == -1) continue;
                    Num norm_sq;
                    for(int i=0;i<=n;++i){
                        //norm_sq += T[i][j]*T[i][j];
                        Num::addmul_long(norm_sq, T[i][j], T[i][j]);
                    }
                    Fraction slope_j (T[x][j] * Num::abs(T[x][j]), norm_sq);
                    if(std::make_pair(slope_j, nonbasic[j]) < std::make_pair(slope, enter==-1 ? -1 : nonbasic[enter])){
                        enter = j;
                        slope = move(slope_j);
                    }
                }
                if(T[x][enter].sign() >= 0) return -1;
                int leave = -1;
                Fraction ratio = Fraction::inf();
                for(int i=0;i<n;++i){
                    if(T[i][enter].sign() > 0){
                        Fraction ratio_i = Fraction(T[i][d+1], T[i][enter]);
                        if(std::make_pair(ratio_i, basic[i])  < std::make_pair(ratio, leave==-1 ? -1 : basic[leave])){
                            leave = i;
                            ratio = move(ratio_i);
                        }
                    }
                }
                if(leave == -1) return enter;
                pivot(T, scale, basic, nonbasic, enter, leave);
            }
        }
        Lp_Result tableau_simplex(vector<vector<Num> > T){
            const int n = T.size()-2, d = T[0].size()-2;
            vector<int> basic(n); std::iota(basic.begin(), basic.end(), d);
            vector<int> nonbasic(d+1, -1); std::iota(nonbasic.begin(), prev(nonbasic.end()), 0);
            int leave = 0;
            for(int i=1;i<n;++i){
                if(T[i].back() < T[leave].back()) leave = i;
            }
            Num scale(1);
            if(T[leave][d+1].sign() < 0){
                pivot(T, scale, basic, nonbasic, d, leave);
                const int feasible_fail = run_phase(T, scale, basic, nonbasic, 1);
                if(feasible_fail != -1){
                    return Lp_Result::infeasible_result();
                }
                for(int i=0;i<n;++i) if(basic[i] == -1){
                    int enter = 0;
                    for(int j=1;j<=d;++j){
                        if(std::make_pair(T[i][j], nonbasic[j]) < std::make_pair(T[i][enter], nonbasic[enter])){
                            enter = j;
                        }
                    }
                    pivot(T, scale, basic, nonbasic, enter, i);
                }
            }
            const int bounded_fail = run_phase(T, scale, basic, nonbasic, 2);
            vector<Num> x(d+1, 0);
            x.back() = scale;
            for(int i=0;i<n;++i) if(basic[i] < d){
                x[basic[i]] = T[i].back();
            }
            if(bounded_fail != -1){
                const int bf = bounded_fail;
                // unbounded ray is needed for clarkson to work
                vector<Num> ray(d+1, 0);
                ray.back() = 0;
                for(int i=0;i<n;++i) if(basic[i] < d){
                    ray[basic[i]] = -T[i][bf];
                }
                if(nonbasic[bf] < d){
                    ray[nonbasic[bf]] = scale;
                }
                return Lp_Result(Lp_Status::UNBOUNDED, move(x), move(ray), Fraction::inf());
            }
            return Lp_Result(Lp_Status::OPTIMAL, move(x), {}, Fraction(T[n].back(), scale));
        }
    }

    vector<Num> transform_back(vector<Num> const&v){
        assert(v.size()%2 == 1);
        const int d0 = v.size()/2;
        vector<Num> ret(d0+1);
        for(int i=0;i<d0;++i){
            ret[i] = v[2*i] - v[2*i+1];
        }
        ret.back() = v.back();
        return ret;
    }

    Lp_Result solve_simplex(Lp_Instance lp){
        const int d0 = lp.d();
        const int d = 2*d0;
        const int n = lp.n();
        vector<vector<Num> > T(n+2, vector<Num>(d+2));
        for(int i=0;i<n;++i){
            for(int j=0;j<d0;++j){
                T[i][2*j] = lp.get_A()[i][j];
                T[i][2*j+1] = -lp.get_A()[i][j];
            }
            T[i][d] = Num(-1);
            T[i][d+1] = -lp.get_A()[i].back();
        }
        for(int j=0;j<d0;++j){
            T[n][2*j] = -lp.get_c()[j];
            T[n][2*j+1] = lp.get_c()[j];
        }
        T[n+1][d] = Num(1);
        Lp_Result ret = detail::tableau_simplex(move(T));
        if(ret.is_feasible()){
            ret.set_x(move(transform_back(ret.get_x())));
        }
        if(!ret.is_bounded()){
            ret.set_ray(move(transform_back(ret.get_ray())));
        }
        ret.reduce_all();
        return ret;
    }

} }
#endif // SIMPLEX_HPP
