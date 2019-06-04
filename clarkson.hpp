#ifndef CLARKSON_HPP
#define CLARKSON_HPP

#include "fraction.hpp"
#include "lp_instance.hpp"
#include "lp_result.hpp"
#include "num.hpp"
#include "seidel.hpp"
#include "simplex.hpp"
#include "util.hpp"

namespace dacin{ namespace lp{

    using Backend = Lp_Result(*)(Lp_Instance);
    namespace detail{
        template<typename T>
        T randint(T l, T r){
            return std::uniform_int_distribution<T>(l, r)(rng);
        }
        Lp_Result clarkson_1(vector<vector<Num> > const&A, vector<Num> const&c, Backend backend){
            const int n = A.size(), d = c.size();
            const int k = 6*d*d;
            vector<int64_t> weight(n, 1);
            if(n <= k) return backend(std::move(Lp_Instance(A, c)));
            auto get_sublp = [&A, &c, &weight, n, d, k](){
                int64_t total_weight = std::accumulate(weight.begin(), weight.end(), int64_t{0});
                vector<int64_t> samples;
                while((int)samples.size() < k){
                    int64_t new_sample = randint<int64_t>(0, total_weight-1);
                    if(!std::count(samples.begin(), samples.end(), new_sample)){
                        samples.push_back(new_sample);
                    }
                }
                std::sort(samples.begin(), samples.end());
                int64_t weight_pre = 0;
                vector<vector<Num> > A_sub;
                for(int i=0, j=0; i<n && j<k; ++i){
                    weight_pre += weight[i];
                    if(samples[j] < weight_pre){
                        A_sub.push_back(A[i]);
                        ++j;
                    }
                }
                return Lp_Instance(std::move(A_sub), c);
            };
            vector<bool> is_violated(n);
            for(size_t iter = 0;;++iter){
                auto res = backend(std::move(get_sublp()));
                if(!res.is_feasible()){
                    return res;
                }
                int64_t violated_weight = 0;
                int64_t total_weight = 0;
                for(int i=0;i<n;++i){
                    is_violated[i] = res.violates(A[i]);
                    total_weight+=weight[i];
                    if(is_violated[i]){
                        violated_weight+=weight[i];
                    }
                }
                if(violated_weight == 0){
                    // std::cerr << "Clarkson 1 iter: " << iter << "\n";
                    return res;
                }
                if(violated_weight * 3 * d <= total_weight){
                    for(int i=0;i<n;++i){
                        if(is_violated[i]) weight[i]*=2;
                    }
                }
            }
        }
        Lp_Result clarkson_2(vector<vector<Num> > const&A, vector<Num> const&c, Backend backend){
            const int n = A.size(), d = c.size();
            if(n <= 9*d*d){
                return clarkson_1(A, c, backend);
            }
            const int root_n = llround(sqrt(n));
            const int k = d * root_n;
            vector<vector<Num> > A_sub;
            for(;;){
                const int s = A_sub.size();
                for(int i=0;i<k;++i){
                    A_sub.push_back(A[randint<int>(0, n-1)]);
                }
                auto res = clarkson_1(A_sub, c, backend);
                A_sub.erase(A_sub.begin()+s, A_sub.end());
                if(!res.is_feasible()){
                    return res;
                }
                vector<int> violators;
                for(int i=0;i<n;++i){
                    if(res.violates(A[i])) violators.push_back(i);
                }
                if(violators.empty()){
                    return res;
                }
                if((int)violators.size() <= 2*root_n){
                    for(auto &e:violators){
                        A_sub.push_back(A[e]);
                    }
                }
            }
        }
    }
    Lp_Result solve_clarkson(Lp_Instance const&lp, Backend backend){
        Lp_Result res = detail::clarkson_2(lp.get_A(), lp.get_c(), backend);
        return res;
    }

    template<bool move_to_front = false>
    Lp_Result solve_clarkson_seidel(Lp_Instance const&lp){
        return solve_clarkson(lp, solve_seidel<move_to_front>);
    }
    Lp_Result solve_clarkson_simplex(Lp_Instance const&lp){
        return solve_clarkson(lp, solve_simplex);
    }

} }
#endif // CLARKSON_HPP
