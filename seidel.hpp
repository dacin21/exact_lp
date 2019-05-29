#ifndef SEIDEL_HPP
#define SEIDEL_HPP

#include "fraction.hpp"
#include "lp_instance.hpp"
#include "lp_result.hpp"
#include "num.hpp"

namespace dacin{ namespace lp{
    namespace detail{

        std::pair<std::vector<Num>, int> make_projection(std::vector<Num> constraint){
            const int k = std::find_if(constraint.begin(), constraint.end(), [](Num const&x){ return x.sign() != 0; }) - constraint.begin();
            assert(k+1 < (int)constraint.size());
            if(constraint[k].sign() == -1){
                for(auto &e:constraint){
                    e = -e;
                }
            }
            return std::make_pair(std::move(constraint), k);
        }
        std::vector<Num> project_down(std::vector<Num> const&vec, std::vector<Num> const&plane, int const i){
            const size_t n = vec.size();
            assert(n <= plane.size() && plane.size() <= n+1);
            assert(plane[i].sign() > 0);
            std::vector<Num> ret(n-1);
            for(int j=0;j<i;++j) ret[j] = vec[j]*plane[i] - vec[i]*plane[j];
            for(int j=i+1;j<(int)n;++j) ret[j-1] = vec[j]*plane[i] - vec[i]*plane[j];
            return ret;
        }
        std::vector<Num> project_up(std::vector<Num> const&vec, std::vector<Num> const&plane, int const i){
            const size_t n = vec.size();
            assert(plane.size() == n+1);
            assert(plane[i].sign() > 0);
            std::vector<Num> ret(n+1);
            for(int j=0;j<i;++j){
                ret[j] = vec[j] * plane[i];
                ret[i] -= vec[j] * plane[j];
            }
            for(int j=i;j<(int)n;++j){
                ret[j+1] = vec[j] * plane[i];
                ret[i] -= vec[j] * plane[j+1];
            }
            return ret;
        }

        Lp_Result seidel_rec(Lp_Instance lp){
            const int n = lp.n(), d = lp.d();
            if(d == 0){
                Lp_Result ret(Lp_Status::OPTIMAL, {Num(1)}, {Num(0)}, Fraction(0));
                for(auto const&e:lp.get_A()){
                    if(ret.violates(e)) return Lp_Result::infeasible_result();
                }
                return ret;
            }
            auto get_base_result = [&lp, &n, &d](){
                std::vector<Num> x(d+1), ray(d+1);
                x.back() = Num(1);
                for(int i=0;i<d;++i){
                    ray[i] = Num(lp.get_c()[i].sign());
                }
                bool obj_unbounded = scal_affine(lp.get_c(), ray).sign() > 0;
                Lp_Result ret(obj_unbounded ? Lp_Status::UNBOUNDED : Lp_Status::OPTIMAL, move(x), move(ray), Fraction(0));
                return ret;
            };
            if(d == 1){
                // solve in deterministic O(n)
                Lp_Result ret = get_base_result();
                for(auto &e:lp.get_A()){
                    if(ret.violates(e)){
                        switch(e[0].sign()){
                            case -1:
                                ret.reset_ray();
                                ret.set_x({e[1], -e[0]});
                                break;
                            case 1:
                                ret.reset_ray();
                                ret.set_x({-e[1], e[0]});
                                break;
                            case 0:
                                return Lp_Result::infeasible_result();
                        }
                        ret.recalc_objective(lp.get_c());
                    }
                }
                for(auto &e:lp.get_A()){
                    if(ret.violates(e)){
                        return Lp_Result::infeasible_result();
                    }
                }
                return ret;
            } else {
                {
                    auto &A = const_cast<std::vector<std::vector<Num>>&>(lp.get_A());
                    std::shuffle(A.begin(), A.end(), rng);
                }
                Lp_Result ret = get_base_result();
                for(int i=0;i<n;++i){
                    auto const&e = lp.get_A()[i];
                    if(ret.violates(e)){
                        // project down, recurse, project up
                        auto projection = make_projection(e);
                        auto const plane = projection.first;
                        const int k = projection.second;
                        std::vector<std::vector<Num> > A_sub;
                        for(int j=0;j<i;++j){
                            A_sub.push_back(std::move(project_down(lp.get_A()[j], plane, k)));
                        }
                        std::vector<Num> c_sub = project_down(lp.get_c(), plane, k);
                        Lp_Result sub_result = seidel_rec(std::move(Lp_Instance(std::move(A_sub), std::move(c_sub))));
                        if(!sub_result.is_feasible()) return sub_result;
                        ret = sub_result;
                        ret.set_x(project_up(sub_result.get_x(), plane, k));
                        ret.set_ray(project_up(sub_result.get_ray(), plane, k));
                    }
                }
                return ret;
            }
        }
    }
    Lp_Result solve_seidel(Lp_Instance lp){
        auto const c = lp.get_c();
        auto res = detail::seidel_rec(std::move(lp));
        res.reduce_all();
        res.recalc_objective(c);
        return res;
    }

} }
#endif // SEIDEL_HPP
