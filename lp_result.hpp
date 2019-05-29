#ifndef LP_RESULT_HPP
#define LP_RESULT_HPP

#include "fraction.hpp"
#include "num.hpp"
#include "util.hpp"

namespace dacin{ namespace lp{

    enum class Lp_Status{
        INFEASIBLE, OPTIMAL, UNBOUNDED, ERROR
    };

    class Lp_Result{
    public:
        Lp_Result() : status(Lp_Status::ERROR), x(), ray(), objective(0) {}
        Lp_Result(Lp_Status status_, vector<Num> x_, vector<Num> ray_, Fraction objective_) : status(status_), x(std::move(x_)), ray(std::move(ray_)), objective(std::move(objective_)) {}

        static Lp_Result infeasible_result() {
            return Lp_Result(Lp_Status::INFEASIBLE, {}, {}, -Fraction::inf());
        }

        bool is_feasible() const {
            assert(status != Lp_Status::ERROR);
            return status != Lp_Status::INFEASIBLE;
        }
        bool is_bounded() const {
            assert(status != Lp_Status::ERROR);
            return status != Lp_Status::UNBOUNDED;
        }
        vector<Num> const& get_x() const { return x; }
        vector<Num> const& get_ray() const { return ray; }
        void set_x(vector<Num> const& x_){ x = x_; }
        void set_ray(vector<Num> const& ray_){ ray = ray_; }
        Fraction const& get_objective() const { return objective; }

        void reset_ray() {
            assert(is_feasible());
            if(!is_bounded()) status = Lp_Status::OPTIMAL;
            std::fill(ray.begin(), ray.end(), Num(0));
        }
        void recalc_objective(vector<Num> const&c){
            switch(status){
                case Lp_Status::OPTIMAL:
                    objective = Fraction(scal_affine(x, c), x.back());
                    break;
                case Lp_Status::UNBOUNDED:
                    objective = Fraction::inf();
                    break;
                case Lp_Status::INFEASIBLE:
                    objective = -Fraction::inf();
                    break;
                case Lp_Status::ERROR:
                    assert(0);
            }
        }
        void reduce_all(){ reduce_by_gcd(x); reduce_by_gcd(ray); }
        bool violates(vector<Num> const&row) const {
            if(!is_bounded()){
                const int sign = Num::cmp(scal(row, ray), Num(0));
                if(sign != 0) return sign == 1 ? true : false;
            }
            return scal(row, x) > Num(0);
        }

    private:
        Lp_Status status;
        vector<Num> x, ray;
        Fraction objective;
    };

} }
#endif // LP_RESULT_HPP
