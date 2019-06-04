#ifndef TESTS_HPP
#define TESTS_HPP

#include <chrono>

#include "simplex.hpp"

#include "clarkson.hpp"
#include "fraction.hpp"
#include "lp_instance.hpp"
#include "lp_result.hpp"
#include "num.hpp"
#include "seidel.hpp"
#include "util.hpp"

namespace dacin{ namespace lp{
    using namespace std;
    using Solver = Lp_Result (*)(Lp_Instance);

    struct Timer{
        using Clock = chrono::high_resolution_clock;
        using Timepoint = Clock::time_point;
        Timer(string const&name_) : name(name_), start(Clock::now()) {}
        ~Timer() {
            Timepoint end = Clock::now();
            auto duration = chrono::duration_cast<chrono::nanoseconds>(end-start);
            cerr << "[TIMER] " << name << " took " << fixed << setprecision(6) << (duration.count() * 1e-9) << "\n";
        }
        string name;
        Timepoint start;
    };
    template<typename Fun, typename... Args>
    auto execute_timed(string const&name, Fun fun, Args... args) -> decltype(fun(args...)) {
        Timer timer(name);
        return fun(args...);
    }
    void test_from_file(string filename, Solver solver){
        reset_seed();
        cerr << "Test  file=" << filename << ", solver=" << (uintptr_t)solver << "\n";
        ifstream in(filename);
        in.exceptions(ios::badbit | ios::eofbit | ios::failbit);
        string s;
        in >> s;
        assert(s == "DACIN_LP");
        Lp_Instance lp;
        in >> lp;
        in >> s;
        assert(s == "SOL");
        auto read_num = [&in, &s](Num&f){
            in >> s;
            f = move(Num(s.c_str()));
        };
        Num numer, denom;
        read_num(numer);
        read_num(denom);
        Fraction objective(move(numer), move(denom));
        Lp_Status status;
        vector<Num> x;
        if(objective == -Fraction::inf()){
            status = Lp_Status::INFEASIBLE;
        } else if(objective == Fraction::inf()){
            status = Lp_Status::UNBOUNDED;
        } else {
            status = Lp_Status::OPTIMAL;
            x.resize(lp.d()+1);
            for(auto &e:x){
                read_num(e);
            }
        }

        Lp_Result sol = execute_timed("Lp", solver, lp);
        auto error = [&](string const& message){
            cerr << "[ERROR] " << message << "\n";
        };
        if(status == Lp_Status::INFEASIBLE){
            if(sol.is_feasible()){
                error("Lp is infeasible but solver reports feasible or unbounded.");
                return;
            }
        } else {
            if(!sol.is_feasible()){
                error("Lp is feasible but solver reports infeasible.");
                return;
            }
        }
        if(status == Lp_Status::UNBOUNDED){
            if(sol.is_bounded()){
                error("Lp is unbounded but solver report optimal.");
                return;
            }
        } else {
            if(!sol.is_bounded()){
                error("Lp is bounded but solver report unbounded.");
                return;
            }
        }
        if(status == Lp_Status::OPTIMAL){
            if(objective != sol.get_objective()){
                error("Lp objectives differ");
                cerr << objective << " vs " << sol.get_objective() << "\n";
            }
            if(reduced(x) != reduced(sol.get_x())){
                error("Solution differs (lp might be degenerate)");
            }
        }
        cerr << "    Test Passed.\n\n";
    }

    void run_tests_annulus(Solver solver){
        if(0){
        test_from_file("examples/enclosing_annulus_2_50.lp", solver);
        test_from_file("examples/enclosing_annulus_2_100.lp", solver);
        test_from_file("examples/enclosing_annulus_2_1000.lp", solver);
        test_from_file("examples/enclosing_annulus_3_50.lp", solver);
        test_from_file("examples/enclosing_annulus_4_50.lp", solver);
        }
        if(0){
        test_from_file("examples/enclosing_annulus_spherical_2_50.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_2_100.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_2_1000.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_3_50.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_4_50.lp", solver);
        }
        if(0){
        test_from_file("examples/enclosing_annulus_2_1000.lp", solver);
        test_from_file("examples/enclosing_annulus_2_10000.lp", solver);
        test_from_file("examples/enclosing_annulus_3_1000.lp", solver);
        test_from_file("examples/enclosing_annulus_3_10000.lp", solver);
        }
        if(0){
        test_from_file("examples/enclosing_annulus_spherical_2_1000.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_2_10000.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_3_1000.lp", solver);
        test_from_file("examples/enclosing_annulus_spherical_3_10000.lp", solver);
        }
        if(1){
            for(int it=0;it<5;++it){
                test_from_file("examples/enclosing_annulus_3_10000.lp", solver);
                test_from_file("examples/enclosing_annulus_spherical_3_10000.lp", solver);
            }
        }
    }
    void run_tests(){
        //run_tests_annulus(solve_seidel);
        run_tests_annulus([](Lp_Instance lp){return solve_clarkson(move(lp), solve_seidel);});
        run_tests_annulus([](Lp_Instance lp){return solve_clarkson(move(lp), solve_seidel<true>);});
        //run_tests_annulus(solve_simplex);
        //run_tests_annulus([](Lp_Instance lp){return solve_clarkson(move(lp), solve_simplex);});
    }

} }
#endif // TESTS_HPP
