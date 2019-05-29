#ifndef LP_INSTANCE_HPP
#define LP_INSTANCE_HPP

#include <cassert>
#include <istream>
#include <ostream>

#include "fraction.hpp"
#include "num.hpp"

namespace dacin{ namespace lp{

    class Lp_Instance{
    public:
        Lp_Instance() {}
        Lp_Instance(std::vector<std::vector<Num> > A_, std::vector<Num> c_) : A(std::move(A_)), c(std::move(c_)) { for(auto const&e:A) assert(e.size() == c.size()+1); }
        Lp_Instance(std::vector<std::vector<Num> > A_, std::vector<Num> b, std::vector<Num> c_) : A(std::move(A_)), c(std::move(c_)) {
            assert(A.size() == b.size());
            for(size_t i=0;i<A.size();++i){
                assert(A[i].size() == c.size());
                A[i].push_back(-b[i]);
            }
        }
        friend std::istream& operator>>(std::istream&in, Lp_Instance &lp){
            int n, d;
            in >> n >> d;
            std::vector<std::vector<Num> > A(n, std::vector<Num>(d+1));
            std::vector<Num> c(d);
            std::string s;
            auto read_num = [&in, &s](Num&f){
                in >> s;
                f = std::move(Num(s.c_str()));
            };
            for(auto &e:A) for(auto &f:e) read_num(f);
            for(auto &e:c) read_num(e);
            lp.set_A(std::move(A));
            lp.set_c(std::move(c));
            return in;
        }
        friend std::ostream& operator<<(std::ostream&o, Lp_Instance const&lp){
            o << lp.A.size() << " " << lp.c.size() << "\n";
            auto write_num = [&o](Num const&f){
                int g; assert(f.can_convert_to_int(&g));
                o << g;
            };
            for(auto &e:lp.A){
                for(auto &f:e) {
                    write_num(f);
                    o << " ";
                }
                o << "\n";
            }
            for(auto &e:lp.c){
                write_num(e);
                o << " ";
            }
            o << "\n";
            return o;
        }

        std::vector<std::vector<Num> > const& get_A() const { return A; }
        std::vector<Num> const& get_c() const { return c; }
        void set_A(std::vector<std::vector<Num> > A_) { A = std::move(A_); }
        void set_c(std::vector<Num> c_) { c = std::move(c_); }
        int n() const { return A.size(); }
        int d() const { return c.size(); }
    private:
        /*
         * maximize    c * x
         * subject to  A * (x|1) <= 0
         */
        std::vector<std::vector<Num> > A;
        std::vector<Num> c;
    };

} }
#endif // LP_INSTANCE_HPP
