#ifndef FRACTION_HPP
#define FRACTION_HPP

#include "num.hpp"

namespace dacin{ namespace lp{

    class Fraction{
    public:
        Fraction() : a(0), b(1) {}
        explicit Fraction(Num x) : a(std::move(x)), b(1) {}
        Fraction(Num numerator, Num denominator) : a(std::move(numerator)), b(std::move(denominator)) { fix_sign(); }

        static Fraction inf(){return Fraction(1, 0);}

        Fraction& operator+=(Fraction const&o){
            a *= o.b;
            a += b*o.a;
            b *= o.b;
            return *this;
        }
        Fraction operator+(Fraction const&o) const {
            Fraction ret(*this);
            ret+=o;
            return ret;
        }
        Fraction& operator-=(Fraction const&o){
            a *= o.b;
            a += b*o.a;
            b *= o.b;
            return *this;
        }
        Fraction operator-(Fraction const&o) const {
            Fraction ret(*this);
            ret-=o;
            return ret;
        }

        Fraction& operator*=(Fraction const&o){
            a*=o.a;
            b*=o.b;
            return *this;
        }
        Fraction operator*(Fraction const&o) const {
            Fraction ret(*this);
            ret*=o;
            return ret;
        }
        Fraction& operator/=(Fraction const&o){
            a*=o.b;
            b*=o.a;
            return *this;
        }
        Fraction operator/(Fraction const&o) const {
            Fraction ret(*this);
            ret/=o;
            return ret;
        }
        Fraction operator-() const {
            return Fraction(-a, b);
        }

        int cmp(Fraction const&o) const {
            if(a.neg != o.a.neg){ // fix -inf < inf
                return a.neg ? -1 : 1;
            }
            return Num::cmp(a * o.b, b * o.a);
        }
        #define DECLARE_CMP_OP(op)\
        bool operator op (Fraction const&o) const{\
            return cmp(o) op 0;\
        }
        DECLARE_CMP_OP(<)
        DECLARE_CMP_OP(>)
        DECLARE_CMP_OP(<=)
        DECLARE_CMP_OP(>=)
        DECLARE_CMP_OP(==)
        DECLARE_CMP_OP(!=)
        #undef DECLARE_CMP_OP


        friend std::ostream& operator<<(std::ostream&o, Fraction const&f){
            return o << f.a << "/" << f.b;
        }

    private:
        void fix_sign(){
            if(a == Num(0)){
                if(b != Num(0)){
                    b = Num(1);
                }
            } else {
                if(b < Num(0)){
                    a.set_neg(-a.neg);
                    b.set_neg(false);
                }
            }
        }
        Num a, b;
    };

} }
#endif // FRACTION_HPP
