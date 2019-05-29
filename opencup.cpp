// problem B "Basirovich Maxim" from gp of moscow 2019
#include "clarkson.hpp"

#include <bits/stdc++.h>
using namespace std;
using ll = long long;
template<typename S, typename T>
void xmin(S&a, T const&b){if(b<a) a=b;}
template<typename S, typename T>
void xmax(S&a, T const&b){if(b>a) a=b;}

signed gen(int T){
    mt19937 rng(43151);
    auto get_rand = [&](int64_t l, int64_t r){
        return uniform_int_distribution<int64_t>(l, r)(rng);
    }; (void) get_rand;
    auto get_double = [&](double l, double r){
        return uniform_real_distribution<double>(l, r)(rng);
    }; (void) get_double;
    ofstream o("gen.txt");
    o << T << "\n";
    for(int cas=0;cas<T;++cas){
        int n = get_rand(1, 5000);
        int k = get_rand(2, 4);
        o << n << " " << k << "\n";
        for(int i=0;i<n;++i){
            o << get_rand(1, 100) << " ";
        }
        o << "\n";
        o << 0 << " ";
        for(int i=1;i<n;++i){
            o << get_rand(0, k-1) << " ";
        }
        o << "\n";


        o << "\n";
    }
    o << endl;
    o.close();
    return 0;
}

signed main()
{
    #ifdef LOCAL_RUN
    freopen("in.txt", "r", stdin);
    //freopen("out.txt", "w", stdout);
    cin.tie(0); cout.tie(0); ios_base::sync_with_stdio(false);
    int TTT; cin >> TTT;
	if(TTT < 0) return gen(-TTT);
	while(TTT--){
    #else
    cin.tie(0); cout.tie(0); ios_base::sync_with_stdio(false);
    #endif // LOCAL_RUN

    int n, k;
    cin >> n >> k;

    const int v = k;

    vector<vector<Num> > A;
    vector<Num> b, c(v);
    c[0] = 1;
    for(int i=1;i<k;++i){
        A.emplace_back(v);
        A.back()[i] = -1;
        b.emplace_back(0);
    }

    A.emplace_back(v, -1);
    A.back()[0] = 0;
    b.emplace_back(-1);

    vector<Num> pre(k, 0);

    vector<int> X(n), Y(n);
    for(auto &e:X) cin >> e;
    for(auto &e:Y) cin >> e;

    for(int i=0;i<n;++i){
        pre[Y[i]]+=X[i];
        A.emplace_back(v, 0);
        for(int j=0;j<k;++j){
            A.back()[j] = pre[j];
        }
        b.emplace_back(0);
    }
    auto res = dacin::lp::solve_clarkson_seidel(dacin::lp::Lp_Instance(move(A), move(b), move(c))); // 255 ms
    //auto res = dacin::lp::solve_clarkson_simplex(dacin::lp::Lp_Instance(move(A), move(b), move(c))); // 2.5s
    assert(res.is_bounded());
    if(!res.is_feasible()){
        cerr << "Infeasible\n";
        assert(0);
    } else {
        long double out = 0;
        out-=(long double)res.get_x()[0].to_double() / (long double)res.get_x().back().to_double();
        /*for(auto &e:res.get_x()){
            cerr << e << " ";
        }
        cerr << "\n";*/
        //cerr << out << "\n";
        assert(!isnan(out));
        //assert(0 <= out);
        cout << fixed << setprecision(20) << out << "\n";
    }





    #ifdef LOCAL_RUN
    }
    #endif // LOCAL_RUN
    return 0;
}

