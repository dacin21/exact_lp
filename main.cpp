#ifdef RELEASE
#undef _GLIBCXX_DEBUG
#endif // RELEASE


#include "clarkson.hpp"
#include "seidel.hpp"
#include "tests.hpp"

#include <bits/stdc++.h>

using namespace std;

int main()
{
    /*Num a(0);
    auto b = -a;
    cout << (a<b) << " " << (a==b) << " " << a << " " << b << "\n";*/
    dacin::lp::run_tests();
    return 0;
}
