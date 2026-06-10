// reads cases from stdin:
//   ncases
//   n ml md mh mw
//   l d h w   (x n)
// prints per case: p / k / o / ok on one line each
#include "gbp4d.h"
#include "gbp1d.h"

#include <armadillo>
#include <iostream>
#include <iomanip>

using namespace arma;
using namespace std;

int main() {
  int ncases;
  cin >> ncases;

  cout << std::setprecision(15);

  for (int t = 0; t < ncases; t++) {
    int n;
    double ml, md, mh, mw;
    cin >> n >> ml >> md >> mh >> mw;

    mat ldhw(4, n);
    for (int i = 0; i < n; i++) {
      cin >> ldhw(0, i) >> ldhw(1, i) >> ldhw(2, i) >> ldhw(3, i);
    }

    vec m = {ml, md, mh, mw};

    vec p = gbp4d_solver_dpp_prep_create_p(ldhw, m);
    gbp4d sn = gbp4d_solver_dpp(p, ldhw, m);

    { wall_clock timer; timer.tic(); vec p2 = gbp4d_solver_dpp_prep_create_p(ldhw, m); gbp4d sn2 = gbp4d_solver_dpp(p2, ldhw, m); cout << "case " << t << " cpp_ms: " << timer.toc()*1000 << endl; }
    cout << "p:"; for (uword i = 0; i < p.size(); i++) cout << " " << p(i); cout << endl;
    cout << "k:"; for (uword i = 0; i < sn.k.size(); i++) cout << " " << sn.k(i); cout << endl;
    cout << "o: " << sn.o << endl;
    cout << "ok: " << (sn.ok ? 1 : 0) << endl;
  }

  return 0;
}
