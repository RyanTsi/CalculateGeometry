#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

using Real = double;
const Real eps = 1e-8;
const Real INF = 1e18;
const Real PI = acos(static_cast<Real>(-1));

struct Point;   // 点
struct Line;    // 直线
struct Seg;     // 线段
struct Polygon; // 多边形
struct Circle;  // 圆
typedef Point Vec;
using Points = std::vector<Point>;
using Lines  = std::vector<Line>;
// 位置关系
enum Position { CCW = 1, CW = -1, BACK = 2, FORNT = -2, ON = 0 };
// 包含关系
enum Contain { in = 2, on = 1, out = 0 };

// 判断实数的符号
inline int sgn(const Real& r) {
    return r < -eps ? -1 : r > eps ? 1 : 0;
}
// 两个实数判断大小
inline int cmp(const Real &a, const Real &b) {
    return abs(a - b) <= eps ? 0 : a > b ? 1 : -1;
}
// 判断两个实数是否相等
inline bool eq(const Real &a, const Real &b) {
    return cmp(a, b) == 0;
}
// 距离
Real distance(const Point& a,   const Point& b);
Real distance(const Line& l,    const Point& p);
Real distance(const Line& l,    const Line& m);
Real distance(const Line& l,    const Seg& s);
Real distance(const Seg& s,     const Point& p);
Real distance(const Seg& s,     const Seg& t);
// 相交判断
bool intersect(const Line& l,   const Point& p);
bool intersect(const Line& l,   const Line& m);
bool intersect(const Line& l,   const Seg& s);
bool intersect(const Seg& s,    const Point& p);
bool intersect(const Seg& s,    const Seg& t);
bool intersect(const Circle& O, const Line& l);
int  intersect(const Circle& O, const Seg& s);
int intersect(const Circle& c,  const Polygon& ps);
// 交点
Point crosspoint(const Line& l, const Line& m);
Points crosspoint(const Circle& c, const Line& l);
Points crosspoint(const Circle& c, const Seg& s);
Points crosspoint(const Circle& c1, const Circle& c2);

// 根据极坐标得到直角坐标系中的点
Point polar(const Real& r, const Real& theta);
// 点在直线上的投影点
Point proj(const Point& p, const Line& l);
// 点在直线上的对称点
Point refl(const Point& p, const Line& l);
// 中垂线
Line midline(const Point& p, const Point& q);
// 判断两直线平行
bool LPL(const Line& a, const Line& b);
// 判断两直线垂直
bool LOL(const Line& a, const Line& b);
// c 和 ab 的位置关系
int ccw(const Point& a, Point b, Point c);
// b 和 oa 的位置关系
int ccw(const Point& a,const Point& b);
// 海伦公式
Real Area(const Point& a, const Point& b, const Point& c);
// 以点 O 为原点进行极角排序 (CCW)，欧式距离为次关键字
void CCWwhitO(Point O, Points &ps);
// 平面最近点
Real closest_pair(Points ps);
// 判断点与多边形位置关系
int iscontain(const Point &a, const Points &ps);
// 求平面上所有点的凸包
Polygon Convexhell(Points &ps);
// 求凸多边形直径
Real Diameter(const Polygon& ps);
// 求凸多边形最小矩形覆盖 resp 返回矩形四角
Real RectangleArea(const Polygon& ps, Points& resp);
// 直线切凸包的左半边
Points Convex_cut(const Points &ps, const Line &l);
// 两个圆的公切线个数
int count_tangent(const Circle& c1, const Circle& c2);
// 两个圆的公切线
Lines common_tangent(const Circle& c1, const Circle& c2);
// 三角形外接圆
Circle circumcircle(Point a, Point b, const Point& c);
// 三角形内切圆
Circle incircle(const Point& a, const Point& b, const Point& c);
// 切点
Points tangent_to_circle(const Circle& c, const Point& p);
// 两圆公共面积
Real commonarea(Circle c1, Circle c2);
// 圆上两点 + 半径求圆心位置
Points center_given_radius(const Point& a, const Point& b, const Real& r);
// 平面最小圆覆盖
Circle get_min_circle(Points ps);
// 两个凸包的闵可夫斯基和
Polygon Minkowski(Polygon &A, Polygon &B);
// 判断点是否在一个凸多边形内
bool check_in(Point& p, Polygon& ps);

struct Point {
    Real x, y;
    Point() {}
    Point(Real a, Real b) : x(a), y(b) {}
    friend std::ostream& operator<<(std::ostream& os, const Point& v) {
        return os << v.x << " " << v.y;
    }
    friend std::istream& operator>>(std::istream& is, Point& v) {
        return is >> v.x >> v.y;
    }
    Point& operator+=(const Point& p) {
        x += p.x, y += p.y;
        return *this;
    }
    Point& operator-=(const Point& p) {
        x -= p.x, y -= p.y;
        return *this;
    }
    Point& operator*=(const Real& k) {
        x *= k, y *= k;
        return *this;
    }
    Point& operator/=(const Real& k) {
        x /= k, y /= k;
        return *this;
    }
    Point operator-() {
        return Point(-x, -y);
    }
    Point operator+(const Point& p) const {
        return Point(*this) += p;
    }
    Point operator-(const Point& p) const {
        return Point(*this) -= p;
    }
    Point operator*(const Real& k) const {
        return Point(*this) *= k;
    }
    Point operator/(const Real& k) const {
        return Point(*this) /= k;
    }
    bool operator==(const Point& p) {
        return eq(x, p.x) && eq(y, p.y);
    }
    bool operator!=(const Point& p) {
        return !(*this == p);
    }
    bool operator<(const Point& p) const { 
        return cmp(x, p.x) < 0 || (eq(x, p.x) && cmp(y, p.y) < 0);
    }
    bool operator>(const Point& p) const {
        return cmp(x, p.x) > 0 || (eq(x, p.x) && cmp(y, p.y) > 0);
    }
    // 叉积
    friend Real crs(const Point &a, const Point &b) {
        return a.x * b.y - b.x * a.y;
    }
    // 点积
    friend Real dot(const Point &a, const Point &b) {
        return a.x * b.x + a.y * b.y;
    }
    Real norm() const {
        return x * x + y * y;
    }
    // 模长
    Real abs() const {
        return std::hypot(x, y);
    }
    // 极角
    Real arg() const {
        return std::atan2(y, x);
    }
    // 逆时针旋转 90
    Point normal() const {
        return Point(-y, x);
    }
    Point unit() const {
        return *this / abs();
    }
    // 旋转
    Point rotate(Real theta) const {
        return Point(x * std::cos(theta) - y * std::sin(theta),
                        x * std::sin(theta) + y * std::cos(theta));
    }
};

struct Line {
    Point a, b;
    Line() {}
    Line(const Point &a, const Point &b) : a(a), b(b) {}
    Line(const Real& A, const Real& B, const Real& C) {  // Ax + By = c
        if (eq(A, 0)) {
            assert(!eq(B, 0));
            a = Point(0, C / B), b = Point(1, C / B);
        } else if (eq(B, 0)) {
            a = Point(C / A, 0), b = Point(C / A, 1);
        } else {
            a = Point(0, C / B), b = Point(C / A, 0);
        }
    }
    friend std::ostream& operator<<(std::ostream& os, const Line& l) {
        return os << l.a << " " << l.b << '\n';
    }
    friend std::istream& operator>>(std::istream& is, Line& l) {
        return is >> l.a >> l.b;
    }
    Point diff() const {
        return b - a;
    }
};

struct Seg : Line {
    Seg() {}
    Seg(Point a, Point b) : Line(a, b) {}
    Real len() const { return diff().abs(); }
};

struct Circle {
    Point c;
    Real r;
    Circle() {}
    Circle(const Point& c, const Real& r) : c(c), r(r) {}
    friend std::ostream& operator<<(std::ostream& os, const Circle& c) {
        return os << c.c << " " << c.r << '\n';
    }
    friend std::istream& operator>>(std::istream& is, Circle& c) {
        return is >> c.c >> c.r;
    }
};

Point polar(const Real& r, const Real& theta) {
    return Point(cos(theta), sin(theta)) * r;
}

Point proj(const Point& p, const Line& l) {
    Point v = l.diff().unit();
    return l.a + v * dot(v, p - l.a);
}

Point refl(const Point& p, const Line& l) {
    Point h = proj(p, l);
    return h + (h - p);
}

Line midline(const Point& p, const Point& q) {
    Point c = (p + q) * 0.5;
    Point v = (q - p).normal();
    return Line(c - v, c + v);
}

bool LPL(const Line& a, const Line& b) {
    return eq(crs(a.diff(), b.diff()), 0);
}

bool LOL(const Line& a, const Line& b) {
    return eq(dot(a.diff(), b.diff()), 0);
}

int ccw(const Point& a, Point b, Point c) {
    b -= a, c -= a;
    if(sgn(crs(b, c)) == +1) return CCW;
    if(sgn(crs(b, c)) == -1) return CW;
    if(sgn(dot(b, c)) == -1) return BACK;
    if(cmp(b.norm(), c.norm()) == -1) return FORNT;
    return ON;
}
int ccw(const Point& a,const Point& b) {
    return ccw(Point(0, 0), a, b);
}

Real distance(const Point& a, const Point& b) {
    return (a - b).abs();
}
Real distance(const Line& l, const Point& p) {
    return distance(p, proj(p, l));
}
Real distance(const Point &p, const Line& l) {
    return distance(l, p);
}
Real distance(const Line& l, const Line& m) {
    return intersect(l, m) ? 0 : distance(l, m.a);
}
Real distance(const Line& l, const Seg& s) {
    return intersect(l, s) ? 0 : std::min(distance(l, s.a), distance(l, s.b));
}
Real distance(const Seg& s, const Line& l) {
    return distance(l, s);
}
Real distance(const Seg& s, const Point& p) {
    Point h = proj(p, s);
    return intersect(s, h) ? distance(p, h) : std::min(distance(p, s.a), distance(p, s.b));
}
Real distance(const Point& p, const Seg& s) {
    return distance(p, s);
}
Real distance(const Seg& s, const Seg& t) {
    if (intersect(s, t)) return 0;
    return std::min({distance(s, t.a), distance(s, t.b), distance(t, s.a), distance(t, s.b)});
}

bool intersect(const Line& l, const Point& p) {
    return abs(ccw(l.a, l.b, p)) != 1;
}

bool intersect(const Line& l, const Line& m) {
    Real A = crs(l.diff(), m.diff()), B = crs(l.diff(), l.b - m.a);
    if (eq(A, 0) && eq(B, 0)) return true;      // check if equal
    return !LPL(l, m);
}

bool intersect(const Line& l, const Seg& s) {
    return sgn(crs(l.diff(), s.a - l.a)) * sgn(crs(l.diff(), s.b - l.a)) <= 0;
}

bool intersect(const Seg& s, const Point& p) {
    return ccw(s.a, s.b, p) == 0;
}

bool intersect(const Seg& s, const Seg& t) {
    return ccw(s.a, s.b, t.a) * ccw(s.a, s.b, t.b) <= 0 && ccw(t.a, t.b, s.a) * ccw(t.a, t.b, s.b) <= 0;
}

bool intersect(const Circle& O, const Line& l) {
    return cmp(O.r, distance(l, O.c)) >= 0;
}

bool intersect(const Line& l, const Circle& O) {
    return intersect(O, l);
}

int intersect(const Circle& c, const Seg& s) {
    Point h = proj(c.c, s);
    if (cmp(distance(c.c, h), c.r) > 0) return 0;
    Real d1 = (c.c - s.a).abs(), d2 = (c.c - s.b).abs();
    if (cmp(c.r, d1) >= 0 && cmp(c.r, d2) >= 0) return 0;
    if (cmp(c.r, d1) * cmp(c.r, d2) < 0) return 1;
    if (sgn(dot(s.a - h, s.b - h)) < 0) return 2;
    return 0;
}

Point crosspoint(const Line& l, const Line& m) {
    assert(intersect(l, m));
    Real A = crs(l.diff(), m.diff()), B = crs(l.diff(), l.b - m.a);
    if (eq(A, 0) && eq(B, 0)) return m.a;
    return m.a + m.diff() * B / A;
}

Point LTL(const Points &ps) {
    Point res = ps[0];
    for(auto p : ps) {
        if((p.y < res.y) || (eq(p.y, res.y) && p.x < res.x)) {
            res = p;
        }
    }
    return res;
}

Real Area(const Point& a, const Point& b, const Point& c) {
    return abs(crs(b - a, c - a)) / 2;
}

void CCWwhitO(Point O, Points &ps) {
    auto upp = [O](Point a) {
        if(a.y > O.y || (eq(a.y, O.y) && a.x >= O.x)) return true;
        return false;
    };
    auto cmp = [O, upp](Point a, Point b) {
        if(upp(a) != upp(b)) return upp(a);
        else {
            int s = sgn(crs(a - O, b - O));
            if(s == 0) {
                return (a-O).norm() < (b - O).norm();
            } else {
                return s > 0;
            }
        }
    };
    sort(ps.begin(), ps.end(), cmp);
}

Real closest_pair(Points ps) {
    int n = ps.size();
    if (n == 1) return 0;
    sort(ps.begin(), ps.end());
    auto cmp_y = [&](const Point& p, const Point& q) { return p.y < q.y; };
    std::vector<Point> cand(n);
    auto dfs = [&](auto self, int l, int r) -> Real {
        if (r - l <= 1) return INF;
        int mid = (l + r) >> 1;
        auto x_mid = ps[mid].x;
        auto res = std::min(self(self, l, mid), self(self, mid, r));
        inplace_merge(ps.begin() + l, ps.begin() + mid, ps.begin() + r, cmp_y);
        for (int i = l, cur = 0; i < r; i++) {
            if (abs(ps[i].x - x_mid) >= res) continue;
            for (int j = cur - 1; j >= 0; j--) {
                Point diff = {ps[i].x - cand[j].x, ps[i].y - cand[j].y};
                if (diff.y >= res) break;
                res = std::min(res, diff.abs());
            }
            cand[cur++] = ps[i];
        }
        return res;
    };
    return dfs(dfs, 0, n);
}

struct Polygon : Points {
    using Points::vector;
    Polygon() {}
    Polygon(int n) : Points(n) {}
    Real len () {
        Real res = 0;
        for(int i = 0; i < size(); i ++) {
            int j = (i + 1) % size();
            res += ((*this)[j] - (*this)[i]).abs();
        }
        return res;
    }
    Real Area() {
        Real res = 0;
        for(int i = 0; i < size(); i ++) {
            int j = (i + 1) % size();
            res += ((*this)[i].x + (*this)[j].x) * ((*this)[i].y - (*this)[j].y);
        }
        return abs(res) / 2;
    }
    bool is_convex(bool accept_on_segment = false) const {
        int n = size();
        for (int i = 0; i < n; i++) {
            if (accept_on_segment) {
                if (ccw((*this)[i], (*this)[(i + 1) % n], (*this)[(i + 2) % n]) == CW) {
                    return false;
                }
            } else {
                if (ccw((*this)[i], (*this)[(i + 1) % n], (*this)[(i + 2) % n]) != CCW) {
                    return false;
                }
            }
        }
        return true;
    }
};

int iscontain(const Point &a, const Points &ps) {
    int cnt = 0;
    for(int i = 0; i < ps.size(); i ++) {
        int j = (i + 1) % ps.size();
        Point u = ps[i] - a, v = ps[j] - a;
        if(u.y > v.y) std::swap(u, v);
        if(u.y <= 0 && v.y > 0 && crs(u, v) < 0) cnt ++;
        if(crs(u, v) == 0 && dot(u, v) <= 0) return on;
    }
    return cnt % 2 ? in : out;
}

Polygon Convexhell(Points &ps) {
    sort(ps.begin(), ps.end(), [&] (Point a, Point b) {
        return a < b;
    });
    int n = ps.size(), k = 0;
    Polygon res(2 * n);
    for(int i = 0; i < n; res[k ++] = ps[i ++]) {
        while(k >= 2 && ccw(res[k - 2], res[k - 1], ps[i]) != CCW) {
            k --;
        }
    }
    for(int i = n - 2, t = k + 1; i >= 0; res[k ++] = ps[i --]) {
        while(k >= t && ccw(res[k - 2], res[k - 1], ps[i]) != CCW) {
            k --;
        }
    }
    res.resize(k - 1);
    return res;
}

Real Diameter(const Polygon& ps) {
    int n = ps.size();
    if(n < 3) {
        return (ps[1] - ps[0]).abs();
    }
    Real res = 0;
    for(int i = 0, j = 1, k = 2 % n; i < n; ++ i, (++ j) %= n) {
        while(crs(ps[j] - ps[i], ps[(k + 1) % n] - ps[k]) >= 0) {
            (++ k) %= n;
        }
        res = std::max(res, std::max((ps[i] - ps[k]).abs(), (ps[j]- ps[k]).abs()));
    }
    return res;
}

Real RectangleArea(const Polygon& ps, Points& resp) {
    int n = ps.size();
    Real res = 1e18;
    int pr = 1, pl = 1;
    for(int i = 0, j = 1, k = 2 % n; i < n; ++ i, (++ j) %= n) {
        while(crs(ps[j] - ps[i], ps[(k + 1) % n] - ps[k]) >= 0) {
            (++ k) %= n;
        }
        while(dot(ps[j] - ps[i], ps[pr] - ps[i])
            <= dot(ps[j] - ps[i], ps[(pr + 1) % n] - ps[i])
        ) {
            (++ pr) %= n;
        }
        if(!i) pl = pr;
        while(dot(ps[j] - ps[i], ps[pl] - ps[i])
            >= dot(ps[j] - ps[i], ps[(pl + 1) % n] - ps[i])
        ) {
            (++ pl) %= n;
        }
        auto S= (dot(ps[j] - ps[i], ps[pr] - ps[i]) - 
                dot(ps[j] - ps[i], ps[pl] - ps[i])) /
                (ps[i]- ps[j]).abs() *
                distance(Line(ps[i], ps[j]), ps[k]);
        if(res > S) {
            res = S;
            auto l1 = Line(ps[i], ps[j]), l2 = Line(ps[k], ps[k] + l1.diff());
            resp[0] = proj(ps[pl], l1), resp[1] = proj(ps[pr], l1);
            resp[2] = proj(ps[pl], l2), resp[3] = proj(ps[pr], l2);
        }
    }
    return res;
}

Points Convex_cut(const Points &ps, const Line &l) {
    Points res;
    int n = ps.size();
    for(int i = 0, j = 1; i < n; ++ i, (++ j) %= n) {
        Real crs1 = crs(l.diff(), ps[i] - l.a);
        Real crs2 = crs(l.diff(), ps[j] - l.a);
        if(crs1 > 0.) res.push_back(ps[i]);
        if(crs1 > 0. && crs2 <= 0.) res.push_back(crosspoint(Line (ps[i], ps[j]), l));
        else if(crs1 <= 0. && crs2 > 0.) res.push_back(crosspoint(Line (ps[i], ps[j]), l));
    }
    return res;
}
    
int count_tangent(const Circle& c1, const Circle& c2) {
    Real r1 = c1.r, r2 = c2.r;
    if (r1 < r2) return count_tangent(c2, c1);
    Real d = distance(c1.c, c2.c);
    if (cmp(d, r1 + r2) > 0) return 4;
    if (cmp(d, r1 + r2) == 0) return 3;
    if (cmp(d, r1 - r2) > 0) return 2;
    if (cmp(d, r1 - r2) == 0) return 1;
    return 0;
}

Lines common_tangent(const Circle& c1, const Circle& c2) {
    if (c1.r < c2.r) return common_tangent(c2, c1);
    Lines res;
    Real g = distance(c1.c, c2.c);
    if (eq(g, 0)) return res;
    Point u = (c2.c - c1.c) / g, v = u.normal();
    for (int s : {-1, 1}) {
        Real h = (c1.r + c2.r * s) / g;
        if (eq(1 - h * h, 0))
            res.emplace_back(c1.c + u * c1.r, c1.c + (u + v) * c1.r);
        else if (cmp(1 - h * h, 0) > 0) {
            Point U = u * h, V = v * std::sqrt(1 - h * h);
            res.emplace_back(c1.c + (U + V) * c1.r, c2.c - (U + V) * c2.r * s);
            res.emplace_back(c1.c + (U - V) * c1.r, c2.c - (U - V) * c2.r * s);
        }
    }
    return res;
}

Circle circumcircle(Point a, Point b, const Point& c) {
    Point O = crosspoint(midline(a, c), midline(b, c));
    return Circle(O, distance(O, c));
}

Circle incircle(const Point& a, const Point& b, const Point& c) {
    Real A = (b - c).abs(), B = (c - a).abs(), C = (a - b).abs();
    Point O = (a * A + b * B + c * C) / (A + B + C);
    return Circle(O, distance(Seg(a, b), O));
}

Points crosspoint(const Circle& c, const Line& l) {
    Point h = proj(c.c, l);
    Real d = c.r * c.r - (c.c - h).norm();
    if (sgn(d) < 0) return {};
    if (sgn(d) == 0) return {h};
    Point v = l.diff().unit() * sqrt(d);
    return {h - v, h + v};
}

Points crosspoint(const Circle& c, const Seg& s) {
    int num = intersect(c, s);
    if (num == 0) return {};
    auto res = crosspoint(c, Line(s.a, s.b));
    if (num == 2) return res;
    if (sgn(dot(s.a - res[0], s.b - res[0])) > 0) std::swap(res[0], res[1]);
    return {res[0]};
}

Points crosspoint(const Circle& c1, const Circle& c2) {
    Real r1 = c1.r, r2 = c2.r;
    assert(!(eq(r1, r2) && sgn(distance(c1.c, c2.c) == 0)));
    if (r1 < r2) return crosspoint(c2, c1);
    Real d = distance(c1.c, c2.c);
    if (cmp(d, r1 + r2) > 0 || cmp(d, r1 - r2) < 0) return {};
    Real alpha = std::acos((r1 * r1 + d * d - r2 * r2) / (2 * r1 * d));
    Real theta = (c2.c - c1.c).arg();
    Point p = c1.c + polar(r1, theta + alpha);
    Point q = c1.c + polar(r1, theta - alpha);
    if (p == q) return {p};
    return {p, q};
}

Points tangent_to_circle(const Circle& c, const Point& p) {
    return crosspoint(c, Circle(p, sqrt((c.c - p).norm() - c.r * c.r)));
}

Real commonarea(Circle c1, Circle c2) {
    Real r1 = c1.r, r2 = c2.r;
    Real d = (c1.c - c2.c).abs();
    if (cmp(r1 + r2, d) <= 0) return 0;
    if (cmp(std::fabs(r1 - r2), d) >= 0) return PI * std::min(r1, r2) * std::min(r1, r2);
    Real res = 0;
    for (int _ = 0; _ < 2; _++) {
        r1 = c1.r, r2 = c2.r;
        Real cosine = (d * d + r1 * r1 - r2 * r2) / (2 * d * r1);
        Real theta = std::acos(cosine) * 2;
        res += (theta - std::sin(theta)) * r1 * r1 / 2;
        std::swap(c1, c2);
    }
    return res;
}

Points center_given_radius(const Point& a, const Point& b, const Real& r) {
    Point m = (b - a) * 0.5;
    Real d1 = m.abs();
    if (eq(d1, 0) || d1 > r) return {};
    Real d2 = sqrt(r * r - d1 * d1);
    Point n = m.normal() * d2 / d1;
    Point p = a + m - n, q = a + m + n;
    if (p == q) return {p};
    return {p, q};
}

int intersect(const Circle& c, const Polygon& ps) {
    Real d = 1e18;
    int n = ps.size();
    for(int i = 0; i < n; i ++) {
        d = std::min(d, distance(c.c, Seg(ps[i], ps[(i + 1) % n])));
    }
    int s = sgn(d - c.r);
    return s > 0 ? out : s < 0 ? in : on;
}

int CCP(Circle C, Point p) {
    Real dis = distance(C.c, p);
    int s = cmp(dis, C.r);
    return s == 0 ? on : s > 0 ? out : in;
}

Circle get_min_circle(Points ps) {
    int n = ps.size();
    random_shuffle(ps.begin(), ps.end()); // 随机化
    Circle res(Point(0, 0), 0);
    for(int i = 0; i < n; i ++) {
        // 圆是否包含点 i
        if(CCP(res, ps[i]) == out) {
            res = Circle(ps[i], 0);
            for(int j = 0; j < i; j ++) {
                if(CCP(res, ps[j]) == out) {
                    // 两点确定一个圆
                    res = Circle((ps[i] + ps[j]) / 2, distance(ps[i], ps[j]) / 2);
                    for(int k = 0; k < j; k ++) {
                        if(CCP(res, ps[k]) == out) {
                            // 三点确定一个圆
                            res = circumcircle(ps[i], ps[j], ps[k]);
                        }
                    }
                }
            }
        }
    }
    return res;
}

Polygon Minkowski(Polygon &A, Polygon &B) {
    int n = A.size(), m = B.size();
    Polygon res (n + m + 1);
    res[0] = A[0] + B[0];
    int k = 1, i = 0, j = 0;
    int cA = 0, cB = 0;
    while(cA < n && cB < m) {
        int nxi = (i + 1) % n, nxj = (j + 1) % m;
        Point x = A[nxi] - A[i];
        Point y = B[nxj] - B[j];
        if(ccw(x, y) == CCW) {
            res[k ++] = res[k - 1] + x;
            i = nxi, cA ++;
        } else if(ccw(x, y) == CW) {
            res[k ++] = res[k - 1] + y;
            j = nxj, cB ++;
        } else {
            res[k ++] = res[k - 1] + x + y;
            j = nxj, i = nxi;
            cA ++, cB ++;
        }
    }
    while(cA < n) {
        int nxi = (i + 1) % n;
        res[k ++] = res[k - 1] + A[nxi] - A[i];
        cA ++, i = nxi;
    }
    while(cB < m) {
        int nxj = (j + 1) % m;
        res[k ++] = res[k - 1] + B[nxj] - B[j];
        cB ++, j = nxj;
    }
    res.resize(k - 1);
    return res;
}

bool check_in(Point& p, Polygon& ps) {
    int n = ps.size();
    int O = 0, l = 1, r = n - 1;
    if(ccw(ps[O], ps[r], p) == CCW || ccw(ps[O], ps[l], p) == CW) {
        return false;
    }
    while(l + 1 < r) {
        int m = l + r >> 1;
        if(ccw(ps[O], ps[m], p) == CCW) {
            l = m;
        } else {
            r = m;
        }
    }
    if(ccw(ps[l], ps[r], p) == CW) {
        return false;
    } else {
        return true;
    }
}