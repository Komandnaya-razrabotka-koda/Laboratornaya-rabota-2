
# -*- coding: utf-8 -*-
"""
JM (Jelinski–Moranda) calculator — методичка + красивый вывод
"""

from math import isfinite, log

def _sum_ix(X):
    return sum((i)*x for i, x in enumerate(X, start=1))

def _f_B(B, X):
    n = len(X)
    S1 = sum(X)
    S2 = _sum_ix(X)
    if B <= n-1:
        return float('inf')
    denom_right = (B+1)*S1 - S2
    if denom_right <= 0:
        return float('inf')
    left = sum(1.0/(B - i + 1.0) for i in range(1, n+1))
    right = n*S1 / denom_right
    return left - right

def _solve_B(X, Bmax_initial=1e6, max_expand=8, tol=1e-12):
    n = len(X)
    Bmin = (n-1) + 1e-8
    fmin = _f_B(Bmin, X)
    Bmax = Bmax_initial
    fmax = _f_B(Bmax, X)
    expansions = 0
    while (not isfinite(fmax) or fmin*fmax > 0) and expansions < max_expand:
        Bmax *= 10.0
        fmax = _f_B(Bmax, X)
        expansions += 1
    info = {"Bmin": Bmin, "f(Bmin)": fmin, "Bmax": Bmax, "f(Bmax)": fmax, "expansions": expansions}
    if not (isfinite(fmin) and isfinite(fmax) and fmin*fmax < 0):
        info["diverges"] = True
        return None, info
    a, b, fa, fb = Bmin, Bmax, fmin, fmax
    for _ in range(200):
        m = 0.5*(a+b)
        fm = _f_B(m, X)
        if not isfinite(fm):
            a, fa = m, fm
            continue
        if abs(fm) < tol or (b-a) < 1e-9:
            return m, info
        if fa*fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5*(a+b), info

def _harmonic(m: int) -> float:
    if m <= 0:
        return 0.0
    if m < 100000:
        return sum(1.0/j for j in range(1, m+1))
    gamma = 0.5772156649015329
    return log(m) + gamma + 1.0/(2*m) - 1.0/(12*m*m)

def solve_jm_sheet(X):
    n = len(X)
    S1 = sum(X)
    S2 = _sum_ix(X)
    res = {"X": X, "n": n, "sum_X": S1, "sum_iXi": S2}
    B, info = _solve_B(X)
    res["solver_info"] = info
    if B is None:
        res["B"] = None
        return res
    res["B"] = B
    D = (B + 1) * S1 - S2
    K = n / D
    res["K"] = K
    res["denominator"] = D
    Rn = K * (B - n)
    res["R_n"] = Rn
    res["X_next"] = 1.0 / Rn if Rn > 0 else float('inf')
    remaining = B - n
    res["remaining_errors_est"] = remaining
    m = round(remaining)
    res["remaining_errors_rounded"] = m
    res["T_end"] = (1.0 / K) * _harmonic(m) if m > 0 else 0.0
    return res

if __name__ == "__main__":
    X = [5,4,11,13,6,2,5,5,8,7,1,4,2,7,6,2,3,1,4,78,25,10,7,16,3,1,2]
    res = solve_jm_sheet(X)
    print("\n=== Расчёт по модели Джелинского–Моранды ===\n")
    if res["B"] is None:
        print("❌ Не удалось найти конечное значение B (расходимость).")
    else:
        print(f"Число ошибок (наблюдено): n = {res['n']}")
        print(f"Сумма интервалов: ΣXᵢ = {res['sum_X']} ч")
        print(f"Взвешенная сумма: Σ(i·Xᵢ) = {res['sum_iXi']}")
        print(f"\nНайденное значение B: {res['B']:.6f} (общее число ошибок)")
        print(f"Осталось ошибок ≈ {res['remaining_errors_est']:.2f} "
              f"(округлённо {res['remaining_errors_rounded']})")
        print(f"Знаменатель (B+1)ΣXᵢ - Σ(i·Xᵢ) = {res['denominator']:.4f}")
        print(f"K = {res['K']:.10f} 1/ч")
        print(f"Интенсивность после {res['n']}-й ошибки: R(tₙ) = {res['R_n']:.6f} 1/ч")
        print(f"Ожидаемое время до следующей ошибки: Xₙ₊₁ = {res['X_next']:.2f} ч")
        print(f"Время до окончания тестирования: tₖ ≈ {res['T_end']:.2f} ч "
              f"(≈ {res['T_end']/24:.1f} суток)")
    print("\n=== Конец отчёта ===\n")
