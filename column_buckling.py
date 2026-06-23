import numpy as np
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa
    
    Return: העומס P בניוטון (float)
    """
    P_euler = (np.pi**2 * E * A) / (L / r)**2
    def f(P):
        if P == 0:
            return -sigma_allow
            
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        sec_val = 1.0 / np.cos(theta)
        sigma_max = (P / A) * (1 + (e * c / r**2) * sec_val)
        return sigma_max - sigma_allow

    P_critical = bisect(f, 0, 0.9999 * P_euler)
    
    return float(P_critical)
