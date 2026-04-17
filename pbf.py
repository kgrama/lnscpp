"""
PBF: Parameterized Bounded Format (Table-Free Edition)

All transcendental operations via continued fractions.
No lookup tables. Pure integer ALU + CF evaluation.
"""

from math import pi, log10
import random


def clz(x: int, n_bits: int) -> int:
    """Count leading zeros in n_bits-wide integer."""
    if x == 0:
        return n_bits
    return n_bits - x.bit_length()


def kdiv(a: int, b: int) -> int:
    """Kuttaka division: round to nearest."""
    if b == 0:
        return a
    if b > 0:
        return (a + b // 2) // b
    else:
        return (a - b // 2) // b


class PBF:
    """
    Parameterized Bounded Format with CF-based arithmetic.
    
    No tables — all logarithms and exponentials computed
    via continued fractions at runtime.
    """
    
    # Fixed-point scale (Q32)
    Q: int = 1 << 32
    
    # Precomputed constants in Q format
    LN2: int = 2977044472  # ln(2) * 2^32
    
    def __init__(self, n_bits: int = 8, v_min: float = 1e-4, v_max: float = 1e4):
        self.n_bits = n_bits
        self.v_min = v_min
        self.v_max = v_max
        
        # Encoding structure
        self.max_code = (1 << n_bits) - 1
        self.mid = 1 << (n_bits - 1)
        self.n_levels = self.mid - 1
        
        # Scale factor: maps [0, n_levels] to [ln(v_min), ln(v_max)]
        ln_range = self._ln_cf(int(v_max / v_min * self.Q))
        self.scale = ln_range // self.n_levels
    
    # ===== CONTINUED FRACTIONS =====
    
    def _ln1p_cf(self, x: int) -> int:
        """
        ln(1 + x/Q) * Q via continued fraction.
        """
        if x == 0:
            return 0
        if x <= -self.Q:
            return -self.Q * 100
        
        Q = self.Q
        h = 8 * Q
        h = 7 * Q + kdiv(16 * x * Q, h)
        h = 6 * Q + kdiv(9 * x * Q, h)
        h = 5 * Q + kdiv(9 * x * Q, h)
        h = 4 * Q + kdiv(4 * x * Q, h)
        h = 3 * Q + kdiv(4 * x * Q, h)
        h = 2 * Q + kdiv(x * Q, h)
        h = 1 * Q + kdiv(x * Q, h)
        return kdiv(x * Q, h)
    
    def _ln_cf(self, x: int) -> int:
        """
        ln(x/Q) * Q via range reduction + ln1p.
        """
        if x <= 0:
            return -self.Q * 100
        if x == self.Q:
            return 0
        
        Q = self.Q
        
        k = 0
        m = x
        while m >= 2 * Q:
            m //= 2
            k += 1
        while m < Q // 2:
            m *= 2
            k -= 1
        
        ln_m = self._ln1p_cf(m - Q)
        return k * self.LN2 + ln_m
    
    def _expm1_cf(self, x: int) -> int:
        """
        (exp(x/Q) - 1) * Q via continued fraction.
        """
        if x == 0:
            return 0
        
        Q = self.Q
        
        if x > 20 * Q:
            return self.Q * (1 << 20)
        if x < -20 * Q:
            return -Q
        
        h = 2 * Q
        h = 7 * Q + kdiv(-x * Q, h)
        h = 2 * Q + kdiv(x * Q, h)
        h = 5 * Q + kdiv(-x * Q, h)
        h = 2 * Q + kdiv(x * Q, h)
        h = 3 * Q + kdiv(-x * Q, h)
        h = 2 * Q + kdiv(x * Q, h)
        h = Q + kdiv(-x * Q, h)
        
        return kdiv(x * Q, h)
    
    def _exp_cf(self, x: int) -> int:
        """exp(x/Q) * Q via range reduction + expm1."""
        if x == 0:
            return self.Q
        
        Q = self.Q
        
        k = kdiv(x, self.LN2)
        r = x - k * self.LN2
        
        exp_r = Q + self._expm1_cf(r)
        
        if k >= 0:
            return exp_r << k
        else:
            return exp_r >> (-k)
    
    def _softplus_cf(self, d: int) -> int:
        """
        ln(1 + exp(d/Q)) * Q — the softplus function.
        """
        Q = self.Q
        
        if d > 10 * Q:
            return d
        if d < -10 * Q:
            return self._exp_cf(d)
        
        if d > 0:
            exp_neg = self._exp_cf(-d)
            return d + self._ln1p_cf(exp_neg)
        else:
            # ln(1 + exp(d)) where exp(d) <= 1
            # exp_d is in Q format: exp_d/Q = exp(d/Q)
            exp_d = self._exp_cf(d)
            return self._ln1p_cf(exp_d)  # ln(1 + exp_d/Q) * Q
    
    def _softminus_cf(self, d: int) -> int:
        """
        ln(1 - exp(d/Q)) * Q for d < 0.
        """
        Q = self.Q
        
        if d >= 0:
            return -100 * Q
        if d < -10 * Q:
            return -self._exp_cf(d)
        
        exp_d = self._exp_cf(d)
        return self._ln1p_cf(-exp_d)
    
    # ===== ENCODING / DECODING =====
    
    def encode(self, value: float) -> int:
        """Encode float to PBF byte."""
        if value == 0:
            return 0
        if abs(value) == float('inf'):
            return self.max_code
        
        sign = 1 if value > 0 else -1
        mag = max(self.v_min, min(abs(value), self.v_max))
        
        ratio_q = int(mag / self.v_min * self.Q)
        ln_ratio = self._ln_cf(ratio_q)
        
        level = kdiv(ln_ratio, self.scale)
        level = max(0, min(level, self.n_levels - 1))
        
        if sign > 0:
            return self.mid + level
        else:
            return self.mid - 1 - level
    
    def decode(self, code: int) -> float:
        """Decode PBF byte to float."""
        if code == 0:
            return 0.0
        if code == self.max_code:
            return float('inf')
        
        if code >= self.mid:
            level = code - self.mid
            sign = 1.0
        else:
            level = self.mid - 1 - code
            sign = -1.0
        
        ln_ratio = level * self.scale
        exp_ratio = self._exp_cf(ln_ratio)
        
        return sign * self.v_min * exp_ratio / self.Q
    
    # ===== ARITHMETIC =====
    
    def neg(self, a: int) -> int:
        """Negate: ones' complement."""
        if a == 0 or a == self.max_code:
            return a
        return self.max_code - a
    
    def mul(self, a: int, b: int) -> int:
        """Multiply: add log magnitudes with v_min offset."""
        if a == 0 or b == 0:
            return 0
        if a == self.max_code or b == self.max_code:
            return self.max_code
        
        a_neg = a < self.mid
        b_neg = b < self.mid
        a_mag = (self.mid - 1 - a) if a_neg else (a - self.mid)
        b_mag = (self.mid - 1 - b) if b_neg else (b - self.mid)
        
        # a_val = v_min * r^a_mag, b_val = v_min * r^b_mag
        # a * b = v_min² * r^(a_mag + b_mag)
        # To express as v_min * r^c_mag: c_mag = a_mag + b_mag + ln(v_min)/scale
        # ln(v_min) in Q format:
        ln_vmin = self._ln_cf(int(self.v_min * self.Q))
        offset = kdiv(ln_vmin, self.scale)
        
        c_mag = a_mag + b_mag + offset
        if c_mag >= self.n_levels:
            return self.max_code
        if c_mag < 0:
            return 0
        
        c_neg = a_neg != b_neg
        return (self.mid - 1 - c_mag) if c_neg else (self.mid + c_mag)
    
    def div(self, a: int, b: int) -> int:
        """Divide: subtract log magnitudes with v_min offset."""
        if a == 0:
            return 0
        if b == 0 or a == self.max_code:
            return self.max_code
        if b == self.max_code:
            return 0
        
        a_neg = a < self.mid
        b_neg = b < self.mid
        a_mag = (self.mid - 1 - a) if a_neg else (a - self.mid)
        b_mag = (self.mid - 1 - b) if b_neg else (b - self.mid)
        
        # a / b = (v_min * r^a_mag) / (v_min * r^b_mag) = r^(a_mag - b_mag)
        # To express as v_min * r^c_mag: c_mag = a_mag - b_mag - ln(v_min)/scale
        ln_vmin = self._ln_cf(int(self.v_min * self.Q))
        offset = kdiv(ln_vmin, self.scale)
        
        c_mag = a_mag - b_mag - offset
        if c_mag < 0:
            return 0
        if c_mag >= self.n_levels:
            return self.max_code
        
        c_neg = a_neg != b_neg
        return (self.mid - 1 - c_mag) if c_neg else (self.mid + c_mag)
    
    def add(self, a: int, b: int) -> int:
        """Add using continued fractions."""
        if a == 0:
            return b
        if b == 0:
            return a
        if a == self.max_code or b == self.max_code:
            return self.max_code
        
        a_neg = a < self.mid
        b_neg = b < self.mid
        
        if a_neg == b_neg:
            return self._add_same_sign(a, b, a_neg)
        else:
            return self._add_diff_sign(a, b, a_neg, b_neg)
    
    def _add_same_sign(self, a: int, b: int, is_neg: bool) -> int:
        """Add magnitudes when signs match."""
        a_mag = (self.mid - 1 - a) if is_neg else (a - self.mid)
        b_mag = (self.mid - 1 - b) if is_neg else (b - self.mid)
        
        if a_mag < b_mag:
            a_mag, b_mag = b_mag, a_mag
        
        d_levels = b_mag - a_mag
        d = d_levels * self.scale
        
        correction = self._softplus_cf(d)
        correction_levels = kdiv(correction, self.scale)
        
        c_mag = a_mag + correction_levels
        c_mag = max(0, min(c_mag, self.n_levels - 1))
        
        return (self.mid - 1 - c_mag) if is_neg else (self.mid + c_mag)
    
    def _add_diff_sign(self, a: int, b: int, a_neg: bool, b_neg: bool) -> int:
        """Subtract magnitudes when signs differ."""
        a_mag = (self.mid - 1 - a) if a_neg else (a - self.mid)
        b_mag = (self.mid - 1 - b) if b_neg else (b - self.mid)
        
        if a_mag == b_mag:
            return 0
        
        if a_mag > b_mag:
            big_mag, small_mag, result_neg = a_mag, b_mag, a_neg
        else:
            big_mag, small_mag, result_neg = b_mag, a_mag, b_neg
        
        d_levels = small_mag - big_mag
        d = d_levels * self.scale
        
        correction = self._softminus_cf(d)
        correction_levels = kdiv(correction, self.scale)
        
        c_mag = big_mag + correction_levels
        c_mag = max(0, min(c_mag, self.n_levels - 1))
        
        return (self.mid - 1 - c_mag) if result_neg else (self.mid + c_mag)
    
    def sub(self, a: int, b: int) -> int:
        """Subtract: add negation."""
        return self.add(a, self.neg(b))
    
    # ===== SIMILARITY =====
    
    def similarity(self, a: int, b: int) -> int:
        """CLZ(XOR) similarity metric."""
        return clz(a ^ b, self.n_bits)


# ===== PRESETS =====

def PBF4(v_min=0.1, v_max=10):
    return PBF(n_bits=4, v_min=v_min, v_max=v_max)

def PBF8(v_min=1e-4, v_max=1e4):
    return PBF(n_bits=8, v_min=v_min, v_max=v_max)

def PBF16(v_min=1e-8, v_max=1e8):
    return PBF(n_bits=16, v_min=v_min, v_max=v_max)


# ===== SNR TESTING =====

def snr_test(p, n_samples=1000, seed=42):
    """Compute SNR for PBF operations."""
    random.seed(seed)
    
    results = {}
    
    # Use the middle 80% of the range to avoid edge effects
    log_min = log10(p.v_min)
    log_max = log10(p.v_max)
    log_range = log_max - log_min
    safe_min = log_min + log_range * 0.1
    safe_max = log_max - log_range * 0.1
    
    # Roundtrip SNR
    signal_sq, noise_sq = 0.0, 0.0
    for _ in range(n_samples):
        v = 10 ** random.uniform(safe_min, safe_max)
        v *= random.choice([-1, 1])
        code = p.encode(v)
        back = p.decode(code)
        signal_sq += v * v
        noise_sq += (back - v) ** 2
    results['roundtrip'] = 10 * log10(signal_sq / max(noise_sq, 1e-100))
    
    # Addition SNR (same sign)
    signal_sq, noise_sq = 0.0, 0.0
    mid = (safe_min + safe_max) / 2
    add_range = (safe_max - safe_min) / 3
    for _ in range(n_samples):
        a = 10 ** random.uniform(mid - add_range, mid + add_range)
        b = 10 ** random.uniform(mid - add_range, mid + add_range)
        expected = a + b
        result = p.decode(p.add(p.encode(a), p.encode(b)))
        signal_sq += expected * expected
        noise_sq += (result - expected) ** 2
    results['add'] = 10 * log10(signal_sq / max(noise_sq, 1e-100))
    
    # Subtraction SNR (a > b, both positive)
    signal_sq, noise_sq = 0.0, 0.0
    for _ in range(n_samples):
        a = 10 ** random.uniform(mid, mid + add_range)
        b = 10 ** random.uniform(mid - add_range, mid)
        if a <= b:
            a, b = b * 1.1, a
        expected = a - b
        result = p.decode(p.sub(p.encode(a), p.encode(b)))
        signal_sq += expected * expected
        noise_sq += (result - expected) ** 2
    results['sub'] = 10 * log10(signal_sq / max(noise_sq, 1e-100))
    
    # Multiplication SNR - keep products in range
    signal_sq, noise_sq = 0.0, 0.0
    mul_range = min(add_range / 2, log_range / 4)
    for _ in range(n_samples):
        a = 10 ** random.uniform(mid - mul_range, mid + mul_range)
        b = 10 ** random.uniform(mid - mul_range, mid + mul_range)
        expected = a * b
        result = p.decode(p.mul(p.encode(a), p.encode(b)))
        signal_sq += expected * expected
        noise_sq += (result - expected) ** 2
    results['mul'] = 10 * log10(signal_sq / max(noise_sq, 1e-100))
    
    # Division SNR - keep quotients in range
    signal_sq, noise_sq = 0.0, 0.0
    for _ in range(n_samples):
        a = 10 ** random.uniform(mid - mul_range, mid + mul_range)
        b = 10 ** random.uniform(mid - mul_range, mid + mul_range)
        expected = a / b
        result = p.decode(p.div(p.encode(a), p.encode(b)))
        signal_sq += expected * expected
        noise_sq += (result - expected) ** 2
    results['div'] = 10 * log10(signal_sq / max(noise_sq, 1e-100))
    
    return results


def error_distribution(p, n_samples=1000, seed=42):
    """Analyze error distribution for add operation."""
    random.seed(seed)
    
    log_min = log10(p.v_min) + 1
    log_max = log10(p.v_max) - 1
    
    errors = []
    for _ in range(n_samples):
        a = 10 ** random.uniform(log_min, log_max)
        b = 10 ** random.uniform(log_min, log_max)
        expected = a + b
        result = p.decode(p.add(p.encode(a), p.encode(b)))
        if expected != 0:
            rel_err = (result - expected) / expected * 100
            errors.append(rel_err)
    
    errors.sort()
    n = len(errors)
    
    return {
        'min': errors[0],
        'p5': errors[int(n * 0.05)],
        'p25': errors[int(n * 0.25)],
        'median': errors[n // 2],
        'p75': errors[int(n * 0.75)],
        'p95': errors[int(n * 0.95)],
        'max': errors[-1],
        'mean': sum(errors) / n,
        'rms': (sum(e*e for e in errors) / n) ** 0.5
    }


# ===== DEMO =====

if __name__ == '__main__':
    print("=" * 65)
    print("PBF: Table-Free Continued Fraction Edition — SNR Analysis")
    print("=" * 65)
    
    print("\n" + "-" * 65)
    print("SNR by Format (higher = better, ~6 dB per bit)")
    print("-" * 65)
    
    formats = [
        ("PBF8 ", PBF8()),
        ("PBF16", PBF16()),
    ]
    
    print(f"\n{'Format':<8} {'Roundtrip':>12} {'Add':>12} {'Sub':>12} {'Mul':>12} {'Div':>12}")
    print("-" * 65)
    
    for name, fmt in formats:
        snr = snr_test(fmt, n_samples=2000)
        print(f"{name:<8}", end="")
        for op in ['roundtrip', 'add', 'sub', 'mul', 'div']:
            db = snr[op]
            bits = db / 6.02
            print(f" {db:5.1f}dB/{bits:.1f}b", end="")
        print()
    
    print("\n" + "-" * 65)
    print("Error Distribution for Addition (PBF8)")
    print("-" * 65)
    
    p8 = PBF8()
    dist = error_distribution(p8, n_samples=5000)
    
    print(f"\n  Min:    {dist['min']:+7.2f}%")
    print(f"  5th %:  {dist['p5']:+7.2f}%")
    print(f"  25th %: {dist['p25']:+7.2f}%")
    print(f"  Median: {dist['median']:+7.2f}%")
    print(f"  75th %: {dist['p75']:+7.2f}%")
    print(f"  95th %: {dist['p95']:+7.2f}%")
    print(f"  Max:    {dist['max']:+7.2f}%")
    print(f"  Mean:   {dist['mean']:+7.2f}%")
    print(f"  RMS:    {dist['rms']:7.2f}%")
    
    print("\n" + "-" * 65)
    print("CF Primitive Accuracy")
    print("-" * 65)
    
    import math
    
    p = PBF8()
    Q = p.Q
    
    tests = [
        ("ln(2)", p._ln_cf(2 * Q) / Q, math.log(2)),
        ("ln(e)", p._ln_cf(int(math.e * Q)) / Q, 1.0),
        ("ln(10)", p._ln_cf(10 * Q) / Q, math.log(10)),
        ("exp(1)", p._exp_cf(Q) / Q, math.e),
        ("exp(-1)", p._exp_cf(-Q) / Q, 1/math.e),
        ("exp(2)", p._exp_cf(2 * Q) / Q, math.e ** 2),
    ]
    
    print(f"\n  {'Function':<10} {'CF Result':>12} {'Exact':>12} {'Error':>10}")
    print("  " + "-" * 48)
    for name, result, exact in tests:
        err = abs(result - exact) / exact * 100
        print(f"  {name:<10} {result:12.6f} {exact:12.6f} {err:9.4f}%")
    
    print("\n" + "=" * 65)
    print("Pure CF. No tables. No floats in core ALU.")
    print("=" * 65)
