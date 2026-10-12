"""Minimal, dependency-free Jordan-Wigner mapping used ONLY for the audit cross-checks.

Conventions (stated explicitly because Pauli strings depend on them; energies do not):
  * Spin-orbital index  k = 2*p + s, with s = 0 (alpha), 1 (beta)  -> "interleaved" ordering.
  * Qubit k <-> spin-orbital k. Occupied = |1>, empty = |0>.
  * a_k  = Z_0 ... Z_{k-1} (X_k + i Y_k)/2,   a_k^dag = Z_0 ... Z_{k-1} (X_k - i Y_k)/2.
  * A Pauli string is stored as a pair of bitmasks (x, z), meaning P = i^{popcount(x&z)} X^x Z^z,
    so that (1,1) on one qubit is Y. Every such P is Hermitian.
  * Dense matrices use qubit 0 as the LEFT-most Kronecker factor.
Hamiltonian (chemists' notation for the spatial integrals (pq|rs)):
  H = E_nuc + sum_{pq,s} h_pq a+_{ps} a_{qs} + 1/2 sum_{pqrs,s,t} (pq|rs) a+_{ps} a+_{rt} a_{st} a_{qs}
"""
from collections import defaultdict
import numpy as np


def _mul(a, b):
    """Product of Pauli strings a=(x1,z1), b=(x2,z2). Returns (phase_exponent_mod4, (x, z))."""
    x1, z1 = a
    x2, z2 = b
    x, z = x1 ^ x2, z1 ^ z2
    e = ((x1 & z1).bit_count() + (x2 & z2).bit_count() - (x & z).bit_count()
         + 2 * (z1 & x2).bit_count()) % 4
    return e, (x, z)


_IPOW = (1, 1j, -1, -1j)


def ladder(k, dagger):
    zs = (1 << k) - 1
    xk = 1 << k
    s = -0.5j if dagger else 0.5j
    return [(0.5, (xk, zs)), (s, (xk, zs | xk))]


def product(ops):
    """Multiply a list of operators, each a list of (coef, (x, z)). Returns dict."""
    terms = {(0, 0): 1.0 + 0j}
    for op in ops:
        new = defaultdict(complex)
        for p, c in terms.items():
            for c2, q in op:
                e, r = _mul(p, q)
                new[r] += c * c2 * _IPOW[e]
        terms = new
    return terms


def molecular_qubit_hamiltonian(ecore, h1, eri, tol=1e-12):
    """h1[p,q] (spatial MO basis), eri[p,q,r,s] = (pq|rs) chemists' notation (full 4-index array).
    Returns dict {(x, z): complex coefficient} for the JW-mapped Hamiltonian on 2*norb qubits."""
    norb = h1.shape[0]
    H = defaultdict(complex)
    H[(0, 0)] += ecore
    for p in range(norb):
        for q in range(norb):
            if abs(h1[p, q]) < tol:
                continue
            for s in (0, 1):
                for key, c in product([ladder(2 * p + s, True), ladder(2 * q + s, False)]).items():
                    H[key] += h1[p, q] * c
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for s_ in range(norb):
                    v = eri[p, q, r, s_]
                    if abs(v) < tol:
                        continue
                    for a in (0, 1):
                        for b in (0, 1):
                            P, Q, R, S = 2 * p + a, 2 * q + a, 2 * r + b, 2 * s_ + b
                            if P == R or S == Q:
                                continue  # a+_P a+_P = 0 and a_S a_S = 0
                            ops = [ladder(P, True), ladder(R, True), ladder(S, False), ladder(Q, False)]
                            for key, c in product(ops).items():
                                H[key] += 0.5 * v * c
    return {k: c for k, c in H.items() if abs(c) > tol}


def pauli_label(key, nq):
    x, z = key
    out = []
    for k in range(nq):
        bx, bz = (x >> k) & 1, (z >> k) & 1
        if bx or bz:
            out.append({(1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(bx, bz)] + str(k))
    return " ".join(out) if out else "I"


_P1 = {
    (0, 0): np.eye(2, dtype=complex),
    (1, 0): np.array([[0, 1], [1, 0]], dtype=complex),
    (0, 1): np.array([[1, 0], [0, -1]], dtype=complex),
    (1, 1): np.array([[0, -1j], [1j, 0]], dtype=complex),
}


def pauli_matrix(key, nq):
    x, z = key
    m = np.ones((1, 1), dtype=complex)
    for k in range(nq):
        m = np.kron(m, _P1[((x >> k) & 1, (z >> k) & 1)])
    return m


def dense_from_paulis(H, nq):
    M = np.zeros((2 ** nq, 2 ** nq), dtype=complex)
    for key, c in H.items():
        M += c * pauli_matrix(key, nq)
    return M


def dense_fermion_ops(nq):
    """Independent construction: explicit JW annihilation matrices built from 2x2 blocks."""
    Z = np.diag([1.0, -1.0]).astype(complex)
    I2 = np.eye(2, dtype=complex)
    lower = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1| : annihilates an occupied mode
    ops = []
    for k in range(nq):
        m = np.ones((1, 1), dtype=complex)
        for j in range(nq):
            m = np.kron(m, Z if j < k else (lower if j == k else I2))
        ops.append(m)
    return ops


def dense_hamiltonian_direct(ecore, h1, eri):
    """Build H directly from matrix products of annihilation/creation matrices (no Pauli algebra)."""
    norb = h1.shape[0]
    nq = 2 * norb
    a = dense_fermion_ops(nq)
    ad = [m.conj().T for m in a]
    H = ecore * np.eye(2 ** nq, dtype=complex)
    for p in range(norb):
        for q in range(norb):
            for s in (0, 1):
                H += h1[p, q] * ad[2 * p + s] @ a[2 * q + s]
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for s_ in range(norb):
                    for x in (0, 1):
                        for y in (0, 1):
                            H += 0.5 * eri[p, q, r, s_] * (
                                ad[2 * p + x] @ ad[2 * r + y] @ a[2 * s_ + y] @ a[2 * q + x])
    return H


def number_and_sz(nq):
    a = dense_fermion_ops(nq)
    n_ops = [m.conj().T @ m for m in a]
    N = sum(n_ops)
    Sz = 0.5 * (sum(n_ops[0::2]) - sum(n_ops[1::2]))
    return N, Sz
