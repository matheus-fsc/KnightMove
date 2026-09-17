"""symmetry.py — Decomposição D₄ de tours fechados.

8 elementos de D₄: identidade, rot 90/180/270, 2 reflexões eixos, 2 reflexões diagonais.
canonical_form(tour) = menor representação lex entre as 8 imagens (cada uma normalizada
para iniciar em 0 e seguir direção com vizinho menor).
"""

import numpy as np


def _apply_transform(r, c, n, trans):
    if trans == 0:
        return r, c
    if trans == 1:
        return c, n - 1 - r
    if trans == 2:
        return n - 1 - r, n - 1 - c
    if trans == 3:
        return n - 1 - c, r
    if trans == 4:
        return r, n - 1 - c
    if trans == 5:
        return n - 1 - r, c
    if trans == 6:
        return c, r
    if trans == 7:
        return n - 1 - c, n - 1 - r
    raise ValueError(f"trans inválido: {trans}")


def d4_symmetries(n):
    """Retorna lista de 8 permutações np.int32 de tamanho V=n²."""
    V = n * n
    perms = []
    for trans in range(8):
        perm = np.zeros(V, dtype=np.int32)
        for v in range(V):
            r, c = divmod(v, n)
            r2, c2 = _apply_transform(r, c, n, trans)
            perm[v] = r2 * n + c2
        perms.append(perm)
    return perms


def _normalize_closed_tour(tour):
    """Canônica: inicia em vértice 0, segue direção com 2º vértice menor."""
    tour = np.asarray(tour, dtype=np.int32)
    V = len(tour)
    pos = np.where(tour == 0)[0]
    if len(pos) == 0:
        raise ValueError("tour não contém vértice 0")
    idx0 = int(pos[0])
    forward = np.concatenate([tour[idx0:], tour[:idx0]])
    rev = np.empty(V, dtype=np.int32)
    rev[0] = 0
    rev[1:] = forward[:0:-1]
    if int(forward[1]) <= int(rev[1]):
        return tuple(int(x) for x in forward)
    return tuple(int(x) for x in rev)


def canonical_form(tour, syms):
    """Menor representação lex entre as 8 imagens D₄ (cada uma normalizada)."""
    tour_arr = np.asarray(tour, dtype=np.int32)
    best = None
    for perm in syms:
        transformed = perm[tour_arr]
        norm = _normalize_closed_tour(transformed)
        if best is None or norm < best:
            best = norm
    return best


def is_canonical(tour_repr, syms):
    """True se tour_repr já está na sua forma canônica."""
    arr = np.asarray(tour_repr, dtype=np.int32)
    return tuple(int(x) for x in arr) == canonical_form(arr, syms)


def expand_d4(canonical_tours, n):
    """Expande tours canônicos por D₄, normalizando e deduplicando.

    Retorna set de tuplas (normalizadas). Cada órbita gera ≤ 8 tours
    (menos se o tour tem simetria interna)."""
    syms = d4_symmetries(n)
    out = set()
    for ct in canonical_tours:
        ct_arr = np.asarray(ct, dtype=np.int32)
        for perm in syms:
            transformed = perm[ct_arr]
            norm = _normalize_closed_tour(transformed)
            out.add(norm)
    return out


def verify_d4_decomposition(n, expected_total=None, expected_canonical=None,
                            K=None, verbose=True):
    """Enumera tours, computa formas canônicas, verifica expansão D₄.

    Retorna (canonical_set, unique_set)."""
    from .tours import knight_tours
    if K is None:
        K = 20000 if n == 6 else 200
    if verbose:
        print(f"Enumerando tours de {n}×{n} (K={K})...")
    tours = knight_tours(n, K, seed=0)
    unique = set(tuple(int(x) for x in t) for t in tours)
    if verbose:
        print(f"  {len(unique)} tours únicos encontrados")

    syms = d4_symmetries(n)
    canonical = set()
    for t in unique:
        canonical.add(canonical_form(np.asarray(t, dtype=np.int32), syms))
    if verbose:
        print(f"  {len(canonical)} formas canônicas")

    expanded = expand_d4(canonical, n)
    if verbose:
        print(f"  Expansão D₄: {len(expanded)} tours")

    if expected_total is not None:
        assert len(unique) == expected_total, (
            f"esperava {expected_total} tours, obteve {len(unique)}")
    if expected_canonical is not None:
        assert len(canonical) == expected_canonical, (
            f"esperava {expected_canonical} canônicos, obteve {len(canonical)}")
    assert expanded == unique, (
        f"expand_d4(canonical) tem {len(expanded)} elementos, "
        f"mas há {len(unique)} tours únicos — decomposição inconsistente")
    return canonical, unique


def knight_tours_canonical(n, K, seed=None):
    """Gera K tours canônicos (1 por órbita D₄), via enumeração + dedup."""
    from .tours import knight_tours
    syms = d4_symmetries(n)
    canonical_seen = set()
    canonical_list = []
    # Enumera com K_total alto e deduplica
    K_total = max(K * 12, 1000)
    tours = knight_tours(n, K_total, seed=seed)
    for t in tours:
        cf = canonical_form(t, syms)
        if cf not in canonical_seen:
            canonical_seen.add(cf)
            canonical_list.append(np.asarray(cf, dtype=np.int32))
            if len(canonical_list) >= K:
                break
    return canonical_list
