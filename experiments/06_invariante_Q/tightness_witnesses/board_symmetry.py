#!/usr/bin/env python3
"""
board_symmetry.py
=================
Simetrias do tabuleiro aplicadas a flip-realizabilidade.

O grupo de automorfismos do grafo do cavalo num retangulo R x C e' o grupo de
Klein (ordem 4) quando R != C, e D4 (ordem 8) quando R == C. Como sao
automorfismos, **flip-realizabilidade e' constante em cada orbita** -- e o
certificado se transporta: se (tau, tau XOR C) certifica C, entao
(g.tau, g.tau XOR g.C) certifica g.C.

Isso rende duas economias, deliberadamente separadas:

  (1) ANTES da busca -- `orbit_representatives`: quocientar os hexagonos do
      bulk pelo grupo e testar so' um representante por orbita. Fator ~4 em
      retangulo, ~8 em quadrado.

  (2) DEPOIS da bola -- `close_certificates`: fechar o conjunto certificado
      sob o grupo antes de mandar o residuo para a busca exaustiva. A bola e'
      assimetrica (ela cobre o que calha de cobrir), entao o residuo dela
      quase sempre contem hexagonos cujo espelho ja' esta' certificado.

O custo de nao fazer (2) foi medido: no 7x8 o residuo de 1 hexagono consumiu
**4.781 s de busca exaustiva e voltou INCONCLUSIVO**, enquanto o transporte
por simetria o resolveu em milissegundos com certificado verificado.

E (1) teria pago antes: os 12 excepcionais do 6x6 sao 2 orbitas D4 (4+8) e os
4 do 6x7 sao 1 orbita de Klein.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

from flip_graph import cycle_mask

Perm = Callable[[int, int], Tuple[int, int]]


def board_group(R: int, C: int) -> Dict[str, Perm]:
    """Automorfismos do tabuleiro: Klein se R != C, D4 se R == C."""
    g: Dict[str, Perm] = {
        "id": lambda r, c: (r, c),
        "flip-h": lambda r, c: (r, C - 1 - c),
        "flip-v": lambda r, c: (R - 1 - r, c),
        "rot180": lambda r, c: (R - 1 - r, C - 1 - c),
    }
    if R == C:
        n = R
        g.update({
            "rot90": lambda r, c: (c, n - 1 - r),
            "rot270": lambda r, c: (n - 1 - c, r),
            "diag": lambda r, c: (c, r),
            "anti": lambda r, c: (n - 1 - c, n - 1 - r),
        })
    return g


def apply_mask(mask: int, g: Perm, C: int, edges, eidx) -> int:
    """Imagem de um conjunto de arestas (mascara) sob `g`."""
    out = 0
    x = mask
    while x:
        b = x & -x
        u, w = edges[b.bit_length() - 1]
        a = g(*divmod(u, C))
        d = g(*divmod(w, C))
        a = a[0] * C + a[1]
        d = d[0] * C + d[1]
        out ^= 1 << eidx[(min(a, d), max(a, d))]
        x ^= b
    return out


def orbit_representatives(R: int, C: int, hexes: Sequence, edges, eidx):
    """(reps, orbit_of) -- um representante por orbita, e o mapa
    mascara -> mascara do representante."""
    G = board_group(R, C)
    masks = {cycle_mask(c, eidx): c for c in hexes}
    orbit_of: Dict[int, int] = {}
    reps = []
    for cm, cyc in masks.items():
        if cm in orbit_of:
            continue
        reps.append(cyc)
        for g in G.values():
            im = apply_mask(cm, g, C, edges, eidx)
            if im in masks:
                orbit_of[im] = cm
    return reps, orbit_of


def close_certificates(R: int, C: int, realized: Dict[int, Tuple[int, int]],
                       hex_masks, edges, eidx, ham) -> Dict[int, Tuple[int, int]]:
    """Fecha `realized` (mascara do hexagono -> par (tau, tau XOR C)) sob o
    grupo. Cada certificado transportado e' VERIFICADO, nao assumido:
    os dois tours transportados tem de ser hamiltonianos e o XOR tem de dar
    exatamente o hexagono imagem."""
    G = board_group(R, C)
    out = dict(realized)
    changed = True
    while changed:
        changed = False
        for cm, (t1, t2) in list(out.items()):
            for name, g in G.items():
                if name == "id":
                    continue
                icm = apply_mask(cm, g, C, edges, eidx)
                if icm in out or icm not in hex_masks:
                    continue
                a = apply_mask(t1, g, C, edges, eidx)
                b = apply_mask(t2, g, C, edges, eidx)
                if not ham(a) or not ham(b) or (a ^ b) != icm:
                    continue      # nunca assumir: so' entra se verificado
                out[icm] = (a, b)
                changed = True
    return out
