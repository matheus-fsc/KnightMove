#!/usr/bin/env python3
"""
Protótipo mínimo — estudo de viabilidade de "fusão XOR de loops (GF(2))"
aplicada à busca de sítios ativos em proteínas (estilo GASS).

O que ESTE script mede (e só isto):

  M1. O tamanho real do espaço de busca do casamento template->proteína,
      em função de k = número de resíduos do template, no grafo de
      correspondência Gamma (product graph). Serve para responder à
      pergunta "o limite prático ~5 é combinatório?".

  M2. Quantos candidatos passam no teste PAR-A-PAR (barato, local:
      |d(i,j) - d(a,b)| <= eps) e quantos sobrevivem ao teste GLOBAL
      (RMSD por superposição de Kabsch <= tau). Serve para identificar
      qual é, no caso proteico, o análogo do par "condição linear barata"
      vs "barreira cara" do passeio do cavalo.

  M3. Se o conjunto de soluções é fechado sob XOR (soma simétrica de
      conjuntos de arestas). Isto testa DIRETAMENTE a premissa da ideia:
      para que Z_1 = ker(partial_1) aja sobre as soluções, elas precisam
      ser subgrafos pares e precisam formar um coset de um subespaco.

O que este script NÃO mede: desempenho do GASS, qualidade biológica dos
casamentos, substituições conservativas, sítios inter-cadeia.

Uso:  ../venv/bin/python prototipo.py
Saída: resultados/medicoes.json  (+ resumo no stdout)
"""

from __future__ import annotations

import itertools
import json
import math
import pathlib
import time

import numpy as np

AQUI = pathlib.Path(__file__).resolve().parent
DADOS = AQUI / "dados"
SAIDA = AQUI / "resultados"

# Átomos de cadeia principal — o "last heavy atom" (LHA) usado pelo GASS
# é o último átomo pesado da CADEIA LATERAL.
BACKBONE = {"N", "CA", "C", "O", "OXT"}

EPS_PAR = 1.5      # tolerância par-a-par em Angstrom (critério do CatSId)
TAU_RMSD = 2.0     # limiar global de RMSD (Angstrom)
TAU_RMSD_FROUXO = 5.0  # limiar "GASS-like" (fitness <= 5 A)


# ----------------------------------------------------------------------
# Leitura de PDB (mínima, só ATOM da primeira MODEL)
# ----------------------------------------------------------------------
def ler_residuos(caminho: pathlib.Path, cadeia: str) -> dict:
    """Retorna {(cadeia, num_res, icode): (nome_res, coord_LHA)}."""
    atomos: dict = {}
    ordem: list = []
    with open(caminho) as fh:
        for linha in fh:
            if linha.startswith("ENDMDL"):
                break
            if not linha.startswith("ATOM"):
                continue
            alt = linha[16]
            if alt not in (" ", "A"):
                continue
            if linha[21] != cadeia:
                continue
            nome_atomo = linha[12:16].strip()
            elemento = linha[76:78].strip() or nome_atomo[0]
            if elemento == "H":
                continue
            nome_res = linha[17:20].strip()
            chave = (cadeia, linha[22:26].strip(), linha[26])
            if chave not in atomos:
                atomos[chave] = [nome_res, []]
                ordem.append(chave)
            xyz = (float(linha[30:38]), float(linha[38:46]), float(linha[46:54]))
            atomos[chave][1].append((nome_atomo, xyz))

    residuos = {}
    for chave in ordem:
        nome_res, lista = atomos[chave]
        lateral = [(a, c) for a, c in lista if a not in BACKBONE]
        if lateral:
            coord = lateral[-1][1]          # último átomo pesado da cadeia lateral
        else:
            ca = [c for a, c in lista if a == "CA"]
            if not ca:
                continue
            coord = ca[0]
        residuos[chave] = (nome_res, np.array(coord, dtype=float))
    return residuos


# ----------------------------------------------------------------------
# RMSD por superposição ótima (Kabsch)
# ----------------------------------------------------------------------
def rmsd_kabsch(P: np.ndarray, Q: np.ndarray) -> float:
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    diff = (R @ Pc.T).T - Qc
    return float(math.sqrt((diff ** 2).sum() / len(P)))


def matriz_dist(coords: np.ndarray) -> np.ndarray:
    d = coords[:, None, :] - coords[None, :, :]
    return np.sqrt((d ** 2).sum(axis=-1))


# ----------------------------------------------------------------------
# Grafo de correspondência Gamma
# ----------------------------------------------------------------------
def componentes(nV: int, arestas) -> int:
    pai = list(range(nV))

    def acha(x):
        while pai[x] != x:
            pai[x] = pai[pai[x]]
            x = pai[x]
        return x

    for u, v in arestas:
        ru, rv = acha(u), acha(v)
        if ru != rv:
            pai[ru] = rv
    return len({acha(x) for x in range(nV)})


def construir_gamma(template: list, alvo: dict, eps: float = EPS_PAR,
                    modo: str = "tipo_exato"):
    """
    template: lista de (nome_res, coord)  — k resíduos
    alvo:     {chave: (nome_res, coord)}  — proteína varrida

    modo="tipo_exato":  nós de Gamma são pares (i, a) com o alvo do MESMO tipo
                        de aminoácido (substituição conservativa não modelada).
    modo="geometrico":  nós são todos os pares (i, a) — casamento puramente
                        geométrico. É o regime em que o custo combinatório
                        aparece; serve para medir a explosão em k.
    Arestas: {(i,a),(j,b)} com i != j, a != b e
             | d_template(i,j) - d_alvo(a,b) | <= eps.
    """
    k = len(template)
    chaves = list(alvo.keys())
    idx_alvo = {c: n for n, c in enumerate(chaves)}
    coords_alvo = np.array([alvo[c][1] for c in chaves])
    D_alvo = matriz_dist(coords_alvo)

    coords_tpl = np.array([c for _, c in template])
    D_tpl = matriz_dist(coords_tpl)

    nos = []
    por_indice = [[] for _ in range(k)]
    for i, (nome_i, _) in enumerate(template):
        for c in chaves:
            if modo == "geometrico" or alvo[c][0] == nome_i:
                n = len(nos)
                nos.append((i, c))
                por_indice[i].append(n)

    arestas = set()
    for u in range(len(nos)):
        i, a = nos[u]
        for v in range(u + 1, len(nos)):
            j, b = nos[v]
            if i == j or a == b:
                continue
            if abs(D_tpl[i, j] - D_alvo[idx_alvo[a], idx_alvo[b]]) <= eps:
                arestas.add((u, v))

    return {
        "k": k,
        "nos": nos,
        "por_indice": por_indice,
        "arestas": arestas,
        "coords_alvo": coords_alvo,
        "idx_alvo": idx_alvo,
        "coords_tpl": coords_tpl,
    }


def enumerar_casamentos(G, tau=TAU_RMSD, tau_frouxo=TAU_RMSD_FROUXO, teto=2_000_000):
    """
    Backtracking sobre os índices do template (um nó de Gamma por índice),
    exigindo aresta em Gamma entre todo par já escolhido -> k-clique
    transversal. Depois aplica o teste GLOBAL de RMSD.

    Retorna estatísticas e a lista de casamentos válidos (como tuplas de nós).
    """
    k = G["k"]
    arestas = G["arestas"]
    nos = G["nos"]
    coords_alvo = G["coords_alvo"]
    idx_alvo = G["idx_alvo"]
    coords_tpl = G["coords_tpl"]

    def adj(u, v):
        return (u, v) in arestas or (v, u) in arestas

    est = {"nos_backtracking": 0, "cliques_par_a_par": 0,
           "validos_tau": 0, "validos_tau_frouxo": 0, "estourou": False}
    validos = []
    usados_alvo = set()

    def rec(i, parcial):
        if est["estourou"]:
            return
        if est["nos_backtracking"] > teto:
            est["estourou"] = True
            return
        if i == k:
            est["cliques_par_a_par"] += 1
            P = coords_tpl
            Q = np.array([coords_alvo[idx_alvo[nos[u][1]]] for u in parcial])
            r = rmsd_kabsch(P, Q)
            if r <= tau_frouxo:
                est["validos_tau_frouxo"] += 1
            if r <= tau:
                est["validos_tau"] += 1
                validos.append((tuple(parcial), r))
            return
        for u in G["por_indice"][i]:
            chave = nos[u][1]
            if chave in usados_alvo:
                continue
            if all(adj(u, w) for w in parcial):
                est["nos_backtracking"] += 1
                usados_alvo.add(chave)
                parcial.append(u)
                rec(i + 1, parcial)
                parcial.pop()
                usados_alvo.discard(chave)

    rec(0, [])
    return est, validos


# ----------------------------------------------------------------------
# M3: fechamento sob XOR
# ----------------------------------------------------------------------
def arestas_do_casamento(clique):
    """Conjunto de arestas (como frozenset de pares ordenados) da k-clique."""
    return frozenset(
        (min(u, v), max(u, v)) for u, v in itertools.combinations(sorted(clique), 2)
    )


def testar_fechamento_xor(validos, max_pares=5000):
    """
    Testa se o conjunto de soluções é fechado sob XOR de conjuntos de arestas.
    Para cada par (A,B) de soluções distintas, calcula A ^ B e verifica:
      - é o conjunto de arestas de alguma k-clique transversal? (validade)
      - todos os vértices têm grau PAR? (pertence a Z_1)
    """
    sol_edges = [arestas_do_casamento(c) for c, _ in validos]
    conjunto = set(sol_edges)
    n_par = 0
    n_valido = 0
    n_grau_par = 0
    for A, B in itertools.islice(itertools.combinations(sol_edges, 2), max_pares):
        X = A ^ B
        n_par += 1
        if X in conjunto:
            n_valido += 1
        graus = {}
        for u, v in X:
            graus[u] = graus.get(u, 0) + 1
            graus[v] = graus.get(v, 0) + 1
        if all(g % 2 == 0 for g in graus.values()):
            n_grau_par += 1
    return {"pares_testados": n_par,
            "xor_e_solucao": n_valido,
            "xor_tem_grau_par": n_grau_par}


def paridade_da_clique(k):
    """Uma k-clique tem todos os graus k-1: pertence a Z_1 sse k é ímpar."""
    return (k - 1) % 2 == 0


# ----------------------------------------------------------------------
# Montagem do experimento
# ----------------------------------------------------------------------
# Tríade catalítica das serino-proteases (numeração do quimotripsinogênio).
TRIADE = [("1PPF", "E", "57"), ("1PPF", "E", "102"), ("1PPF", "E", "195")]


def montar_template(res_fonte, k):
    """
    k=3 -> tríade catalítica His57/Asp102/Ser195 de 1PPF cadeia E.
    k>3 -> acrescenta os resíduos mais próximos do centroide da tríade
           (qualquer tipo), imitando um template de sítio ativo maior.
    """
    base = []
    for _, cad, num in TRIADE:
        chave = (cad, num, " ")
        nome, coord = res_fonte[chave]
        base.append((nome, coord, chave))
    if k <= 3:
        base = base[:k]
    else:
        centro = np.mean([c for _, c, _ in base], axis=0)
        usados = {ch for _, _, ch in base}
        cand = [(np.linalg.norm(c - centro), ch, n, c)
                for ch, (n, c) in res_fonte.items() if ch not in usados]
        cand.sort(key=lambda t: t[0])
        for _, ch, n, c in cand[: k - 3]:
            base.append((n, c, ch))
    return [(n, c) for n, c, _ in base], [ch for _, _, ch in base]


def main():
    SAIDA.mkdir(exist_ok=True)
    fonte = ler_residuos(DADOS / "1PPF.pdb", "E")     # elastase leucocitária
    alvos = {
        "1ACB_E (quimotripsina, homóloga)": ler_residuos(DADOS / "1ACB.pdb", "E"),
        "1LYZ_A (lisozima, controle negativo)": ler_residuos(DADOS / "1LYZ.pdb", "A"),
    }

    relatorio = {"eps_par": EPS_PAR, "tau_rmsd": TAU_RMSD,
                 "tau_rmsd_frouxo": TAU_RMSD_FROUXO, "execucoes": []}

    for modo in ("tipo_exato", "geometrico"):
        for nome_alvo, alvo in alvos.items():
            print(f"\n=== modo={modo} | alvo: {nome_alvo} ({len(alvo)} resíduos) ===")
            for k in (3, 4, 5, 6, 7):
                template, chaves_tpl = montar_template(fonte, k)
                t0 = time.time()
                G = construir_gamma(template, alvo, modo=modo)
                est, validos = enumerar_casamentos(G)
                dt = time.time() - t0

                nV = len(G["nos"])
                nE = len(G["arestas"])
                comp = componentes(nV, G["arestas"])
                beta1 = nE - nV + comp
                espaco_bruto = 1
                for lst in G["por_indice"]:
                    espaco_bruto *= max(len(lst), 0)

                xor = testar_fechamento_xor(validos) if len(validos) >= 2 else None

                linha = {
                    "modo": modo,
                    "alvo": nome_alvo,
                    "k": k,
                    "template": [f"{n}{ch[1]}" for (n, _), ch in zip(template, chaves_tpl)],
                    "gamma_V": nV,
                    "gamma_E": nE,
                    "gamma_componentes": comp,
                    "gamma_beta1": beta1,
                    "espaco_bruto_produto": espaco_bruto,
                    "nos_backtracking": est["nos_backtracking"],
                    "cliques_par_a_par": est["cliques_par_a_par"],
                    "validos_rmsd_2A": est["validos_tau"],
                    "validos_rmsd_5A": est["validos_tau_frouxo"],
                    "estourou_teto": est["estourou"],
                    "clique_esta_em_Z1": paridade_da_clique(k),
                    "segundos": round(dt, 3),
                    "xor": xor,
                }
                relatorio["execucoes"].append(linha)
                print(
                    f"  k={k}  |V(Γ)|={nV:5d} |E(Γ)|={nE:7d} β₁={beta1:7d}  "
                    f"produto={espaco_bruto:>15,}  cliques={est['cliques_par_a_par']:7d}  "
                    f"RMSD≤5Å={est['validos_tau_frouxo']:6d}  RMSD≤2Å={est['validos_tau']:5d}  "
                    f"k-clique∈Z₁={paridade_da_clique(k)}  ({dt:.2f}s)"
                )
                if xor:
                    print(f"        XOR: {xor['xor_e_solucao']}/{xor['pares_testados']} pares "
                          f"dão solução; {xor['xor_tem_grau_par']} têm grau par")

    with open(SAIDA / "medicoes.json", "w") as fh:
        json.dump(relatorio, fh, indent=2, ensure_ascii=False)
    print(f"\nEscrito: {SAIDA / 'medicoes.json'}")


if __name__ == "__main__":
    main()
