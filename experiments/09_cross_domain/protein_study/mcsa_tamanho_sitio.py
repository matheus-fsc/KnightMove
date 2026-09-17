#!/usr/bin/env python3
"""
Mede a distribuição do NÚMERO DE RESÍDUOS por sítio catalítico curado no
M-CSA. Serve só para uma coisa: decidir a ambiguidade "limite ~5" —
5 proteínas ou 5 resíduos?

Fonte: https://www.ebi.ac.uk/thornton-srv/m-csa/media/flat_files/curated_data.csv
(baixado em dados/mcsa_curated_data.csv)

Uso: ../venv/bin/python mcsa_tamanho_sitio.py
"""
import csv
import json
import pathlib
import statistics
from collections import Counter, defaultdict

AQUI = pathlib.Path(__file__).resolve().parent
CSV = AQUI / "dados" / "mcsa_curated_data.csv"
SAIDA = AQUI / "resultados" / "mcsa_tamanho_sitio.json"


def main():
    todos = defaultdict(set)      # entrada -> {(cadeia, resid)}
    nao_espectador = defaultdict(set)
    with open(CSV, newline="") as fh:
        leitor = csv.DictReader(fh)
        for r in leitor:
            if r["residue/reactant/product/cofactor"] != "residue":
                continue
            eid = r["M-CSA ID"]
            chave = (r["chain/kegg compound"], r["resid/chebi id"])
            todos[eid].add(chave)
            if r["role type"].strip().lower() != "spectator":
                nao_espectador[eid].add(chave)

    def resumo(d, rotulo):
        tam = sorted(len(v) for v in d.values())
        c = Counter(tam)
        info = {
            "rotulo": rotulo,
            "n_entradas": len(tam),
            "media": round(statistics.mean(tam), 2),
            "mediana": statistics.median(tam),
            "p90": tam[int(0.90 * len(tam)) - 1],
            "p95": tam[int(0.95 * len(tam)) - 1],
            "max": tam[-1],
            "frac_ate_5": round(sum(1 for t in tam if t <= 5) / len(tam), 4),
            "frac_ate_6": round(sum(1 for t in tam if t <= 6) / len(tam), 4),
            "histograma": {str(k): c[k] for k in sorted(c)},
        }
        print(f"\n--- {rotulo} ---")
        print(f"entradas={info['n_entradas']}  média={info['media']}  "
              f"mediana={info['mediana']}  p90={info['p90']}  p95={info['p95']}  "
              f"máx={info['max']}")
        print(f"fração com <= 5 resíduos: {info['frac_ate_5']:.1%}   "
              f"<= 6: {info['frac_ate_6']:.1%}")
        print("histograma (tamanho: nº de sítios):",
              " ".join(f"{k}:{v}" for k, v in info["histograma"].items()))
        return info

    r = {
        "todos_residuos_anotados": resumo(todos, "todos os resíduos anotados"),
        "sem_espectadores": resumo(nao_espectador, "excluindo role type=spectator"),
    }
    SAIDA.parent.mkdir(exist_ok=True)
    with open(SAIDA, "w") as fh:
        json.dump(r, fh, indent=2, ensure_ascii=False)
    print(f"\nEscrito: {SAIDA}")


if __name__ == "__main__":
    main()
