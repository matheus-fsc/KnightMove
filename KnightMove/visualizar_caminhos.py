"""Visualizador gráfico dos caminhos do movimento do cavalo.

Lê o arquivo 'caminho.txt' (mesma pasta) com o formato:

	Caminho 1:
	0,0
	1,2
	...

Cada bloco começa com a linha 'Caminho N:' seguida das coordenadas (x,y) uma por linha
até uma linha em branco ou até começar o próximo bloco.

Funcionalidades:
 - Seleção de caminho via Combobox ou teclado (Setas ↑↓ mudam caminho)
 - Animação passo a passo (Play/Pause) com controle de velocidade
 - Botões: Primeiro, Anterior, Play/Pause, Próximo, Último, Redesenhar
 - Atalhos teclado:
	   Espaço -> Play/Pause
	   ← / →  -> Passo anterior / próximo
	   Home / End -> Primeiro / Último
 - Exibição: Tabuleiro colorido xadrez, numeração da ordem de visita, cavalo atual
 - Linhas conectando os movimentos (opcional desligar via checkbox)

Possíveis melhorias futuras (ver final do arquivo): exportar GIF, salvar PNG, filtro de tours completos, etc.
"""

from __future__ import annotations

import re
import sys
import math
import time
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

try:
	import tkinter as tk
	from tkinter import ttk, messagebox
except ImportError:  # pragma: no cover - ambiente sem Tk
	print("Tkinter não disponível neste ambiente.")
	sys.exit(1)


CAMINHO_ARQ = Path(__file__).with_name("caminho.txt")


Coord = Tuple[int, int]


@dataclass
class KnightPath:
	indice: int
	coords: List[Coord]

	def is_closed_tour(self) -> bool:
		"""Verifica se é um tour fechado (retorna ao início em um salto de cavalo)."""
		if len(self.coords) < 2:
			return False
		x0, y0 = self.coords[0]
		xn, yn = self.coords[-1]
		return (abs(x0 - xn), abs(y0 - yn)) in {(1, 2), (2, 1)}

	def is_full_cover(self, board_size: int) -> bool:
		return len(self.coords) == board_size * board_size and len(set(self.coords)) == len(self.coords)


def parse_paths(path: Path) -> List[KnightPath]:
	if not path.exists():
		raise FileNotFoundError(f"Arquivo não encontrado: {path}")

	content = path.read_text(encoding="utf-8", errors="ignore").splitlines()
	paths: List[KnightPath] = []
	current: List[Coord] = []
	current_index: Optional[int] = None

	header_re = re.compile(r"^Caminho\s+(\d+):\s*$", re.IGNORECASE)
	coord_re = re.compile(r"^\s*(\d+)\s*,\s*(\d+)\s*$")

	def flush():
		nonlocal current, current_index
		if current_index is not None and current:
			paths.append(KnightPath(current_index, current))
		current = []
		current_index = None

	for line in content:
		if not line.strip():
			# separador
			continue
		mh = header_re.match(line)
		if mh:
			# novo caminho
			if current_index is not None:
				flush()
			current_index = int(mh.group(1))
			continue
		mc = coord_re.match(line)
		if mc and current_index is not None:
			x, y = int(mc.group(1)), int(mc.group(2))
			current.append((x, y))
		else:
			# linha inesperada -> flush anterior (robustez) e tenta continuar
			flush()
	# flush final
	if current_index is not None:
		flush()

	# Ordena por índice (caso venham fora de ordem)
	paths.sort(key=lambda p: p.indice)
	return paths


class KnightPathsViewer:

	def selecionar_menos_colisoes_hipotenusa(self):
		"""Seleciona o caminho com menos colisões de hipotenusas válidas (não consecutivas)."""
		def hipotenusa(tri):
			from math import dist
			a, b, c = tri
			dists = [(a, b, dist(a, b)), (a, c, dist(a, c)), (b, c, dist(b, c))]
			dists.sort(key=lambda x: -x[2])
			return (dists[0][0], dists[0][1])
		def segments_cross(seg1, seg2):
			def ccw(A, B, C):
				return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
			A, B = seg1
			C, D = seg2
			return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))
		min_idx = 0
		min_col = None
		for idx, path in enumerate(self.paths):
			triangles = self._optimized_triangles(path)
			hipots = [hipotenusa(tri) for tri, cells, alt_tri, alt_cells, opt_idx in triangles]
			colisoes = 0
			for i in range(len(hipots)):
				for j in range(i+1, len(hipots)):
					if i+1 == j or i == j+1:
						continue
					v1 = set(hipots[i])
					v2 = set(hipots[j])
					if v1 & v2:
						continue
					if segments_cross(hipots[i], hipots[j]):
						colisoes += 1
			if min_col is None or colisoes < min_col:
				min_col = colisoes
				min_idx = idx
		self._load_path(min_idx)
		self.combo.current(min_idx)

	def _optimized_triangles(self, path: KnightPath):
		"""Para cada bloco de 7 saltos, faz busca exaustiva para minimizar colisões locais naquele bloco.
		Retorna lista de tuplas: (tri_escolhido, cells_escolhido, tri_alternativo, cells_alternativo, opt_idx)
		"""
		from itertools import product
		N = 7  # tamanho do bloco
		moves = []
		for i in range(1, len(path.coords)):
			x1, y1 = path.coords[i - 1]
			x2, y2 = path.coords[i]
			dx = x2 - x1
			dy = y2 - y1
			if (abs(dx), abs(dy)) not in {(1, 2), (2, 1)}:
				moves.append(None)
				continue
			# Duas opções de canto intermediário
			if abs(dx) == 2:
				mid1 = (x1 + (2 if dx > 0 else -2), y1)
				mid2 = (x2, y1)
			else:
				mid1 = (x1, y1 + (2 if dy > 0 else -2))
				mid2 = (x1, y2)
			options = []
			for mid in [mid1, mid2]:
				triangle = [(x1, y1), mid, (x2, y2)]
				cells = set()
				for x in range(self.board_size):
					for y in range(self.board_size):
						if self._point_in_triangle((x + 0.5, y + 0.5), triangle):
							cells.add((x, y))
				options.append((triangle, cells))
			moves.append(options)

		seq = []
		i = 0
		while i < len(moves):
			bloco = moves[i:i+N]
			n_moves = len(bloco)
			best_seq = None
			best_col = None
			# Busca exaustiva local para o bloco
			for choices in product([0, 1], repeat=n_moves):
				bloco_seq = []
				covered = set()
				colisoes = 0
				for j, opt in enumerate(choices):
					tri, cells = bloco[j][opt]
					alt_tri, alt_cells = bloco[j][1-opt]
					bloco_seq.append((tri, cells, alt_tri, alt_cells, opt))
					colisoes += len(cells & covered)
					covered.update(cells)
				if best_col is None or colisoes < best_col:
					best_col = colisoes
					best_seq = list(bloco_seq)
			seq.extend(best_seq)
			i += N
		return seq
	def selecionar_mais_colisoes(self):
		"""Seleciona o caminho com mais colisões de triângulos."""
		max_idx = 0
		max_col = -1
		for idx, path in enumerate(self.paths):
			triangles_cells = self._all_triangles_cells(path)
			colisoes = 0
			for i in range(len(triangles_cells)):
				for j in range(i+1, len(triangles_cells)):
					if triangles_cells[i] & triangles_cells[j]:
						colisoes += 1
			if colisoes > max_col:
				max_col = colisoes
				max_idx = idx
		self._load_path(max_idx)
		self.combo.current(max_idx)

	def selecionar_menos_colisoes(self):
		"""Seleciona o caminho com menos colisões de triângulos."""
		min_idx = 0
		min_col = None
		for idx, path in enumerate(self.paths):
			triangles_cells = self._all_triangles_cells(path)
			colisoes = 0
			for i in range(len(triangles_cells)):
				for j in range(i+1, len(triangles_cells)):
					if triangles_cells[i] & triangles_cells[j]:
						colisoes += 1
			if min_col is None or colisoes < min_col:
				min_col = colisoes
				min_idx = idx
		self._load_path(min_idx)
		self.combo.current(min_idx)

	def _triangles_covered_cells(self, path: KnightPath) -> set:
		"""Retorna o conjunto de casas cobertas por todos os triângulos do caminho."""
		covered = set()
		for i in range(1, len(path.coords)):
			x1, y1 = path.coords[i - 1]
			x2, y2 = path.coords[i]
			dx = x2 - x1
			dy = y2 - y1
			if (abs(dx), abs(dy)) not in {(1, 2), (2, 1)}:
				continue
			if abs(dx) == 2:
				mid = (x1 + (2 if dx > 0 else -2), y1)
			elif abs(dy) == 2:
				mid = (x1, y1 + (2 if dy > 0 else -2))
			else:
				continue
			# Triângulo: (x1,y1), mid, (x2,y2)
			# Cobertura: todas as casas dentro do triângulo (aproximação: casas cujos centros estão dentro do triângulo)
			triangle = [(x1, y1), mid, (x2, y2)]
			for x in range(self.board_size):
				for y in range(self.board_size):
					if self._point_in_triangle((x + 0.5, y + 0.5), triangle):
						covered.add((x, y))
			return covered

	def _triangles_covered_cells_countmap(self, path: KnightPath) -> dict:
		"""Retorna um dict {(x,y): n} com o número de vezes que cada casa é coberta por triângulo."""
		countmap = {}
		for tri_cells in self._all_triangles_cells(path):
			for cell in tri_cells:
				countmap[cell] = countmap.get(cell, 0) + 1
		return countmap

	def _all_triangles_cells(self, path: KnightPath):
		"""Retorna uma lista de conjuntos, cada um com as casas cobertas por um triângulo do caminho, usando a escolha otimizada."""
		return [cells for tri, cells, alt_tri, alt_cells, opt_idx in self._optimized_triangles(path)]

	@staticmethod
	def _point_in_triangle(p, tri):
		"""Verifica se ponto p está dentro do triângulo tri (barycentric)."""
		(x, y) = p
		(x1, y1), (x2, y2), (x3, y3) = tri
		den = ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3))
		if den == 0:
			return False
		a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / den
		b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / den
		c = 1 - a - b
		return 0 <= a <= 1 and 0 <= b <= 1 and 0 <= c <= 1

	def _miolo_cells(self):
		"""Retorna o conjunto de casas do miolo (excluindo bordas)."""
		size = self.board_size
		return set((x, y) for x in range(1, size - 1) for y in range(1, size - 1))

	def __init__(self, root: tk.Tk, paths: List[KnightPath]):
		self.play_btn = None
		self.root = root
		self.paths = paths
		if not paths:
			messagebox.showerror("Erro", "Nenhum caminho válido encontrado em caminho.txt")
			root.destroy()
			return

		self.board_size = self._infer_board_size()
		self.square_size = 90 if self.board_size <= 6 else 60
		self.margin = 40
		self.font_cache: Dict[int, tk.Font] = {}

		self.current_path_index = 0  # índice dentro da lista paths
		self.current_step = 0        # índice dentro das coordenadas do path
		self.playing = False
		self.speed_ms = 600  # intervalo entre passos
		self.after_id: Optional[str] = None
		self.draw_lines_var = tk.BooleanVar(value=True)
		self.draw_triangles_var = tk.BooleanVar(value=True)  # novos triângulos entre movimentos

		self.play_btn = None
		self._build_ui()
		self._bind_keys()



	# ----------------- UI Setup -----------------
	def _build_ui(self):
		main = ttk.Frame(self.root)
		main.pack(fill=tk.BOTH, expand=True)

		controls = ttk.Frame(main)
		controls.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

		ttk.Label(controls, text="Caminho:").pack(side=tk.LEFT)
		self.path_var = tk.StringVar()
		values = [f"{p.indice}" for p in self.paths]
		self.combo = ttk.Combobox(controls, textvariable=self.path_var, values=values, width=6, state="readonly")
		self.combo.current(0)
		self.combo.pack(side=tk.LEFT, padx=4)
		self.combo.bind("<<ComboboxSelected>>", lambda e: self._on_select_path())

		# Botões de navegação
		btn_specs = [
			("⏮", self.first_step, "Primeiro"),
			("◀", self.prev_step, "Anterior"),
			("▶", self.next_step, "Próximo"),
			("⏭", self.last_step, "Último"),
		]
		for txt, cmd, tip in btn_specs:
			b = ttk.Button(controls, text=txt, width=3, command=cmd)
			b.pack(side=tk.LEFT, padx=2)
			b.tooltip = tip  # placeholder (poderia criar tooltip real)

		# Botões de seleção de colisão (em um frame abaixo dos botões de navegação)
		col_btn_frame = ttk.Frame(main)
		col_btn_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=(0, 5))
		mais_btn = ttk.Button(col_btn_frame, text="Mais colisões", width=14, command=self.selecionar_mais_colisoes)
		mais_btn.pack(side=tk.LEFT, padx=1)
		menos_btn = ttk.Button(col_btn_frame, text="Menos colisões", width=14, command=self.selecionar_menos_colisoes)
		menos_btn.pack(side=tk.LEFT, padx=1)
		menos_hipot_btn = ttk.Button(col_btn_frame, text="Menos hipotenusa", width=16, command=self.selecionar_menos_colisoes_hipotenusa)
		menos_hipot_btn.pack(side=tk.LEFT, padx=1)
		mais_hipot_btn = ttk.Button(col_btn_frame, text="Mais hipotenusa", width=16, command=self.selecionar_mais_colisoes_hipotenusa)
		mais_hipot_btn.pack(side=tk.LEFT, padx=1)
		mais_resets_btn = ttk.Button(col_btn_frame, text="Mais resets hipotenusa", width=20, command=self.selecionar_mais_resets_hipotenusa)
		mais_resets_btn.pack(side=tk.LEFT, padx=1)

		self.play_btn = ttk.Button(controls, text="Play", command=self.toggle_play)
		self.play_btn.pack(side=tk.LEFT, padx=6)

		ttk.Label(controls, text="Velocidade (ms):").pack(side=tk.LEFT, padx=(12, 2))
		self.speed_var = tk.IntVar(value=self.speed_ms)
		speed_spin = ttk.Spinbox(controls, from_=50, to=3000, textvariable=self.speed_var, width=6, increment=50, command=self._on_speed_change)
		speed_spin.pack(side=tk.LEFT)

		ttk.Checkbutton(controls, text="Linhas", variable=self.draw_lines_var, command=self.redraw).pack(side=tk.LEFT, padx=8)
		ttk.Checkbutton(controls, text="Triângulos", variable=self.draw_triangles_var, command=self.redraw).pack(side=tk.LEFT, padx=4)

		ttk.Button(controls, text="Redesenhar", command=self.redraw).pack(side=tk.LEFT, padx=4)

		self.info_var = tk.StringVar()
		info_label = ttk.Label(controls, textvariable=self.info_var)
		info_label.pack(side=tk.RIGHT, padx=8)

		# Label de estatísticas (abaixo dos controles, acima do tabuleiro)
		self.stats_var = tk.StringVar()
		stats_frame = ttk.Frame(main)
		stats_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=(0, 5))
		stats_label = ttk.Label(stats_frame, textvariable=self.stats_var, font=("Arial", 11, "bold"), foreground="#444")
		stats_label.pack(side=tk.LEFT, anchor="w")

		# Canvas
		canvas_width = self.margin * 2 + self.board_size * self.square_size
		canvas_height = self.margin * 2 + self.board_size * self.square_size
		self.canvas = tk.Canvas(main, bg="white", width=canvas_width, height=canvas_height)
		self.canvas.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

		# Caixa de texto para debug das coordenadas dos triângulos
		debug_frame = ttk.Frame(main)
		debug_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(0, 10))
		ttk.Label(debug_frame, text="Coordenadas dos triângulos (debug):").pack(side=tk.LEFT, anchor="nw")
		self.tri_debug_text = tk.Text(debug_frame, height=6, width=80, font=("Consolas", 10))
		self.tri_debug_text.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 0))
		self.tri_debug_text.configure(state="disabled")

		# Só chama _load_path após todos os widgets existirem
		self._load_path(0)
	def selecionar_mais_colisoes_hipotenusa(self):
		def hipotenusa(tri):
			from math import dist
			a, b, c = tri
			dists = [(a, b, dist(a, b)), (a, c, dist(a, c)), (b, c, dist(b, c))]
			dists.sort(key=lambda x: -x[2])
			return (dists[0][0], dists[0][1])
		def segments_cross(seg1, seg2):
			def ccw(A, B, C):
				return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
			A, B = seg1
			C, D = seg2
			return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))
		max_idx = 0
		max_col = -1
		for idx, path in enumerate(self.paths):
			triangles = self._optimized_triangles(path)
			hipots = [hipotenusa(tri) for tri, cells, alt_tri, alt_cells, opt_idx in triangles]
			colisoes = 0
			for i in range(len(hipots)):
				for j in range(i+1, len(hipots)):
					if i+1 == j or i == j+1:
						continue
					v1 = set(hipots[i])
					v2 = set(hipots[j])
					if v1 & v2:
						continue
					if segments_cross(hipots[i], hipots[j]):
						colisoes += 1
			if colisoes > max_col:
				max_col = colisoes
				max_idx = idx
		self._load_path(max_idx)
		self.combo.current(max_idx)

	def selecionar_mais_resets_hipotenusa(self):
		def estatistica_resets_colisao_hipotenusa(triangles):
			from math import dist
			def hipotenusa(tri):
				a, b, c = tri
				dists = [(a, b, dist(a, b)), (a, c, dist(a, c)), (b, c, dist(b, c))]
				dists.sort(key=lambda x: -x[2])
				return (dists[0][0], dists[0][1])
			def segments_cross(seg1, seg2):
				def ccw(A, B, C):
					return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
				A, B = seg1
				C, D = seg2
				return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))
			resets = 0
			hipots_atuais = []
			for idx, (tri, cells, alt_tri, alt_cells, opt_idx) in enumerate(triangles):
				h = hipotenusa(tri)
				colidiu = False
				for h2 in hipots_atuais:
					if set(h) & set(h2):
						continue
					if segments_cross(h, h2):
						colidiu = True
						break
				if colidiu:
					resets += 1
					hipots_atuais = []
				hipots_atuais.append(h)
			return resets
		max_idx = None
		max_resets = -1
		for idx, path in enumerate(self.paths):
			triangles = self._optimized_triangles(path)
			if not triangles:
				continue
			resets = estatistica_resets_colisao_hipotenusa(triangles)
			if resets > max_resets:
				max_resets = resets
				max_idx = idx
		if max_idx is not None:
			self._load_path(max_idx)
			self.combo.current(max_idx)
		self.root.title("Caminhos do Cavalo")
		self.root.geometry("1000x800")

		main = ttk.Frame(self.root)
		main.pack(fill=tk.BOTH, expand=True)

		controls = ttk.Frame(main)
		controls.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

		ttk.Label(controls, text="Caminho:").pack(side=tk.LEFT)
		self.path_var = tk.StringVar()
		values = [f"{p.indice}" for p in self.paths]
		self.combo = ttk.Combobox(controls, textvariable=self.path_var, values=values, width=6, state="readonly")
		self.combo.current(0)
		self.combo.pack(side=tk.LEFT, padx=4)
		self.combo.bind("<<ComboboxSelected>>", lambda e: self._on_select_path())

		# Botões de navegação
		btn_specs = [
			("⏮", self.first_step, "Primeiro"),
			("◀", self.prev_step, "Anterior"),
			("▶", self.next_step, "Próximo"),
			("⏭", self.last_step, "Último"),
		]
		for txt, cmd, tip in btn_specs:
			b = ttk.Button(controls, text=txt, width=3, command=cmd)
			b.pack(side=tk.LEFT, padx=2)
			b.tooltip = tip  # placeholder (poderia criar tooltip real)

		# Botões de seleção de colisão (após definição de controls)
		col_btn_frame = ttk.Frame(controls)
		col_btn_frame.pack(side=tk.LEFT, padx=6)
		mais_btn = ttk.Button(col_btn_frame, text="Mais colisões", command=self.selecionar_mais_colisoes)
		mais_btn.pack(side=tk.LEFT, padx=1)
		menos_btn = ttk.Button(col_btn_frame, text="Menos colisões", command=self.selecionar_menos_colisoes)
		menos_btn.pack(side=tk.LEFT, padx=1)
		menos_hipot_btn = ttk.Button(col_btn_frame, text="Menos colisão hipotenusa", command=self.selecionar_menos_colisoes_hipotenusa)
		menos_hipot_btn.pack(side=tk.LEFT, padx=1)

		self.play_btn = ttk.Button(controls, text="Play", command=self.toggle_play)
		self.play_btn.pack(side=tk.LEFT, padx=6)

		ttk.Label(controls, text="Velocidade (ms):").pack(side=tk.LEFT, padx=(12, 2))
		self.speed_var = tk.IntVar(value=self.speed_ms)
		speed_spin = ttk.Spinbox(controls, from_=50, to=3000, textvariable=self.speed_var, width=6, increment=50, command=self._on_speed_change)
		speed_spin.pack(side=tk.LEFT)

		ttk.Checkbutton(controls, text="Linhas", variable=self.draw_lines_var, command=self.redraw).pack(side=tk.LEFT, padx=8)
		ttk.Checkbutton(controls, text="Triângulos", variable=self.draw_triangles_var, command=self.redraw).pack(side=tk.LEFT, padx=4)

		ttk.Button(controls, text="Redesenhar", command=self.redraw).pack(side=tk.LEFT, padx=4)

		self.info_var = tk.StringVar()
		info_label = ttk.Label(controls, textvariable=self.info_var)
		info_label.pack(side=tk.RIGHT, padx=8)

		# Label de estatísticas (abaixo dos controles, acima do tabuleiro)
		self.stats_var = tk.StringVar()
		stats_frame = ttk.Frame(main)
		stats_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=(0, 5))
		stats_label = ttk.Label(stats_frame, textvariable=self.stats_var, font=("Arial", 11, "bold"), foreground="#444")
		stats_label.pack(side=tk.LEFT, anchor="w")


		# Canvas
		canvas_width = self.margin * 2 + self.board_size * self.square_size
		canvas_height = self.margin * 2 + self.board_size * self.square_size
		self.canvas = tk.Canvas(main, bg="white", width=canvas_width, height=canvas_height)
		self.canvas.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

		# Caixa de texto para debug das coordenadas dos triângulos
		debug_frame = ttk.Frame(main)
		debug_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(0, 10))
		ttk.Label(debug_frame, text="Coordenadas dos triângulos (debug):").pack(side=tk.LEFT, anchor="nw")
		self.tri_debug_text = tk.Text(debug_frame, height=6, width=80, font=("Consolas", 10))
		self.tri_debug_text.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(8, 0))
		self.tri_debug_text.configure(state="disabled")

		# Só chama _load_path após todos os widgets existirem
		self._load_path(0)


	def _bind_keys(self):
		self.root.bind("<space>", lambda e: self.toggle_play())
		self.root.bind("<Left>", lambda e: self.prev_step())
		self.root.bind("<Right>", lambda e: self.next_step())
		self.root.bind("<Home>", lambda e: self.first_step())
		self.root.bind("<End>", lambda e: self.last_step())
		self.root.bind("<Up>", lambda e: self._change_path(-1))
		self.root.bind("<Down>", lambda e: self._change_path(1))

	# ----------------- Data helpers -----------------
	def _infer_board_size(self) -> int:
		mx = 0
		for p in self.paths:
			for x, y in p.coords:
				mx = max(mx, x, y)
		return mx + 1

	def _on_select_path(self):
		try:
			idx = int(self.path_var.get())
		except ValueError:
			return
		for i, p in enumerate(self.paths):
			if p.indice == idx:
				self._load_path(i)
				break

	def _change_path(self, delta: int):
		new = (self.current_path_index + delta) % len(self.paths)
		self._load_path(new)
		self.combo.current(new)

	def _load_path(self, list_index: int):
		self.stop_play()
		self.current_path_index = list_index
		self.current_step = 0
		# Apenas redesenha, não reconstrói UI
		self.redraw()

	def _on_speed_change(self):
		self.speed_ms = max(10, int(self.speed_var.get()))
		if self.playing:
			self._schedule_next()

	# ----------------- Drawing -----------------
	def redraw(self):
		self.canvas.delete("all")
		self._draw_board()
		self._draw_path(self.paths[self.current_path_index])
		self._update_info()

	def _draw_board(self):
		size = self.board_size
		s = self.square_size
		for y in range(size):
			for x in range(size):
				color = "#EEE" if (x + y) % 2 == 0 else "#999"
				x0 = self.margin + x * s
				y0 = self.margin + y * s
				self.canvas.create_rectangle(x0, y0, x0 + s, y0 + s, fill=color, outline="#555")

		# eixos / labels
		for x in range(size):
			label = str(x)
			x0 = self.margin + x * s + s / 2
			self.canvas.create_text(x0, self.margin / 2, text=label, font=("Arial", 12, "bold"))
		for y in range(size):
			label = str(y)
			y0 = self.margin + y * s + s / 2
			self.canvas.create_text(self.margin / 2, y0, text=label, font=("Arial", 12, "bold"))

	def _draw_path(self, path: KnightPath):
		s = self.square_size
		# Números ordem de visita (até o passo atual para destacar?)
		for i, (x, y) in enumerate(path.coords[: self.current_step + 1]):
			x0 = self.margin + x * s + s / 2
			y0 = self.margin + y * s + s / 2
			fill = "#003366" if i < self.current_step else "#222"
			self.canvas.create_text(x0, y0, text=str(i + 1), fill=fill, font=("Arial", int(s/3), "bold"))

		# Linhas
		if self.draw_lines_var.get():
			for i in range(1, self.current_step + 1):
				x1, y1 = path.coords[i - 1]
				x2, y2 = path.coords[i]
				self._draw_line_move(x1, y1, x2, y2, highlight=(i == self.current_step))

		# Triângulos (contador entre movimentos)
		if self.draw_triangles_var.get():
			self._draw_move_triangles(path)

		# Cavalo atual
		if path.coords:
			x, y = path.coords[self.current_step]
			self._draw_knight_piece(x, y)

	def _draw_line_move(self, x1: int, y1: int, x2: int, y2: int, highlight: bool = False):
		s = self.square_size
		ax = self.margin + x1 * s + s / 2
		ay = self.margin + y1 * s + s / 2
		bx = self.margin + x2 * s + s / 2
		by = self.margin + y2 * s + s / 2
		color = "#d9534f" if highlight else "#444"
		width = 4 if highlight else 2
		self.canvas.create_line(ax, ay, bx, by, fill=color, width=width, arrow=tk.LAST, arrowshape=(12, 15, 6))

	def _draw_knight_piece(self, x: int, y: int):
		s = self.square_size
		cx = self.margin + x * s + s / 2
		cy = self.margin + y * s + s / 2
		r = s * 0.35
		self.canvas.create_oval(cx - r, cy - r, cx + r, cy + r, fill="#ffcc00", outline="#b8860b", width=2)
		self.canvas.create_text(cx, cy, text="♞", font=("Segoe UI Symbol", int(s/2), "bold"))

	def _draw_move_triangles(self, path: KnightPath):
		"""Desenha triângulos otimizados para cada salto, minimizando colisões."""
		if self.current_step == 0:
			return
		s = self.square_size
		triangles = self._optimized_triangles(path)
		for i, (triangle, cells, alt_triangle, alt_cells, opt_idx) in enumerate(triangles[:self.current_step]):
			def cell_center(cx: int, cy: int):
				return (
					self.margin + cx * s + s / 2,
					self.margin + cy * s + s / 2,
				)
			# Só desenha alternativo se for diferente do escolhido
			if alt_triangle != triangle:
				x1a, y1a = alt_triangle[0]
				mida = alt_triangle[1]
				x2a, y2a = alt_triangle[2]
				axa, aya = cell_center(x1a, y1a)
				bxa, bya = cell_center(*mida)
				cxa, cya = cell_center(x2a, y2a)
				self.canvas.create_polygon(
					(axa, aya, bxa, bya, cxa, cya),
					fill='',
					outline="#1e90ff",
					width=2,
					dash=(6, 4),
					joinstyle=tk.MITER
				)
			# Desenha triângulo escolhido (cheio)
			x1, y1 = triangle[0]
			mid = triangle[1]
			x2, y2 = triangle[2]
			ax, ay = cell_center(x1, y1)
			bx, by = cell_center(*mid)
			cx, cy = cell_center(x2, y2)
			fill = "#5f6a72" if i < self.current_step - 1 else "#ff8c1a"
			outline = "#394247" if i < self.current_step - 1 else "#c25e00"
			self.canvas.create_polygon(
				(ax, ay, bx, by, cx, cy),
				fill=fill,
				outline=outline,
				width=2 if i == self.current_step - 1 else 1,
				joinstyle=tk.MITER
			)
			# Número do salto (i+1) no centróide
			centroid_x = (ax + bx + cx) / 3
			centroid_y = (ay + by + cy) / 3
			self.canvas.create_text(
				centroid_x,
				centroid_y,
				text=str(i+1),
				fill="white" if i == self.current_step - 1 else "#f0f0f0",
				font=("Arial", int(s/5), "bold")
			)


	def _update_info(self):
		# Estatística de resets de colisão de hipotenusa
		def estatistica_resets_colisao_hipotenusa(triangles):
			from math import dist
			def hipotenusa(tri):
				a, b, c = tri
				dists = [(a, b, dist(a, b)), (a, c, dist(a, c)), (b, c, dist(b, c))]
				dists.sort(key=lambda x: -x[2])
				return (dists[0][0], dists[0][1])
			def segments_cross(seg1, seg2):
				def ccw(A, B, C):
					return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
				A, B = seg1
				C, D = seg2
				return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))
			resets = 0
			hipots_atuais = []
			for idx, (tri, cells, alt_tri, alt_cells, opt_idx) in enumerate(triangles):
				h = hipotenusa(tri)
				colidiu = False
				for h2 in hipots_atuais:
					if set(h) & set(h2):
						continue
					if segments_cross(h, h2):
						colidiu = True
						break
				if colidiu:
					resets += 1
					hipots_atuais = []
				hipots_atuais.append(h)
			return resets
		p = self.paths[self.current_path_index]
		info = [f"Caminho {p.indice}"]
		if p.is_full_cover(self.board_size):
			info.append("(cobre todo tabuleiro)")
		if p.is_closed_tour():
			info.append("(tour fechado)")
		info.append(f"Passo {self.current_step + 1}/{len(p.coords)}")
		self.info_var.set("  |  ".join(info))

		# Estatísticas organizadas
		countmap = self._triangles_covered_cells_countmap(p)
		triangles = self._optimized_triangles(p)
		triangles_cells = [cells for tri, cells, alt_tri, alt_cells, opt_idx in triangles]
		# Colisão: número de pares de triângulos que se sobrepõem em pelo menos uma casa
		colisoes = 0
		for i in range(len(triangles_cells)):
			for j in range(i+1, len(triangles_cells)):
				if triangles_cells[i] & triangles_cells[j]:
					colisoes += 1

		# Colisão de hipotenusas: apenas se a hipotenusa de um triângulo cruza a de outro
		def hipotenusa(tri):
			# A hipotenusa é o lado entre os dois vértices mais distantes
			from math import dist
			a, b, c = tri
			dists = [(a, b, dist(a, b)), (a, c, dist(a, c)), (b, c, dist(b, c))]
			dists.sort(key=lambda x: -x[2])
			return (dists[0][0], dists[0][1])

		def segments_cross(seg1, seg2):
			# Verifica se os segmentos (p1, p2) e (q1, q2) se cruzam
			def ccw(A, B, C):
				return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
			A, B = seg1
			C, D = seg2
			return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))

		colisoes_hipotenusa = 0
		hipots = [hipotenusa(tri) for tri, cells, alt_tri, alt_cells, opt_idx in triangles]
		for i in range(len(hipots)):
			for j in range(i+1, len(hipots)):
				if segments_cross(hipots[i], hipots[j]):
					colisoes_hipotenusa += 1

		total_casas = self.board_size * self.board_size
		cobertos = set(k for k, v in countmap.items() if v > 0)
		n_nao_coberto = total_casas - len(cobertos)
		percent_nao_coberto = 100 * n_nao_coberto / total_casas if total_casas else 0

		resets_hipot = estatistica_resets_colisao_hipotenusa(triangles)
		stats_text = (
			f"Colisões de triângulos: {colisoes}    "
			f"Colisões de hipotenusas: {colisoes_hipotenusa}    "
			f"Resets de colisão hipotenusa: {resets_hipot}    "
			f"Casas não cobertas: {n_nao_coberto} de {total_casas}  ({percent_nao_coberto:.1f}%)"
		)
		self.stats_var.set(stats_text)

		# Atualiza caixa de debug com coordenadas dos triângulos escolhidos
		lines = []
		for i, (triangle, cells, alt_triangle, alt_cells, opt_idx) in enumerate(triangles[:self.current_step]):
			a, b, c = triangle
			lines.append(f"{i+1}: {a} {b} {c}")
		debug_str = "\n".join(lines)
		self.tri_debug_text.configure(state="normal")
		self.tri_debug_text.delete("1.0", tk.END)
		self.tri_debug_text.insert("1.0", debug_str)
		self.tri_debug_text.configure(state="disabled")

	# ----------------- Navegação passos -----------------
	def first_step(self):
		self.current_step = 0
		self.stop_play()
		self.redraw()

	def last_step(self):
		p = self.paths[self.current_path_index]
		if p.coords:
			self.current_step = len(p.coords) - 1
		self.stop_play()
		self.redraw()

	def prev_step(self):
		if self.current_step > 0:
			self.current_step -= 1
			self.stop_play()
			self.redraw()

	def next_step(self):
		p = self.paths[self.current_path_index]
		if self.current_step < len(p.coords) - 1:
			self.current_step += 1
			self.redraw()
		else:
			self.stop_play()

	# ----------------- Animação -----------------
	def toggle_play(self):
		if self.playing:
			self.stop_play()
		else:
			self.start_play()

	def start_play(self):
		self.playing = True
		if self.play_btn:
			self.play_btn.configure(text="Pause")
		self._schedule_next(immediate=True)

	def stop_play(self):
		self.playing = False
		if self.play_btn:
			self.play_btn.configure(text="Play")
		if self.after_id:
			self.root.after_cancel(self.after_id)
			self.after_id = None

	def _schedule_next(self, immediate: bool = False):
		if not self.playing:
			return
		delay = 10 if immediate else self.speed_ms
		if self.after_id:
			self.root.after_cancel(self.after_id)
		self.after_id = self.root.after(delay, self._advance_animation)

	def _advance_animation(self):
		if not self.playing:
			return
		p = self.paths[self.current_path_index]
		if self.current_step < len(p.coords) - 1:
			self.current_step += 1
			self.redraw()
			self._schedule_next()
		else:
			self.stop_play()


def main():
	try:
		paths = parse_paths(CAMINHO_ARQ)
	except Exception as e:  # pragma: no cover - ambiente runtime
		print(f"Erro ao ler caminhos: {e}")
		sys.exit(1)

	root = tk.Tk()
	viewer = KnightPathsViewer(root, paths)
	root.mainloop()


if __name__ == "__main__":
	main()

"""Notas e melhorias futuras possíveis:
1. Exportar sequência para GIF ou MP4 (usar pillow + imageio) 
2. Salvar imagem estática do caminho atual.
3. Mode 'heatmap' (contagem de visitas se vários caminhos sobrepostos).
4. Filtro para mostrar somente caminhos que cobrem todo o tabuleiro.
5. Mostrar estatísticas agregadas (quantos tours completos, quantos fechados, etc.).
6. Implementar zoom e redimensionamento responsivo do canvas.
7. Adicionar tooltip real nos botões.
"""

