// ============================================================
// Knight Tour Wiki — nav.js
// Renderiza sidebar a partir de PAGES[]; marca a página ativa
// e exporta um botão hamburguer para mobile.
// Vanilla JS, sem dependências.
// ============================================================

const PAGES = [
  { id: 'index',          title: 'Início',                       file: 'index.html' },
  { id: '01-introducao',  title: '1. Introdução',                file: '01-introducao.html' },
  { id: '02-teoria-gf2',  title: '2. GF(2) e Espaço de Ciclos',  file: '02-teoria-gf2.html' },
  { id: '03-teorema-q3',  title: '3. Teorema Q(n) = 3',          file: '03-teorema-q3.html' },
  { id: '04-algoritmo',   title: '4. O Algoritmo',               file: '04-algoritmo.html' },
  { id: '05-fase-local',  title: '5. Fase Local f∞(L)',          file: '05-fase-local.html' },
  { id: '06-benchmark',   title: '6. Benchmark',                 file: '06-benchmark.html' },
  { id: '07-estimativas', title: '7. Estimativas N(n)',          file: '07-estimativas.html' },
  { id: '08-topologia',   title: '8. Análise Topológica',        file: '08-topologia.html' },
  { id: '09-extensoes',   title: '9. Extensões',                 file: '09-extensoes.html' },
  { id: '10-resultados',  title: '10. Todos os Resultados',      file: '10-resultados.html' },
  { id: '11-codigo',      title: '11. Como Usar o Código',       file: '11-codigo.html' },
];

function currentPageId() {
  const file = (location.pathname.split('/').pop() || 'index.html').toLowerCase();
  const match = PAGES.find(p => p.file.toLowerCase() === file);
  return match ? match.id : 'index';
}

function renderSidebar() {
  const sidebar = document.getElementById('w-sidebar');
  if (!sidebar) return;
  const active = currentPageId();

  sidebar.innerHTML = `
    <div class="w-sidebar-title">Knight Tour Wiki</div>
    <ul class="w-sidebar-list">
      ${PAGES.map(p => `
        <li>
          <a href="${p.file}"${p.id === active ? ' class="w-active"' : ''}>
            ${p.title}
          </a>
        </li>
      `).join('')}
    </ul>
    <div class="w-sidebar-footer">
      <a href="../index.html">← Visualizador</a>
    </div>
  `;
}

function setupHamburger() {
  const btn = document.getElementById('w-hamburger');
  const sidebar = document.getElementById('w-sidebar');
  if (!btn || !sidebar) return;
  btn.addEventListener('click', () => {
    sidebar.classList.toggle('w-open');
  });
  document.addEventListener('click', (e) => {
    if (window.innerWidth > 860) return;
    if (!sidebar.contains(e.target) && e.target !== btn) {
      sidebar.classList.remove('w-open');
    }
  });
}

function renderPrevNext() {
  const footer = document.getElementById('w-nav-footer');
  if (!footer) return;
  const id = currentPageId();
  const idx = PAGES.findIndex(p => p.id === id);
  const prev = idx > 0 ? PAGES[idx - 1] : null;
  const next = idx >= 0 && idx < PAGES.length - 1 ? PAGES[idx + 1] : null;

  footer.innerHTML = `
    ${prev ? `<a class="w-prev" href="${prev.file}">← ${prev.title}</a>` : '<span></span>'}
    ${next ? `<a class="w-next" href="${next.file}">${next.title} →</a>` : ''}
  `;
}

document.addEventListener('DOMContentLoaded', () => {
  renderSidebar();
  setupHamburger();
  renderPrevNext();
});
