'use strict';
// Allele — offline DNA report. Analysis core first (pure, tested by test.js), then the interface.

const VERSION = '1.0.2';
const EVF = { strong: 1, moderate: 0.75, limited: 0.5 };
const SIG = { P: 'Pathogenic', LP: 'Likely pathogenic', RF: 'Risk factor', DR: 'Drug response', PR: 'Protective' };
const HUES = [205, 25, 145, 330, 85, 265];
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const hasClinvar = () => typeof CLINVAR !== 'undefined';

// ---------- parsing ----------
function parseLines(lines, from, to, geno, meta) {
  for (let i = from; i < to; i++) {
    const line = lines[i];
    if (!line) continue;
    if (line[0] === '#') { meta.push(line); continue; }
    const p = line.replace(/"/g, '').split(/[\t,]+|\s+/);
    if (!/^(rs|i)\d+$/i.test(p[0])) continue;
    let g = (p.length >= 5 && p[3].length === 1 && p[4].length === 1) ? p[3] + p[4] : (p[3] || '');
    g = g.toUpperCase();
    if (!/^[ACGT]{1,2}$/.test(g)) continue; // no-calls (--, 00), indels (DI)
    geno.set(p[0].toLowerCase(), { g: g.length === 1 ? g + g : g, chr: p[1], pos: +p[2] });
  }
}
function parseRaw(text) {
  const lines = text.split(/\r?\n/), geno = new Map(), meta = [];
  parseLines(lines, 0, lines.length, geno, meta);
  return { geno, meta };
}

// ---------- scoring ----------
// Score of a genotype for a *harmful* effect allele: +1 best ... -1 worst. Beneficial alleles use the negation.
function harmVal(e, k) {
  if (e.v) return e.v[k];
  if (e.model === 'dom') return [1, -1, -1][k];
  if (e.model === 'rec') return [1, 0.5, -1][k];
  return [1, 0, -1][k];
}
const markerVal = (e, w, k) => (w < 0 ? 1 : -1) * harmVal(e, k);

function lookup(e, geno) {
  for (const s of [e, ...(e.proxies || [])]) {
    const r = geno.get(s.id);
    if (!r) continue;
    const ok = [...r.g].every(c => c === s.alt || c === s.oth);
    if (!ok) return { id: s.id, g: r.g, mismatch: true, proxy: s !== e };
    return { id: s.id, g: r.g, k: [...r.g].filter(c => c === s.alt).length, proxy: s !== e, alt: s.alt };
  }
  return null;
}

function genoText(e, k) {
  if (e.fx) return e.fx[k];
  const harmful = Object.values(e.cats).some(w => w < 0);
  if (k === 0) return harmful ? `No copies of the effect allele (${e.alt}) — typical or lower risk for this marker.` : `Does not carry the beneficial ${e.alt} allele.`;
  return (k === 1 ? `One copy of ${e.alt}: ` : `Two copies of ${e.alt}: `) + e.eff;
}

function analyze(geno) {
  const findings = [], cats = {};
  for (const key in CATEGORIES) cats[key] = { key, ...CATEGORIES[key], sum: 0, wsum: 0, total: 0, items: [], missing: [] };
  for (const e of SNP_DB) {
    const r = lookup(e, geno);
    const f = { e, r };
    findings.push(f);
    for (const [c, w] of Object.entries(e.cats)) {
      cats[c].total++;
      if (!r || r.mismatch) { cats[c].missing.push(f); continue; }
      const val = markerVal(e, w, r.k);
      const wt = Math.abs(w) * EVF[e.ev];
      cats[c].sum += wt * val; cats[c].wsum += wt;
      cats[c].items.push({ f, val, w });
    }
  }
  for (const c of Object.values(cats)) {
    c.score = c.wsum ? 50 + 50 * c.sum / c.wsum : null;
    c.items.sort((a, b) => a.val - b.val || Math.abs(b.w) - Math.abs(a.w));
  }
  return { findings, cats: Object.values(cats) };
}

function apoeType(geno) {
  const a = geno.get('rs429358'), b = geno.get('rs7412');
  if (!b) return null;
  // ε2 = T-T (rs429358 T, rs7412 T), ε3 = T-C, ε4 = C-C. Phase is assumed (ε1 is extremely rare).
  const tag = geno.get('rs769449');
  if (!a && !tag) return { text: 'undetermined', note: `rs429358 missing; rs7412 = ${b.g}`, approx: true };
  const n4 = a ? [...a.g].filter(c => c === 'C').length : [...tag.g].filter(c => c === 'A').length;
  const n2 = [...b.g].filter(c => c === 'T').length;
  const alleles = [...Array(n4).fill('ε4'), ...Array(n2).fill('ε2')];
  while (alleles.length < 2) alleles.push('ε3');
  const text = alleles.slice(0, 2).sort().join('/');
  return a ? { text, approx: false } : { text, note: `rs429358 missing — ε4 status from tag SNP rs769449 (${tag.g})`, approx: true };
}

function clinvarHits(geno) {
  const hits = [];
  if (!hasClinvar()) return hits;
  for (const [id, r] of geno) {
    const recs = CLINVAR[id];
    if (!recs) continue;
    for (const [gene, sig, cond, ref, alt, stars, vid] of recs) {
      if (![...r.g].every(c => c === ref || c === alt)) continue;
      const k = [...r.g].filter(c => c === alt).length;
      if (k) hits.push({ id, g: r.g, gene: gene && gene !== '-' ? gene : '', sig, cond, ref, alt, stars, vid, k });
    }
  }
  const order = { P: 0, LP: 1, RF: 2, DR: 3, PR: 4 };
  return hits.sort((a, b) => order[a.sig] - order[b.sig] || b.stars - a.stars || b.k - a.k);
}

// Identity-by-state over every autosomal SNP two files share. IBS0 = opposite homozygotes (no allele shared).
// Sex chromosomes are skipped: a son's X comes only from his mother, so father and son look unrelated there.
// kinship = KING-robust estimator (Manichaikul et al. 2010): 0.5 identical, 0.25 first degree, 0.125 second degree.
const NON_AUTOSOMAL = new Set(['X', 'Y', 'XY', 'MT', 'M', '23', '24', '25', '26']);
function compareGenomes(a, b) {
  let n = 0, same = 0, ibs0 = 0, shared = 0, hetBoth = 0, hetA = 0, hetB = 0;
  const [small, big] = a.size <= b.size ? [a, b] : [b, a];
  for (const [id, r] of small) {
    const s = big.get(id);
    if (!s || NON_AUTOSOMAL.has(String(r.chr).toUpperCase())) continue;
    const x = r.g, y = s.g, hx = x[0] !== x[1], hy = y[0] !== y[1];
    n++; hetA += hx; hetB += hy;
    if ((x[0] === y[0] && x[1] === y[1]) || (x[0] === y[1] && x[1] === y[0])) { same++; shared += 2; hetBoth += hx; continue; }
    const rest = [...y];
    let k = 0;
    for (const c of x) { const i = rest.indexOf(c); if (i >= 0) { k++; rest.splice(i, 1); } }
    shared += k;
    if (!k) ibs0++;
  }
  const kinship = hetA + hetB ? (hetBoth - 2 * ibs0) / (hetA + hetB) : 0;
  return { n, concordance: n ? same / n : 0, ibs0: n ? ibs0 / n : 0, ibs: n ? shared / (2 * n) : 0, kinship };
}

function relationReading(c) {
  if (c.n < 50000) return 'Too few shared SNPs to judge';
  if (c.kinship > 0.354) return 'Same person or identical twin';
  if (c.kinship > 0.177) return c.ibs0 < 0.003 ? 'Parent and child' : 'Full siblings likely';
  if (c.kinship > 0.088) return 'Second-degree relatives likely (grandparent, aunt or uncle, half-sibling)';
  if (c.kinship > 0.044) return 'Third-degree relatives possible (e.g. first cousins)';
  return 'No close relationship detected';
}

// ---------- read check: genotype calls vs the sequencing reads in a CRAM ----------
// Sites are the person's ClinVar hits and curated markers, at GRCh37 positions with a known reference base.
function readSites(p) {
  const sites = new Map();
  for (const h of p.hits) {
    const r = p.geno.get(h.id);
    if (r && r.pos) sites.set(`${r.chr}:${r.pos}`, { id: h.id, chr: r.chr, pos: r.pos, ref: h.ref });
  }
  for (const f of p.found) {
    const snp = [f.e, ...(f.e.proxies || [])].find(s => s.id === f.r.id);
    if (!snp || !snp.b37) continue;
    const [chr, pos, ref] = snp.b37.split(':');
    sites.set(`${chr}:${pos}`, { id: snp.id, chr, pos: +pos, ref });
  }
  return [...sites.values()];
}

// Base a read carries at 1-based reference position pos, and its index in the read for the quality lookup.
// CRAM read features hold 0-based refPos; a position with no feature carries the reference base.
function baseAt(features, readStart0, pos, refBase) {
  let base = refBase, shift = 0;
  for (const ft of features) {
    const fp = ft.refPos + 1;
    if (fp > pos) break;
    if (ft.code === 'D' || ft.code === 'N') { if (pos < fp + ft.data) return null; shift -= ft.data; continue; }
    if (ft.code === 'I' || ft.code === 'S') { shift += ft.data.length; continue; } // inserted bases sit before refPos
    if (ft.code === 'i') { shift += 1; continue; }
    if (fp === pos && ft.code === 'X') base = ft.sub || 'N';
    if (fp === pos && ft.code === 'B') base = ft.data[0];
  }
  return { base, index: pos - readStart0 - 1 + shift };
}

// What the reads say about the file's genotype g at a site with reference base ref; obs = base counts.
function readVerdict(g, ref, obs) {
  const depth = 'ACGT'.split('').reduce((s, b) => s + (obs[b] || 0), 0);
  const alts = [...new Set(g)].filter(b => b !== ref);
  const altSeen = alts.reduce((s, b) => s + (obs[b] || 0), 0), refSeen = obs[ref] || 0;
  const other = depth - altSeen - refSeen;
  if (!depth) return { key: 'none', label: 'No reads', depth };
  if (alts.length) {
    if (altSeen) return { key: 'seen', label: 'Seen in reads', depth, altSeen, refSeen };
    return { key: 'missing', label: `Not seen in ${depth} read${depth > 1 ? 's' : ''}`, depth, altSeen, refSeen };
  }
  if (other) return { key: 'conflict', label: 'Reads show another allele', depth, altSeen: other, refSeen };
  return { key: 'agree', label: 'Matches reads', depth, altSeen, refSeen };
}

async function runReadCheck(blob, sites, CramFile, BlobFile, onProgress = () => {}) {
  const file = new CramFile({
    filehandle: new BlobFile(blob),
    // Only the bases at the checked sites are needed; everything else reads as N.
    fetchReferenceSequence: async (seqId, start, end) => {
      const s = new Array(end - start).fill('N');
      for (const x of bySeq.get(seqId) || []) if (x.pos - 1 >= start && x.pos - 1 < end) s[x.pos - 1 - start] = x.ref;
      return s.join('');
    },
  });
  const sq = (await file.getSamHeader()).filter(l => l.tag === 'SQ').map(l => Object.fromEntries(l.data.map(d => [d.tag, d.value])));
  const chr1 = sq.find(s => /^(chr)?1$/i.test(s.SN));
  if (!chr1) throw new Error('This CRAM has no chromosome 1 in its header, so it cannot be matched to genotype positions.');
  if (+chr1.LN !== 249250621) throw new Error(`This CRAM is aligned to ${+chr1.LN === 248956422 ? 'GRCh38' : 'an unknown reference'}; the read check needs GRCh37 (build 37) alignments, like MyHeritage's.`);
  const bySeq = new Map();
  sq.forEach((s, i) => {
    const list = sites.filter(x => x.chr.toUpperCase() === s.SN.replace(/^chr/i, '').toUpperCase()).sort((a, b) => a.pos - b.pos);
    if (list.length) bySeq.set(i, list);
  });
  const obs = new Map(sites.map(x => [x.id, {}]));
  const anyIn = (list, lo, hi) => { let a = 0, b = list.length; while (a < b) { const m = (a + b) >> 1; if (list[m].pos < lo) a = m + 1; else b = m; } return a < list.length && list[a].pos <= hi ? a : -1; };
  let containers = 0, decoded = 0, reads = 0;
  try {
    for await (const c of file.iterContainers()) {
      const h = await c.getHeader();
      containers++;
      onProgress(Math.min(1, h._endPosition / blob.size));
      if (containers === 1) continue; // the first container holds the SAM header
      if (!h.numRecords || h.refSeqId === -1) break; // EOF marker, or the unmapped reads sorted last
      const list = bySeq.get(h.refSeqId);
      if (h.refSeqId >= 0 && (!list || anyIn(list, h.refSeqStart - 1000, h.refSeqStart + h.alignmentSpan + 1000) < 0)) continue;
      decoded++;
      for (const lm of h.landmarks) {
        for (const r of await c.getSlice(lm).getRecords(() => true)) {
          if (r.flags & 0xF04 || (r.mappingQuality ?? 60) < 20) continue;
          const sl = bySeq.get(r.sequenceId);
          if (!sl) continue;
          let i = anyIn(sl, r.start + 1, r.end);
          if (i < 0) continue;
          reads++;
          const feats = r.readFeatures || [], q = r.qualityScores;
          for (; i < sl.length && sl[i].pos <= r.end; i++) {
            const x = sl[i], hit = baseAt(feats, r.start, x.pos, x.ref);
            if (!hit || !'ACGT'.includes(hit.base) || (q && q[hit.index] < 20)) continue;
            const o = obs.get(x.id);
            o[hit.base] = (o[hit.base] || 0) + 1;
          }
        }
      }
    }
  } catch (err) {
    if (!containers) throw err; // reading past the last container is how the walk ends
  }
  onProgress(1);
  return { obs, containers, decoded, reads };
}

function personName(file) {
  const f = file.replace(/\.[^.]+$/, '');
  const m = f.match(/MyHeritage_raw_dna_data_(.+)/i) || f.match(/^genome_(.+?)_v\d/i);
  return (m ? m[1] : f).replace(/[_-]+/g, ' ').trim() || 'Unnamed';
}

function buildPerson(name, file, data, hue) {
  const method = (data.meta.find(m => /method=/i.test(m)) || '').split('=')[1] || '';
  const build = (data.meta.find(m => /reference=|build/i.test(m)) || '').match(/(build|GRCh)\s*=?\s*3[78]|3[78]/i);
  const a = analyze(data.geno);
  const p = { name, file, geno: data.geno, method, build: build ? (build[0].includes('38') ? 'GRCh38' : 'GRCh37') : '', hue, ...a };
  p.found = a.findings.filter(f => f.r && !f.r.mismatch);
  p.mismatch = a.findings.filter(f => f.r && f.r.mismatch);
  p.apoe = apoeType(data.geno);
  p.notable = [];
  for (const f of p.found) {
    for (const [, w] of Object.entries(f.e.cats)) {
      const val = markerVal(f.e, w, f.r.k);
      if (Math.abs(w) >= 2 && f.e.ev !== 'limited' && ((w < 0 && val <= -0.5) || (w > 0 && f.r.k > 0))) { p.notable.push({ f, val, w }); break; }
    }
  }
  p.notable.sort((x, y) => x.val - y.val);
  refreshClinvar(p);
  return p;
}
function refreshClinvar(p) {
  p.hits = clinvarHits(p.geno);
  p.plp = p.hits.filter(h => h.sig === 'P' || h.sig === 'LP');
  p.genes = new Set([...p.found.map(f => f.e.gene), ...p.hits.map(h => h.gene).filter(Boolean)]);
}

// ======================================================================
// Interface
// ======================================================================
if (typeof window !== 'undefined') (() => {
  const $ = s => document.querySelector(s);
  const reduce = matchMedia('(prefers-reduced-motion: reduce)');
  const state = { people: [], cur: 0, cvSig: 'all', cvStars: 0, cvQuery: '', cmpMode: 'diff', cmpSel: null };
  const P = () => state.people[state.cur];
  const fmt = n => n.toLocaleString('en-US');
  const pctTxt = s => s === null ? '—' : Math.round(s) + '%';

  // ---------- small renderers ----------
  const initial = p => `<span class="initial" style="--h:${p.hue}" aria-hidden="true">${esc(p.name[0].toUpperCase())}</span>`;
  function gtCells(g, alt, diff) {
    return `<span class="gt${diff ? ' diff' : ''}" aria-label="genotype ${g}">${[...g].map(c => `<i class="${c === alt ? 'e' : ''}">${c}</i>`).join('')}</span>`;
  }
  const impactPill = v => v > 0.3 ? '<span class="pill good">Favourable</span>' : v < -0.3 ? '<span class="pill bad">Unfavourable</span>' : '<span class="pill">Intermediate</span>';
  const evidence = ev => ({ strong: 'Strong', moderate: 'Moderate', limited: 'Limited' }[ev]);
  function links(e, id) {
    const l = [`<a href="https://www.ncbi.nlm.nih.gov/snp/${id}" target="_blank" rel="noopener">dbSNP</a>`,
      `<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term=${id}" target="_blank" rel="noopener">ClinVar</a>`,
      `<a href="https://www.ebi.ac.uk/gwas/variants/${id}" target="_blank" rel="noopener">GWAS Catalog</a>`];
    if (e && e.cats.pharma) l.push(`<a href="https://www.pharmgkb.org/variant/${id}" target="_blank" rel="noopener">PharmGKB</a>`);
    if (e) l.push(`<a href="https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(e.gene.split(/[ /(]/)[0])}%5Bsym%5D+AND+human%5Borgn%5D" target="_blank" rel="noopener">Gene</a>`);
    return `<span class="lk">${l.join('')}</span>`;
  }
  const refHtml = ref => esc(ref).replace(/PMID (\d+)/g, '<a href="https://pubmed.ncbi.nlm.nih.gov/$1/" target="_blank" rel="noopener">PMID $1</a>');

  // The person whose tables are being drawn: the current one, or the one a report is exported for.
  let RP = null;
  const rp = () => RP || P();
  const READ_CLS = { seen: 'good', missing: 'warn', conflict: 'bad', agree: '', none: '' };
  function readVerdictFor(p, id, g) {
    const rc = p.readCheck, o = rc && rc.obs.get(id);
    return o ? readVerdict(g, rc.ref.get(id), o) : null;
  }
  function readsTd(id, g) {
    const p = rp();
    if (!p.readCheck) return '';
    const v = readVerdictFor(p, id, g);
    if (!v) return '<td><span class="muted">Not checked</span></td>';
    const o = p.readCheck.obs.get(id), counts = Object.entries(o).sort((a, b) => b[1] - a[1]).map(([b, n]) => `${b} ${n}`).join(' · ');
    return `<td><span class="pill ${READ_CLS[v.key]}">${esc(v.label)}</span>${counts ? `<span class="gt-note mono">${counts}</span>` : ''}</td>`;
  }

  function markerRow(f, val, w) {
    const { e, r } = f;
    const label = val === undefined ? '' : impactPill(w > 0 && r.k > 0 && val >= 0 ? 1 : val);
    return `<tr>
      <td class="gene"><b>${esc(e.gene)}</b><span class="id">${esc(e.id)}</span></td>
      <td>${gtCells(r.g, r.alt)}<span class="gt-note">effect allele ${esc(r.alt)}${r.proxy ? ` · via ${esc(r.id)}` : ''}</span></td>
      <td class="effect"><span class="t">${esc(e.trait)}</span>${esc(genoText(e, r.k))}</td>
      <td>${evidence(e.ev)}</td>
      <td>${label}</td>${readsTd(r.id, r.g)}
      <td class="src">${refHtml(e.ref)}<br>${links(e, r.id)}</td></tr>`;
  }
  const readsTh = () => rp().readCheck ? '<th>Sequencing reads</th>' : '';
  const markerHead = () => `<table class="rows"><thead><tr><th>Gene</th><th>Genotype</th><th>Effect</th><th>Evidence</th><th>Impact</th>${readsTh()}<th>Sources</th></tr></thead><tbody>`;

  // The score ladder: one rung per category; grey ticks = markers, blue dot = category score.
  function ladder(cats, opts = {}) {
    const rows = cats.map((c, i) => {
      // markers sharing a value fan out 4px apart, so a cluster reads as a comb rather than one tick
      const seen = {}, total = {};
      c.items.forEach(it => total[it.val] = (total[it.val] || 0) + 1);
      const ticks = c.items.map(it => {
        const j = seen[it.val] = (seen[it.val] ?? -1) + 1, off = (j - (total[it.val] - 1) / 2) * 4;
        return `<span class="slide" style="--to:${50 + 50 * it.val}%"><span class="tick" style="margin-left:${off}px;transition-delay:${200 + i * 35}ms"></span></span>`;
      }).join('');
      return `<a class="rung" href="#/category/${c.key}" data-cat="${c.key}">
        <span class="name">${esc(c.name)}</span>
        <span class="track">${[25, 75].map(q => `<span class="q" style="left:${q}%"></span>`).join('')}${ticks}
          ${c.score === null ? '' : `<span class="slide" style="--to:${c.score}%"><span class="dot"></span></span>`}</span>
        <span class="pct${c.score === null ? ' na' : ''}" data-count="${c.score === null ? '' : Math.round(c.score)}">${pctTxt(c.score)}</span></a>`;
    }).join('');
    return `<div class="ladder${opts.static ? ' on' : ''}">
      <div class="ladder-axis"><span></span><span class="scale"><span>0% worst possible</span><span>50% balanced</span><span>100% best</span></span><span></span></div>
      ${rows}</div>`;
  }
  function compareLadder(people) {
    const keys = Object.keys(CATEGORIES);
    const rows = keys.map(k => {
      const sc = people.map(p => p.cats.find(c => c.key === k).score);
      const valid = sc.filter(s => s !== null);
      const spread = valid.length > 1 ? Math.round(Math.max(...valid) - Math.min(...valid)) : null;
      const dots = people.map((p, i) => sc[i] === null ? '' : `<span class="slide" style="--to:${sc[i]}%"><span class="pdot" style="--h:${p.hue};margin-left:${sc.slice(0, i).filter(s => s !== null && Math.abs(s - sc[i]) < 3).length * 14}px" title="${esc(p.name)} ${Math.round(sc[i])}%">${esc(p.name[0].toUpperCase())}</span></span>`).join('');
      return `<a class="rung" href="#/category/${k}" data-cmp="${k}">
        <span class="name">${esc(CATEGORIES[k].name)}</span>
        <span class="track">${[25, 75].map(q => `<span class="q" style="left:${q}%"></span>`).join('')}${dots}</span>
        <span class="pct${spread === null ? ' na' : ''}" title="Gap between highest and lowest">${spread === null ? '—' : 'Δ ' + spread}</span></a>`;
    }).join('');
    return `<div class="ladder"><div class="ladder-axis"><span></span><span class="scale"><span>0%</span><span>50%</span><span>100%</span></span><span class="num" style="text-align:right">gap</span></div>${rows}</div>`;
  }

  // ---------- motion: one orchestrated moment per view ----------
  function animate(root) {
    const slides = root.querySelectorAll('.ladder .slide');
    const counts = root.querySelectorAll('[data-count]');
    if (reduce.matches) {
      slides.forEach(s => s.style.setProperty('--x', s.style.getPropertyValue('--to')));
      root.querySelectorAll('.ladder').forEach(l => l.classList.add('on'));
      return;
    }
    root.querySelectorAll('.ladder').forEach(l => l.classList.add('live'));
    requestAnimationFrame(() => requestAnimationFrame(() => {
      root.querySelectorAll('.ladder .rung').forEach((rung, i) => rung.querySelectorAll('.slide').forEach(s => {
        if (!s.querySelector('.tick')) s.style.transitionDelay = `${i * 40}ms`;
        s.style.setProperty('--x', s.style.getPropertyValue('--to'));
      }));
      root.querySelectorAll('.ladder').forEach(l => l.classList.add('on'));
    }));
    counts.forEach(el => {
      const end = parseFloat(el.dataset.count);
      if (isNaN(end)) return;
      const suffix = el.dataset.suffix ?? '%', t0 = performance.now(), dur = 900;
      const step = t => {
        const k = Math.min(1, (t - t0) / dur), e = 1 - Math.pow(1 - k, 3);
        const dec = +(el.dataset.dec || 0);
        el.textContent = (dec ? (end * e).toFixed(dec) : fmt(Math.round(end * e))) + suffix;
        if (k < 1) requestAnimationFrame(step);
      };
      requestAnimationFrame(step);
    });
  }

  // ---------- popover on ladder rungs ----------
  const pop = $('#pop');
  function showPop(rung, x, y) {
    const p = P();
    let html = '';
    if (rung.dataset.cat) {
      const c = p.cats.find(c => c.key === rung.dataset.cat);
      const items = [...c.items].sort((a, b) => Math.abs(b.w) * EVF[b.f.e.ev] - Math.abs(a.w) * EVF[a.f.e.ev]).slice(0, 6);
      html = `<b>${esc(c.name)} — ${pctTxt(c.score)}</b><ul>${items.map(it => `<li><span>${esc(it.f.e.gene)} <span class="mono">${esc(it.f.r.g)}</span></span><span>${it.val > 0.3 ? 'favourable' : it.val < -0.3 ? 'unfavourable' : 'intermediate'}</span></li>`).join('')}</ul>
        <p class="small muted" style="margin:8px 0 0">${c.items.length} of ${c.total} markers found · click for details</p>`;
    } else {
      const k = rung.dataset.cmp;
      html = `<b>${esc(CATEGORIES[k].name)}</b><ul>${cmpPeople().map(q => { const s = q.cats.find(c => c.key === k).score; return `<li><span>${esc(q.name)}</span><span class="num">${pctTxt(s)}</span></li>`; }).join('')}</ul>`;
    }
    pop.innerHTML = html;
    pop.hidden = false;
    const w = 300, left = Math.min(x + 16, innerWidth - w - 12);
    pop.style.left = left + 'px';
    pop.style.top = Math.min(y + 16, innerHeight - pop.offsetHeight - 12) + 'px';
    requestAnimationFrame(() => pop.classList.add('show'));
  }
  function hidePop() { pop.classList.remove('show'); pop.hidden = true; }
  document.addEventListener('pointermove', e => {
    const r = e.target.closest && e.target.closest('.canvas .rung');
    if (r && e.pointerType === 'mouse' && P()) showPop(r, e.clientX, e.clientY); else if (!pop.hidden) hidePop();
  });
  document.addEventListener('scroll', hidePop, true);

  // ---------- views ----------
  function viewLanding(error) {
    return `<div class="view">
      <div class="title-row"><div><h1>Read a raw DNA file</h1>
        <p class="sub">Scores 16 health and trait categories from ${SNP_DB.length} researched markers and checks every SNP against ClinVar. Load two or more files to compare people.</p></div></div>
      <label class="drop" id="drop">
        <input type="file" id="dropFile" accept=".csv,.txt,.tsv" multiple hidden>
        <span id="dropInner"><strong>Drop raw DNA files here</strong><p>.csv or .txt, unzipped. Several files at once open the comparison.</p><span class="btn primary">Choose files</span></span>
      </label>
      ${error ? `<p class="err" role="alert">${esc(error)}</p>` : ''}
      <div class="formats">
        <div class="card"><b>MyHeritage</b><code>MyHeritage_raw_dna_data_*.csv</code></div>
        <div class="card"><b>23andMe</b><code>genome_*_Full_*.txt</code></div>
        <div class="card"><b>AncestryDNA</b><code>AncestryDNA.txt</code></div>
        <div class="card"><b>FamilyTreeDNA</b><code>*_Autosomal_o37_Results_*.csv</code></div>
      </div>
      <section class="card" style="margin-top:var(--gap)" aria-labelledby="h-what"><div class="card-head"><h2 id="h-what">What the report scores</h2><span class="small muted">${SNP_DB.length} markers in ${Object.keys(CATEGORIES).length} categories</span></div>
        <dl class="scored">${Object.values(CATEGORIES).map(c => `<div><dt>${esc(c.name)}</dt><dd>${esc(c.desc.split('. ')[0].replace(/\.$/, ''))}.</dd></div>`).join('')}</dl></section>
      <p class="sub small" style="margin-top:16px">The file is read inside this browser window and never uploaded. ${hasClinvar() ? '' : 'ClinVar data is still loading or missing — run <code>python build_clinvar.py</code> if the gene scan stays empty.'}</p>
    </div>`;
  }

  function titleRow(h1, sub, finding, actions = true, reads = false) {
    return `<div class="title-row"><div><h1>${h1}</h1>${sub ? `<p class="sub">${sub}</p>` : ''}${finding ? `<p class="finding">${finding}</p>` : ''}</div>
      ${actions ? `<div class="actions">${reads ? `<button class="btn" data-act="cram"${state.readJob ? ' disabled' : ''}>${rp().readCheck ? 'Check other reads' : 'Check against reads'}</button>` : ''}<button class="btn" data-act="print">Print or save PDF</button><button class="btn primary" data-act="download">Download report</button></div>` : ''}</div>`;
  }

  // The read check's status line: progress while scanning, the summary when done, or the error.
  function readNotice(p) {
    if (state.readJob && state.readJob.p === p) return `<div class="card readbar" role="status"><div><b>Checking against ${esc(state.readJob.name)}</b><span class="small muted" id="readPct"> 0%</span></div><span class="progress"><i id="readProg" style="--p:0"></i></span></div>`;
    if (p.readError) return `<p class="err" role="alert">${esc(p.readError)}</p>`;
    const rc = p.readCheck;
    if (!rc) return '';
    const t = { seen: 0, missing: 0, conflict: 0, agree: 0, none: 0 };
    for (const [id, ref] of rc.ref) t[readVerdict(p.geno.get(id).g, ref, rc.obs.get(id)).key]++;
    const covered = rc.ref.size - t.none;
    const wrong = rc.n >= 50 && rc.disc / rc.n > 0.03;
    return `<div class="card readbar${wrong ? ' bad' : ''}"><div><b>Read check: ${esc(rc.file)}</b>
      <p class="small muted" style="margin:2px 0 0">${covered} of ${rc.ref.size} checked sites are covered by reads (${fmt(rc.reads)} reads at about 1× coverage). Variants this file carries: ${t.seen} seen in reads, ${t.missing} not seen${t.missing ? ' (at 1–2 reads a single copy is easily missed)' : ''}. ${t.conflict ? `${t.conflict} sites where reads show an allele this file doesn't have. ` : ''}${t.agree} other sites match.</p>
      ${wrong ? `<p class="small" style="margin:6px 0 0"><b>These reads may not be ${esc(p.name)}'s:</b> ${(rc.disc / rc.n * 100).toFixed(1)}% of read bases are not in this file's genotypes (expected under 2% for the same person).</p>` : ''}</div></div>`;
  }
  const lowPass = p => /low-pass/i.test(p.method) ? `<div class="notice" role="note"><span>Low-pass sequencing: most genotypes in this file are imputed, not read directly. Common variants are usually right; rare variants (most ClinVar “pathogenic” calls) are often wrong. Confirm anything important with a clinical test.</span></div>` : '';

  function overviewFinding(p) {
    const scored = p.cats.filter(c => c.score !== null);
    const hi = scored.reduce((a, b) => b.score > a.score ? b : a), lo = scored.reduce((a, b) => b.score < a.score ? b : a);
    const plp = p.plp.length ? ` ${p.plp.length} ClinVar pathogenic variant${p.plp.length > 1 ? 's' : ''} carried${p.plp.every(h => h.k === 1) ? ' (one copy — carrier)' : ''}.` : ' No ClinVar pathogenic variants carried.';
    return `Strongest: ${esc(hi.name)} ${pctTxt(hi.score)}. Weakest: ${esc(lo.name)} ${pctTxt(lo.score)}.${plp}`;
  }

  function glance(p) {
    const g = id => p.found.find(f => f.e.id === id);
    const rows = [];
    if (p.apoe) rows.push(['APOE genotype', `${p.apoe.approx && p.apoe.text !== 'undetermined' ? 'likely ' : ''}${p.apoe.text}`]);
    const lct = g('rs4988235'); if (lct) rows.push(['Lactose', lct.r.k ? 'Likely tolerant' : 'Likely intolerant']);
    const caf = g('rs762551'); if (caf) rows.push(['Caffeine', caf.r.k ? 'Slow metabolizer' : 'Fast metabolizer']);
    const eye = g('rs12913832'); if (eye) rows.push(['Eye colour', ['Brown likely', 'Brown, hazel or green', 'Blue or light likely'][eye.r.k]]);
    const fv = g('rs6025'); if (fv) rows.push(['Factor V Leiden', fv.r.k ? (fv.r.k === 1 ? 'Heterozygous' : 'Homozygous') : 'Not carried']);
    rows.push(['ClinVar pathogenic', p.plp.length ? `${p.plp.length} carried` : 'None found']);
    return `<dl class="kv">${rows.map(([k, v]) => `<dt>${k}</dt><dd>${esc(v)}</dd>`).join('')}</dl>`;
  }

  function attentionList(p, cap = 6) {
    const items = [
      ...p.plp.map(h => ({ gene: h.gene || '—', id: h.id, g: h.g, alt: h.alt, text: h.cond.split(';')[0], pill: `<span class="pill bad">${h.k === 2 ? 'Pathogenic, 2 copies' : 'Pathogenic, carrier'}</span>`, href: '#/genes' })),
      ...p.notable.map(n => ({ gene: n.f.e.gene, id: n.f.e.id, g: n.f.r.g, alt: n.f.r.alt, text: n.f.e.trait, pill: n.w > 0 ? '<span class="pill good">Protective</span>' : '<span class="pill bad">Unfavourable</span>', href: `#/category/${Object.keys(n.f.e.cats)[0]}` })),
    ];
    if (!items.length) return '<p class="empty">No high-impact unfavourable markers or pathogenic variants found.</p>';
    return `<ul class="list">${items.slice(0, cap).map(it => `<li><span><a href="${it.href}" class="link">${esc(it.gene)}</a> ${gtCells(it.g, it.alt)}</span>${it.pill}<span class="meta">${esc(it.text)} · <span class="mono">${esc(it.id)}</span>${(v => v ? ` · reads: ${esc(v.label.toLowerCase())}` : '')(readVerdictFor(p, it.id, it.g))}</span></li>`).join('')}</ul>
      ${items.length > cap ? `<p class="small muted" style="margin:12px 0 0">${items.length - cap} more in Categories and Gene scan.</p>` : ''}`;
  }

  function viewOverview(p) {
    return `<div class="view">
      ${titleRow(esc(p.name), `${esc(p.file)} · ${fmt(p.geno.size)} SNPs${p.method ? ' · ' + esc(p.method) : ''}${p.build ? ' · ' + p.build : ''}`, overviewFinding(p), true, true)}
      ${readNotice(p)}${lowPass(p)}
      <div class="card stats">
        <div class="stat"><span>SNPs read</span><b data-count="${p.geno.size}" data-suffix="">${fmt(p.geno.size)}</b></div>
        <div class="stat"><span>Curated markers found</span><b><span data-count="${p.found.length}" data-suffix="">${p.found.length}</span> <small>of ${SNP_DB.length}</small></b></div>
        <div class="stat"><span>ClinVar variants carried</span><b data-count="${p.hits.length}" data-suffix="">${p.hits.length}</b></div>
        <div class="stat"><span>Genes with known-function variants</span><b data-count="${p.genes.size}" data-suffix="">${p.genes.size}</b></div>
      </div>
      <div class="row2">
        <section class="card" aria-labelledby="h-scores">
          <div class="card-head"><h2 id="h-scores">Category scores</h2><span class="small muted">ticks are markers · dot is the score</span></div>
          ${ladder(p.cats.filter(c => c.score !== null))}
        </section>
        <div class="stack">
          <section class="card" aria-labelledby="h-att"><div class="card-head"><h2 id="h-att">Needs attention</h2><span class="small muted num">${p.plp.length + p.notable.length}</span></div>${attentionList(p)}</section>
          <section class="card" aria-labelledby="h-gl"><div class="card-head"><h2 id="h-gl">At a glance</h2></div>${glance(p)}</section>
        </div>
      </div>
    </div>`;
  }

  function viewCategory(p, key) {
    const c = p.cats.find(c => c.key === key) || p.cats.find(c => c.score !== null);
    const chips = p.cats.filter(x => x.total).map(x => `<a class="chip" href="#/category/${x.key}"${x.key === c.key ? ' aria-current="page"' : ''}>${esc(x.name)} <span class="num">${pctTxt(x.score)}</span></a>`).join('');
    const fav = c.items.filter(i => i.val > 0.3).length, unf = c.items.filter(i => i.val < -0.3).length;
    return `<div class="view">
      ${titleRow(esc(c.name), esc(c.desc), `${c.items.length} of ${c.total} markers found in this file · ${fav} favourable, ${unf} unfavourable.`, true, true)}
      ${readNotice(p)}
      <nav class="cats" aria-label="Categories">${chips}</nav>
      <section class="card" style="margin-bottom:var(--gap)">${c.score === null ? '<p class="empty">None of this category’s markers are in this file.</p>' : ladder([c])}</section>
      <section class="card" aria-labelledby="h-mk"><div class="card-head"><h2 id="h-mk">Markers</h2><span class="small muted">sorted from least to most favourable</span></div>
        <div class="tablewrap">${markerHead()}${c.items.map(i => markerRow(i.f, i.val, i.w)).join('')}</tbody></table></div>
        ${c.missing.length ? `<p class="small muted" style="margin:14px 0 0">Not in this file: ${c.missing.map(f => `${esc(f.e.gene)} <span class="mono">${esc(f.e.id)}</span>${f.r && f.r.mismatch ? ' (allele mismatch)' : ''}`).join(', ')}.</p>` : ''}
      </section>
    </div>`;
  }

  function cvRows(hits) {
    return hits.map(x => `<tr>
      <td class="gene"><b>${esc(x.gene || '—')}</b><span class="id">${x.id}</span></td>
      <td>${gtCells(x.g, x.alt)}<span class="gt-note">${x.ref}&gt;${x.alt} · ${x.k === 2 ? 'two copies' : 'one copy'}</span></td>
      <td class="effect">${esc(x.cond)}</td>
      <td class="num" aria-label="${x.stars} of 4 review stars">${x.stars ? `${x.stars} of 4` : '—'}</td>
      <td><span class="pill ${x.sig === 'P' || x.sig === 'LP' ? 'bad' : x.sig === 'RF' ? 'warn' : x.sig === 'PR' ? 'good' : ''}">${SIG[x.sig]}</span></td>${readsTd(x.id, x.g)}
      <td class="src"><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/${x.vid}/" target="_blank" rel="noopener">ClinVar ${x.vid}</a>${links(null, x.id)}</td></tr>`).join('');
  }
  const cvHead = () => `<table class="rows"><thead><tr><th>Gene</th><th>Genotype</th><th>Condition or drug</th><th class="num">Review</th><th>Classification</th>${readsTh()}<th>Sources</th></tr></thead><tbody>`;
  function cvFiltered(p) {
    const q = state.cvQuery.toLowerCase();
    return p.hits.filter(h => (state.cvSig === 'all' || h.sig === state.cvSig) && h.stars >= state.cvStars && (!q || (h.gene + ' ' + h.cond + ' ' + h.id).toLowerCase().includes(q)));
  }
  function viewGenes(p) {
    if (!hasClinvar()) return `<div class="view">${titleRow('Gene scan', '', '', false)}<p class="empty">ClinVar data isn’t loaded. Run <code>python build_clinvar.py</code> next to index.html, then reload.</p></div>`;
    const segs = [['all', 'All'], ...Object.entries(SIG)].map(([k, v]) => `<button data-sig="${k}" aria-pressed="${state.cvSig === k}">${v} <span class="num muted">${k === 'all' ? p.hits.length : p.hits.filter(h => h.sig === k).length}</span></button>`).join('');
    return `<div class="view">
      ${titleRow('Gene scan', `Every SNP in ${esc(p.name)}’s file checked against NCBI ClinVar (${esc(CLINVAR_DATE)}, ${fmt(Object.keys(CLINVAR).length)} classified variants). Listed: variants where the classified allele is carried.`, p.plp.length ? `${p.plp.length} pathogenic or likely pathogenic, ${p.plp.filter(h => h.k === 1).length} of them as a single copy. For recessive conditions one copy means carrier status only.` : 'No pathogenic or likely pathogenic variants carried.', true, true)}
      ${readNotice(p)}${lowPass(p)}
      <section class="card" aria-labelledby="h-cv">
        <div class="card-head"><h2 id="h-cv">Variants carried</h2><span class="small muted num" id="cvCount"></span></div>
        <div class="toolbar"><div class="seg" role="group" aria-label="Classification">${segs}</div>
          <select class="search" id="cvStars" aria-label="Minimum review stars">${[0, 1, 2, 3].map(n => `<option value="${n}"${state.cvStars === n ? ' selected' : ''}>${n ? `${n}+ review stars` : 'Any review status'}</option>`).join('')}</select>
          <input class="search" id="cvQuery" type="search" placeholder="Search gene or condition" aria-label="Search gene or condition" value="${esc(state.cvQuery)}"></div>
        <div class="tablewrap" id="cvTable"></div>
      </section></div>`;
  }
  function drawCv() {
    const hits = cvFiltered(P()), shown = hits.slice(0, 400);
    $('#cvTable').innerHTML = hits.length ? `${cvHead()}${cvRows(shown)}</tbody></table>` : '<p class="empty">No carried variants match these filters.</p>';
    $('#cvCount').textContent = `${shown.length} of ${P().hits.length}`;
  }

  function viewTraits(p) {
    const traits = p.found.filter(f => !Object.keys(f.e.cats).length);
    return `<div class="view">${titleRow('Traits', 'Characteristics with no better or worse genotype. Not scored.', '')}
      <section class="card"><div class="tablewrap">${traits.length ? markerHead().replace('<th>Impact</th>', '<th></th>') + traits.map(f => markerRow(f)).join('') + '</tbody></table>' : '<p class="empty">None of the trait markers are in this file.</p>'}</div></section></div>`;
  }

  function cmpPeople() { return state.people.filter((_, i) => !state.cmpSel || state.cmpSel.has(i)); }
  function viewCompare() {
    if (state.people.length < 2) return `<div class="view">${titleRow('Compare', 'Load a second raw DNA file to see two or more people side by side.', '', false)}
      <label class="drop" id="drop"><input type="file" id="dropFile" accept=".csv,.txt,.tsv" multiple hidden><span id="dropInner"><strong>Add another person</strong><p>Drop a raw DNA file here or choose one.</p><span class="btn primary">Choose files</span></span></label></div>`;
    const ps = cmpPeople();
    const pairs = [];
    for (let i = 0; i < ps.length; i++) for (let j = i + 1; j < ps.length; j++) pairs.push([ps[i], ps[j]]);
    const rel = pairs.map(([a, b]) => ({ a, b, c: (a.ibsCache ??= {})[b.file + b.name] ??= compareGenomes(a.geno, b.geno) }));
    const profileSim = pairs.map(([a, b]) => {
      const d = a.cats.map((c, i) => c.score !== null && b.cats[i].score !== null ? Math.abs(c.score - b.cats[i].score) : null).filter(x => x !== null);
      return 100 - d.reduce((s, x) => s + x, 0) / d.length;
    });
    const avg = xs => xs.reduce((s, x) => s + x, 0) / xs.length;
    const diffMarkers = SNP_DB.filter(e => { const g = ps.map(p => lookup(e, p.geno)); return g.every(r => r && !r.mismatch) && new Set(g.map(r => r.k)).size > 1; });
    const cvKey = h => h.id + h.alt + h.sig;
    const cvAll = new Map();
    ps.forEach((p, i) => p.hits.forEach(h => { if (!cvAll.has(cvKey(h))) cvAll.set(cvKey(h), { h, who: new Set() }); cvAll.get(cvKey(h)).who.add(i); }));
    const cvUnique = [...cvAll.values()].filter(v => v.who.size < ps.length).sort((x, y) => ['P', 'LP', 'RF', 'DR', 'PR'].indexOf(x.h.sig) - ['P', 'LP', 'RF', 'DR', 'PR'].indexOf(y.h.sig) || y.h.stars - x.h.stars);
    const sel = state.people.map((p, i) => `<button data-pick="${i}" aria-pressed="${!state.cmpSel || state.cmpSel.has(i)}">${esc(p.name)}</button>`).join('');
    const many = pairs.length > 1 ? ' (average of pairs)' : '';

    return `<div class="view">
      ${titleRow('Compare', ps.map(p => esc(p.name)).join(', '), `Genotypes match at ${(avg(rel.map(r => r.c.concordance)) * 100).toFixed(1)}% of shared SNPs${many}; health-profile similarity ${Math.round(avg(profileSim))}%. ${diffMarkers.length} of ${SNP_DB.length} curated markers differ.`)}
      ${state.people.length > 2 ? `<div class="toolbar"><div class="seg" role="group" aria-label="People in comparison">${sel}</div></div>` : ''}
      <div class="card stats">
        <div class="stat"><span>Genotype match${many}</span><b><span data-count="${(avg(rel.map(r => r.c.concordance)) * 100).toFixed(1)}" data-suffix="%" data-dec="1">${(avg(rel.map(r => r.c.concordance)) * 100).toFixed(1)}%</span></b></div>
        <div class="stat"><span>Health-profile similarity</span><b data-count="${Math.round(avg(profileSim))}">${Math.round(avg(profileSim))}%</b></div>
        <div class="stat"><span>Curated markers that differ</span><b><span data-count="${diffMarkers.length}" data-suffix="">${diffMarkers.length}</span> <small>of ${SNP_DB.length}</small></b></div>
        <div class="stat"><span>ClinVar variants not shared</span><b data-count="${cvUnique.length}" data-suffix="">${cvUnique.length}</b></div>
      </div>
      <div class="row2">
        <section class="card" aria-labelledby="h-cl"><div class="card-head"><h2 id="h-cl">Category scores</h2><span class="legend small muted">${ps.map(p => `<span>${initial(p)}${esc(p.name)}</span>`).join('')}</span></div>${compareLadder(ps)}</section>
        <div class="stack">
          <section class="card" aria-labelledby="h-rel"><div class="card-head"><h2 id="h-rel">Relatedness</h2><span class="small muted">identity by state</span></div>
            <ul class="list">${rel.map(r => `<li><span>${esc(r.a.name)} and ${esc(r.b.name)}</span><span class="num">${(r.c.concordance * 100).toFixed(1)}%</span>
              <span class="meta">${esc(relationReading(r.c))} · kinship ${r.c.kinship.toFixed(3)} · ${(r.c.ibs0 * 100).toFixed(2)}% opposite homozygotes · ${fmt(r.c.n)} autosomal SNPs</span></li>`).join('')}</ul>
            <p class="small muted" style="margin:12px 0 0">Estimated from autosomal SNPs. Reliable for first- and second-degree relatives; more distant links need a genealogy service's segment matching.</p></section>
          <section class="card" aria-labelledby="h-lead"><div class="card-head"><h2 id="h-lead">Who leads where</h2></div>
            <dl class="kv">${ps.map(p => { const lead = Object.keys(CATEGORIES).filter(k => { const s = ps.map(q => q.cats.find(c => c.key === k).score); return s.every(x => x !== null) && p.cats.find(c => c.key === k).score === Math.max(...s) && new Set(s).size > 1; }); return `<dt class="who">${initial(p)}${esc(p.name)}</dt><dd class="num">${lead.length} categories</dd>`; }).join('')}</dl></section>
        </div>
      </div>
      <section class="card" style="margin-top:var(--gap)" aria-labelledby="h-md">
        <div class="card-head"><h2 id="h-md">Markers</h2><span class="small muted num" id="mkCount"></span></div>
        <div class="toolbar"><div class="seg" role="group" aria-label="Marker filter">
          <button data-cmp-mode="diff" aria-pressed="${state.cmpMode === 'diff'}">Different</button>
          <button data-cmp-mode="risk" aria-pressed="${state.cmpMode === 'risk'}">Shared unfavourable</button>
          <button data-cmp-mode="all" aria-pressed="${state.cmpMode === 'all'}">All</button></div></div>
        <div class="tablewrap" id="mkTable"></div>
      </section>
      <section class="card" style="margin-top:var(--gap)" aria-labelledby="h-cvd">
        <div class="card-head"><h2 id="h-cvd">ClinVar variants not shared</h2><span class="small muted num">${cvUnique.length}</span></div>
        ${cvUnique.length ? `<div class="tablewrap"><table class="rows"><thead><tr><th>Gene</th><th>Condition or drug</th><th>Classification</th>${ps.map(p => `<th>${esc(p.name)}</th>`).join('')}</tr></thead><tbody>
          ${cvUnique.slice(0, 200).map(({ h }) => `<tr><td class="gene"><b>${esc(h.gene || '—')}</b><span class="id">${h.id}</span></td><td class="effect">${esc(h.cond)}</td>
            <td><span class="pill ${h.sig === 'P' || h.sig === 'LP' ? 'bad' : h.sig === 'RF' ? 'warn' : h.sig === 'PR' ? 'good' : ''}">${SIG[h.sig]}</span></td>
            ${ps.map(p => { const r = p.geno.get(h.id); return `<td>${r ? gtCells(r.g, h.alt) : '<span class="muted">—</span>'}</td>`; }).join('')}</tr>`).join('')}
        </tbody></table></div>` : '<p class="empty">Everyone carries the same ClinVar-classified variants.</p>'}
      </section></div>`;
  }
  function markersTable() {
    const ps = cmpPeople();
    const rows = [];
    for (const e of SNP_DB) {
      const rs = ps.map(p => lookup(e, p.geno));
      if (!rs.every(r => r && !r.mismatch)) continue;
      const differs = new Set(rs.map(r => r.k)).size > 1;
      const w = Object.values(e.cats)[0];
      const allBad = w !== undefined && rs.every(r => markerVal(e, w, r.k) < -0.3);
      if (state.cmpMode === 'diff' && !differs) continue;
      if (state.cmpMode === 'risk' && !allBad) continue;
      const word = r => { const v = markerVal(e, w, r.k); return v > 0.3 ? 'favourable' : v < -0.3 ? 'unfavourable' : 'intermediate'; };
      rows.push(`<tr><td class="gene"><b>${esc(e.gene)}</b><span class="id">${esc(e.id)}</span></td><td class="effect"><span class="t">${esc(e.trait)}</span>${w === undefined ? 'Trait' : esc(Object.keys(e.cats).map(k => CATEGORIES[k].name).join(', '))}</td>
        ${rs.map((r, i) => `<td>${gtCells(r.g, r.alt, i > 0 && r.k !== rs[0].k)}${w === undefined ? '' : `<span class="gt-note">${word(r)}</span>`}</td>`).join('')}</tr>`);
    }
    const html = rows.length ? `<table class="rows"><thead><tr><th>Gene</th><th>Trait and category</th>${ps.map(p => `<th>${esc(p.name)}</th>`).join('')}</tr></thead><tbody>${rows.join('')}</tbody></table>` : '<p class="empty">No markers match this filter.</p>';
    return { html, n: rows.length };
  }
  function drawMarkers() {
    const t = markersTable();
    $('#mkTable').innerHTML = t.html;
    $('#mkCount').textContent = `${t.n} markers`;
  }

  function viewMethod() {
    return `<div class="view">${titleRow('Method and sources', '', '', false)}
      <section class="card method">
        <p><b>Markers.</b> ${SNP_DB.length} variants from peer-reviewed genome-wide association studies, the NHGRI-EBI GWAS Catalog, CPIC pharmacogenomic guidelines and PharmGKB, and ClinVar. Each row links to its primary paper and to dbSNP, ClinVar and the GWAS Catalog. Effect alleles are on the GRCh37 forward strand and were cross-checked against the GWAS Catalog and ClinVar. When a marker is missing, a tag SNP in linkage disequilibrium stands in and is labelled.</p>
        <p><b>Scores.</b> Each marker scores +1 for its best genotype and −1 for its worst; heterozygotes fall in between according to the variant’s model (additive, dominant or recessive). Markers are weighted by effect size (1–3) times evidence (strong 1, moderate 0.75, limited 0.5). Category score = 50 + 50 × the weighted mean, so 100% is the best combination possible for the markers tested, 50% an even balance and 0% the worst possible. A typical person can score low where most people lack a rare protective variant.</p>
        <p><b>Gene scan.</b> Every SNP is checked against NCBI ClinVar${hasClinvar() ? ` (${esc(CLINVAR_DATE)})` : ''}: pathogenic and likely pathogenic variants with at least one review star, plus risk factors, protective variants and drug responses. Refresh with <code>python build_clinvar.py --refresh</code>.</p>
        <p><b>Comparison.</b> Genotype match is the share of SNPs present in both files where both people have the same genotype. Relatedness uses autosomal SNPs only and the KING-robust kinship estimator (Manichaikul et al. 2010): about 0.25 for parent–child and full siblings, 0.125 for second-degree relatives; opposite homozygotes (e.g. AA against GG), near zero only between parent and child, separate those two. Health-profile similarity is 100 minus the mean gap between category scores.</p>
        <p><b>Read check.</b> A MyHeritage CRAM file holds the sequencing reads the genotype file was imputed from. “Check against reads” scans it in the browser and counts, at each ClinVar hit and curated marker, which bases the reads carry (mapping quality and base quality 20 or more). A variant “seen in reads” is confirmed directly; at about 1× coverage most sites have 0–2 reads, so “not seen” in one or two reads is weak evidence against a single copy. The CRAM must be aligned to GRCh37. If the reads disagree with the genotype file at more than 3% of bases, the check warns that they may belong to someone else.</p>
        <p><b>Limits.</b> Most traits depend on hundreds of variants plus lifestyle and environment; these scores cover a few known variants with mostly small effects (odds ratios 1.1–1.5), measured mainly in people of European ancestry. Consumer arrays and low-pass sequencing miss most rare mutations (most BRCA1/2 variants), copy-number changes and CYP2D6 duplications. This is education and research, not diagnosis — talk to a doctor or genetic counsellor before acting on a result.</p>
      </section></div>`;
  }

  // ---------- exported / printed report ----------
  function reportHTML(p) {
    RP = p;
    try { return reportBody(p); } finally { RP = null; }
  }
  function reportBody(p) {
    const traits = p.found.filter(f => !Object.keys(f.e.cats).length);
    return `<div class="report">
      ${titleRow(`Genetic report — ${esc(p.name)}`, `${esc(p.file)} · ${fmt(p.geno.size)} SNPs · generated ${new Date().toLocaleString()} by Allele ${VERSION}`, overviewFinding(p), false)}
      ${readNotice(p)}${lowPass(p)}
      <section><h2>Category scores</h2><div class="card">${ladder(p.cats.filter(c => c.score !== null), { static: true })}</div></section>
      <section><h2>Needs attention</h2><div class="card">${attentionList(p, 50)}</div></section>
      <section><h2>At a glance</h2><div class="card">${glance(p)}</div></section>
      ${p.cats.filter(c => c.items.length).map(c => `<section><h2>${esc(c.name)} — ${pctTxt(c.score)}</h2><p class="sub small">${esc(c.desc)}</p><div class="card">${markerHead()}${c.items.map(i => markerRow(i.f, i.val, i.w)).join('')}</tbody></table></div></section>`).join('')}
      ${traits.length ? `<section><h2>Traits</h2><div class="card">${markerHead()}${traits.map(f => markerRow(f)).join('')}</tbody></table></div></section>` : ''}
      ${p.hits.length ? `<section><h2>ClinVar gene scan</h2><div class="card">${cvHead()}${cvRows(p.hits)}</tbody></table></div></section>` : ''}
      <section>${viewMethod().replace(/<div class="view">|<\/div>$/g, '')}</section></div>`;
  }
  function staticize(html) {
    const d = document.createElement('div');
    d.innerHTML = html;
    d.querySelectorAll('.slide').forEach(s => s.style.setProperty('--x', s.style.getPropertyValue('--to')));
    d.querySelectorAll('.ladder').forEach(l => l.classList.add('on'));
    d.querySelectorAll('.view').forEach(v => v.classList.remove('view'));
    return d.innerHTML;
  }
  function exportHTML() {
    const d = document.createElement('div');
    if (isCompare()) {
      d.innerHTML = `<div class="report">${viewCompare()}</div>`;
      const save = state.cmpMode; state.cmpMode = 'diff';
      const t = markersTable(); state.cmpMode = save;
      d.querySelector('#mkTable').innerHTML = t.html;
      d.querySelector('#mkCount').textContent = `${t.n} markers that differ`;
      d.querySelectorAll('.toolbar').forEach(n => n.remove());
    } else d.innerHTML = reportHTML(P());
    d.querySelectorAll('.actions').forEach(n => n.remove());
    return staticize(d.innerHTML);
  }
  const isCompare = () => location.hash.startsWith('#/compare') && cmpPeople().length > 1;
  function doDownload() {
    const name = isCompare() ? 'DNA comparison ' + cmpPeople().map(p => p.name).join(' ') : 'DNA report ' + P().name;
    const html = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>${esc(name)}</title><style>${document.getElementById('css').textContent}</style></head><body>${exportHTML()}</body></html>`;
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([html], { type: 'text/html' }));
    a.download = name.replace(/\s+/g, '_') + '.html';
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 5000);
  }
  function doPrint() {
    $('#print').innerHTML = exportHTML();
    window.print();
  }

  // ---------- router ----------
  const main = $('#main');
  function render(focus) {
    hidePop();
    const [, view = 'overview', arg] = location.hash.split('/');
    document.querySelectorAll('#nav a').forEach(a => {
      a.toggleAttribute('aria-current', a.dataset.view === view && !!P());
      if (a.dataset.view === view && P()) a.setAttribute('aria-current', 'page');
      a.classList.toggle('off', !P());
    });
    drawPeople();
    if (!P()) { main.innerHTML = viewLanding(state.error); wireDrop(); return; }
    const p = P();
    main.innerHTML = view === 'category' ? viewCategory(p, arg)
      : view === 'genes' ? viewGenes(p)
      : view === 'traits' ? viewTraits(p)
      : view === 'compare' ? viewCompare()
      : view === 'method' ? viewMethod()
      : viewOverview(p);
    if (view === 'genes' && hasClinvar()) drawCv();
    if (view === 'compare' && state.people.length > 1) drawMarkers();
    wireDrop();
    animate(main);
    if (focus) main.focus({ preventScroll: true });
    window.scrollTo({ top: 0 });
  }
  addEventListener('hashchange', () => { $('#side').classList.remove('open'); $('#menuBtn').setAttribute('aria-expanded', 'false'); render(true); });

  main.addEventListener('click', e => {
    const t = e.target.closest('button');
    if (!t) return;
    if (t.dataset.act === 'download') doDownload();
    else if (t.dataset.act === 'print') doPrint();
    else if (t.dataset.act === 'cram') $('#cramFile').click();
    else if (t.dataset.sig) { state.cvSig = t.dataset.sig; main.querySelectorAll('[data-sig]').forEach(b => b.setAttribute('aria-pressed', b === t)); drawCv(); }
    else if (t.dataset.cmpMode) { state.cmpMode = t.dataset.cmpMode; main.querySelectorAll('[data-cmp-mode]').forEach(b => b.setAttribute('aria-pressed', b === t)); drawMarkers(); }
    else if (t.dataset.pick) {
      const i = +t.dataset.pick, s = state.cmpSel || new Set(state.people.map((_, j) => j));
      s.has(i) ? s.delete(i) : s.add(i);
      if (s.size >= 2) { state.cmpSel = s.size === state.people.length ? null : s; render(); }
    }
  });
  main.addEventListener('input', e => {
    if (e.target.id === 'cvQuery') { state.cvQuery = e.target.value; drawCv(); }
    if (e.target.id === 'cvStars') { state.cvStars = +e.target.value; drawCv(); }
  });

  // ---------- people ----------
  function drawPeople() {
    $('#people').innerHTML = state.people.map((p, i) => `<li><button class="person" data-i="${i}" aria-current="${i === state.cur}">${initial(p)}<span>${esc(p.name)}</span></button><button class="remove" data-rm="${i}" aria-label="Remove ${esc(p.name)}">×</button></li>`).join('');
  }
  $('#people').addEventListener('click', e => {
    const b = e.target.closest('button');
    if (!b) return;
    if (b.dataset.rm !== undefined) {
      state.people.splice(+b.dataset.rm, 1);
      state.cur = Math.min(state.cur, Math.max(0, state.people.length - 1));
      state.cmpSel = null;
    } else state.cur = +b.dataset.i;
    if (location.hash.startsWith('#/compare') && state.people.length < 2) location.hash = '#/overview';
    render();
  });
  $('#menuBtn').addEventListener('click', () => { const o = $('#side').classList.toggle('open'); $('#menuBtn').setAttribute('aria-expanded', o); });

  // ---------- read check ----------
  let cramLib;
  const loadCramLib = () => cramLib ??= new Promise((ok, fail) => {
    const s = document.createElement('script');
    s.src = 'vendor/cram.js';
    s.onload = () => ok(window.CramLib);
    s.onerror = () => fail(new Error('Could not load vendor/cram.js, the CRAM reader.'));
    document.head.appendChild(s);
  });
  async function checkReads(file) {
    const p = P();
    if (!p || state.readJob) return;
    p.readError = '';
    state.readJob = { p, name: file.name };
    render();
    try {
      const lib = await loadCramLib();
      const sites = readSites(p);
      let last = 0;
      const res = await runReadCheck(file, sites, lib.CramFile, lib.BlobFile, k => {
        const now = performance.now();
        if (now - last < 100 && k < 1) return;
        last = now;
        const bar = $('#readProg'), pct = $('#readPct');
        if (bar) bar.style.setProperty('--p', k);
        if (pct) pct.textContent = ` ${Math.round(k * 100)}%`;
      });
      let n = 0, disc = 0;
      for (const x of sites) for (const [b, c] of Object.entries(res.obs.get(x.id))) { n += c; if (!p.geno.get(x.id).g.includes(b)) disc += c; }
      p.readCheck = { file: file.name, obs: res.obs, ref: new Map(sites.map(x => [x.id, x.ref])), n, disc, reads: res.reads };
    } catch (err) {
      p.readError = `${file.name}: ${err.message}`;
    }
    state.readJob = null;
    render();
  }
  $('#cramFile').onchange = e => { if (e.target.files[0]) checkReads(e.target.files[0]); e.target.value = ''; };

  async function readFiles(files) {
    files = [...files].filter(f => f.size);
    const cram = files.find(f => /\.cram$/i.test(f.name));
    if (cram) {
      if (P()) checkReads(cram); else { state.error = 'Load the person’s raw DNA file (.csv or .txt) first, then check it against the CRAM reads.'; render(); }
      files = files.filter(f => f !== cram);
    }
    if (!files.length) return;
    state.error = '';
    const before = state.people.length;
    for (const file of files) {
      setProgress(0, `Reading ${file.name}`);
      const text = await file.text();
      const lines = text.split(/\r?\n/), geno = new Map(), meta = [];
      for (let i = 0; i < lines.length; i += 40000) {
        parseLines(lines, i, Math.min(lines.length, i + 40000), geno, meta);
        setProgress(i / lines.length, `Reading ${file.name}`);
        await new Promise(r => setTimeout(r));
      }
      if (geno.size < 1000) { state.error = `${file.name}: no genotype rows found. Expected columns rsid, chromosome, position, genotype.`; continue; }
      let name = personName(file.name);
      if (state.people.some(p => p.name === name)) name += ` (${state.people.length + 1})`;
      setProgress(1, `Scoring ${name}`);
      await new Promise(r => setTimeout(r));
      state.people.push(buildPerson(name, file.name, { geno, meta }, HUES[state.people.length % HUES.length]));
    }
    if (state.people.length > before) {
      state.cur = before;
      state.cmpSel = null;
      const target = state.people.length > 1 && (files.length > 1 || location.hash.startsWith('#/compare')) ? '#/compare' : '#/overview';
      if (location.hash === target) render(true); else location.hash = target;
    } else render();
  }
  function setProgress(k, label) {
    const inner = $('#dropInner');
    if (!inner) return;
    inner.innerHTML = `<strong>${esc(label)}…</strong><p class="num">${Math.round(k * 100)}%</p><span class="progress"><i style="--p:${k}"></i></span>`;
  }
  function wireDrop() {
    const drop = $('#drop'), input = $('#dropFile');
    if (!drop) return;
    input.onchange = () => readFiles(input.files);
    drop.ondragover = e => { e.preventDefault(); drop.classList.add('over'); };
    drop.ondragleave = () => drop.classList.remove('over');
    drop.ondrop = e => { e.preventDefault(); e.stopPropagation(); drop.classList.remove('over'); readFiles(e.dataTransfer.files); };
  }
  $('#addFile').onchange = e => { readFiles(e.target.files); e.target.value = ''; };
  addEventListener('dragover', e => e.preventDefault());
  addEventListener('drop', e => { e.preventDefault(); if (e.dataTransfer.files.length) readFiles(e.dataTransfer.files); });

  // clinvar.js is deferred; it has run by DOMContentLoaded (or failed to load, and the gene scan says so)
  addEventListener('DOMContentLoaded', () => { state.people.forEach(refreshClinvar); render(); });
  $('#ver').textContent = `Allele ${VERSION}`;
  window.allele = { readFiles, state }; // handy for debugging from the console
  render();
})();
