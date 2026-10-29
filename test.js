// Self-check: node test.js  — validates db.js and the scoring, parsing and comparison code in app.js.
const fs = require('fs'), vm = require('vm'), assert = require('assert');
const ctx = vm.createContext({});
vm.runInContext(fs.readFileSync(__dirname + '/db.js', 'utf8'), ctx);
vm.runInContext(fs.readFileSync(__dirname + '/pgx.js', 'utf8'), ctx);
vm.runInContext(fs.readFileSync(__dirname + '/app.js', 'utf8'), ctx);
const { SNP_DB, CATEGORIES } = vm.runInContext('({ SNP_DB, CATEGORIES })', ctx);

// Database integrity
const ids = new Set();
for (const e of SNP_DB) {
  assert(!ids.has(e.id), 'duplicate ' + e.id); ids.add(e.id);
  for (const s of [e, ...(e.proxies || [])]) assert(/^[ACGT]$/.test(s.alt) && /^[ACGT]$/.test(s.oth) && s.alt !== s.oth, 'alleles ' + s.id);
  for (const c in e.cats) assert(CATEGORIES[c], 'unknown category ' + c + ' in ' + e.id);
  assert(['strong', 'moderate', 'limited'].includes(e.ev), 'evidence ' + e.id);
  assert(e.fx ? e.fx.length === 3 : e.eff, 'effect text ' + e.id);
  assert(e.ref, 'reference ' + e.id);
}

// Best possible genome scores 100%, worst scores 0% (categories whose markers all point one way)
const mixed = new Set(SNP_DB.filter(e => new Set(Object.values(e.cats).map(Math.sign)).size > 1).flatMap(e => Object.keys(e.cats)));
const genome = best => new Map(SNP_DB.map(e => {
  const s = Object.values(e.cats)[0] || -1, kGood = s < 0 ? 0 : 2;
  const k = best ? kGood : 2 - kGood;
  return [e.id, { g: e.alt.repeat(k) + e.oth.repeat(2 - k) }];
}));
for (const [best, want] of [[true, 100], [false, 0]]) {
  for (const c of ctx.analyze(genome(best)).cats) if (!mixed.has(c.key)) assert.strictEqual(Math.round(c.score), want, `${c.key} ${best ? 'best' : 'worst'}`);
}

// Parsing: MyHeritage CSV, 23andMe TSV, AncestryDNA split alleles, no-calls
const { geno } = ctx.parseRaw('#c\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","5","AG"\nrs2\t1\t6\tCC\nrs3\t2\t7\tT\tA\nrs4\t1\t8\t--\n'.repeat(1));
assert.deepStrictEqual([...geno.keys()], ['rs1', 'rs2', 'rs3']);
assert.strictEqual(geno.get('rs3').g, 'TA');

// APOE haplotype
assert.strictEqual(ctx.apoeType(new Map([['rs429358', { g: 'TC' }], ['rs7412', { g: 'CC' }]])).text, 'ε3/ε4');
assert.strictEqual(ctx.apoeType(new Map([['rs429358', { g: 'TT' }], ['rs7412', { g: 'CT' }]])).text, 'ε2/ε3');

// Comparison: identical genomes match fully; opposite homozygotes count as IBS0
const g1 = new Map([['rs1', { g: 'AA' }], ['rs2', { g: 'CT' }], ['rs3', { g: 'GG' }], ['rs4', { g: 'AG' }]]);
const g2 = new Map([['rs1', { g: 'GG' }], ['rs2', { g: 'TC' }], ['rs3', { g: 'GA' }], ['rs5', { g: 'AA' }]]);
assert.strictEqual(ctx.compareGenomes(g1, g1).concordance, 1);
const c = ctx.compareGenomes(g1, g2);
assert.deepStrictEqual([c.n, c.concordance, c.ibs0, c.ibs], [3, 1 / 3, 1 / 3, 0.5]);
// Sex chromosomes are ignored (father and son differ on X by design)
assert.strictEqual(ctx.compareGenomes(new Map([['rs9', { g: 'AA', chr: 'X' }], ...g1]), new Map([['rs9', { g: 'GG', chr: 'X' }], ...g2])).n, 3);
assert.strictEqual(ctx.personName('MyHeritage_raw_dna_data_Sample.csv'), 'Sample');
assert.strictEqual(ctx.personName('genome_Anna_Berzina_v5_Full_20240101.txt'), 'Anna Berzina');

// Read check: base at a position from CRAM read features (0-based refPos), and the verdict on a genotype
const feats = [{ code: 'S', refPos: 99, data: 'AA' }, { code: 'X', refPos: 104, sub: 'G' }, { code: 'D', refPos: 109, data: 2 }];
assert.strictEqual(JSON.stringify(ctx.baseAt(feats, 99, 105, 'A')), JSON.stringify({ base: 'G', index: 7 }));
assert.strictEqual(JSON.stringify(ctx.baseAt(feats, 99, 103, 'C')), JSON.stringify({ base: 'C', index: 5 }));
assert.strictEqual(ctx.baseAt(feats, 99, 110, 'C'), null);
assert.strictEqual(JSON.stringify(ctx.baseAt(feats, 99, 113, 'T')), JSON.stringify({ base: 'T', index: 13 }));
assert.strictEqual(ctx.readVerdict('AG', 'G', { A: 1, G: 2 }).key, 'seen');
assert.strictEqual(ctx.readVerdict('AG', 'G', { G: 2 }).key, 'missing');
assert.strictEqual(ctx.readVerdict('GG', 'G', { A: 1 }).key, 'conflict');
assert.strictEqual(ctx.readVerdict('GG', 'G', {}).key, 'none');

// ---------- 1.1.0 engines ----------
const J = x => JSON.stringify(x);
const G = obj => new Map(Object.entries(obj).map(([k, g]) => [k, { g, chr: '1', pos: 1 }]));
// dosage on either strand, null when alleles don't fit
assert.strictEqual(ctx.dosage('AG', 'A', 'G'), 1);
assert.strictEqual(ctx.dosage('TT', 'A', 'G'), 2);           // complement strand
assert.strictEqual(ctx.dosage('AC', 'A', 'G'), null);  // fits neither strand
assert(Math.abs(ctx.normCdf(1.96) - 0.975) < 1e-4 && Math.abs(ctx.normCdf(0) - 0.5) < 1e-6);

// polygenic score: carrying more effect alleles than average gives a percentile above 50
const score = { total: 3, variants: [['rs1', 'A', 'G', 0.2, 0.5], ['rs2', 'C', 'T', 0.1, 0.3], ['rs3', 'G', 'A', 0.3, 0.2]] };
const hi = ctx.scorePGS(G({ rs1: 'AA', rs2: 'CC', rs3: 'GG' }), score), avgish = ctx.scorePGS(G({ rs1: 'AG' }), score);
assert(hi.pct > 95 && hi.used === 3, 'pgs high');
assert(Math.abs(avgish.pct - 50) < 1e-3 && avgish.used === 1, 'pgs average');

// admixture: genotypes drawn from population 0's frequencies are assigned mostly to it
const freqs = {}, sim = new Map();
let seed = 7; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
for (let i = 0; i < 3000; i++) {
  const f0 = rnd(), f1 = rnd();
  freqs['rs' + i] = ['A', 'G', f0, f1];
  sim.set('rs' + i, { g: (rnd() < f0 ? 'A' : 'G') + (rnd() < f0 ? 'A' : 'G') });
}
const adm = ctx.admixture(sim, freqs, ['P0', 'P1'], 100);
assert(adm.props[0][0] === 'P0' && adm.props[0][1] > 0.9, 'admixture ' + J(adm.props));

// pharmacogenomics calls
const PGX_GENES = vm.runInContext('PGX_GENES', ctx);
const gene = n => PGX_GENES.find(g => g.gene === n);
assert.strictEqual(ctx.pgxCall(gene('CYP2C19'), G({ rs4244285: 'GA', rs4986893: 'GG', rs12248560: 'CT' })).phenotype, 'Intermediate metabolizer');
assert.strictEqual(ctx.pgxCall(gene('CYP2C19'), G({ rs4244285: 'GG', rs4986893: 'GG', rs12248560: 'TT' })).phenotype, 'Ultrarapid metabolizer');
const c9 = ctx.pgxCall(gene('CYP2C9'), G({ rs1799853: 'CT', rs1057910: 'AC' }));
assert.strictEqual(c9.activity, 0.5); assert.strictEqual(c9.phenotype, 'Poor metabolizer');
assert.strictEqual(ctx.pgxCall(gene('CYP2D6'), G({ rs3892097: 'CT', rs1065852: 'GA', rs28371725: 'CC' })).phenotype, 'Intermediate metabolizer'); // *1/*4, not *10
assert.strictEqual(ctx.pgxCall(gene('TPMT'), G({ rs1800462: 'CC', rs1800460: 'CT', rs1142345: 'TC' })).phenotype, 'Intermediate metabolizer'); // *1/*3A
assert.strictEqual(ctx.pgxCall(gene('CYP3A5'), G({ rs776746: 'CC', rs10264272: 'CC' })).phenotype, 'Poor metabolizer (non-expresser)');
assert.strictEqual(ctx.pgxCall(gene('G6PD'), G({ rs1050828: 'TT' }), 'male').phenotype, 'Deficient');
assert.strictEqual(ctx.pgxCall(gene('DPYD'), G({})).status, 'Not tested');
const prep = ctx.pgxReport(G({ rs4244285: 'AA', rs4986893: 'GG', rs12248560: 'CC' }), null, null);
assert.strictEqual(prep.drugs.find(d => d.drug === 'Clopidogrel').action, 'avoid');

// shared segments: identical chromosome 1 -> one long segment; parent-like (one allele shared) -> half only
const mk = (fn) => { const a = new Map(), b = new Map(); for (let i = 0; i < 3000; i++) { const [x, y] = fn(i); a.set('s' + i, { g: x, chr: '1', pos: 1000 + i * 80000 }); b.set('s' + i, { g: y, chr: '1', pos: 1000 + i * 80000 }); } return [a, b]; };
const [p1, c1] = mk(i => i % 3 ? ['AG', 'AA'] : ['GG', 'AG']);
const seg = ctx.sharedSegments(p1, c1);
assert(seg.half.length === 1 && seg.total > 200 && seg.full2 === 0, 'parent-child segments ' + J({ n: seg.half.length, t: seg.total, f: seg.full2 }));
const [u1, u2] = mk(i => i % 7 === 0 ? ['AA', 'GG'] : ['AG', 'AG']);
assert.strictEqual(ctx.sharedSegments(u1, u2).half.length, 0);
// trio origin
assert.strictEqual(J(ctx.originAt('AG', 'AA', 'GG')), J({ A: 'mother', G: 'father' }));
assert.strictEqual(ctx.originAt('AG', 'AG', 'AG'), null);
// carrier status
assert.strictEqual(ctx.carrierStatus({ k: 1 }, 'AR', 'male').key, 'carrier');
assert.strictEqual(ctx.carrierStatus({ k: 1 }, 'XL', 'male').key, 'affected');
assert.strictEqual(ctx.carrierStatus({ k: 1 }, 'AD', null).key, 'dominant');

// mtDNA haplogroup on a toy tree: root R, child A (100G), grandchild A1 (100G + 200T), child B (300C)
const ref = 'A'.repeat(400);
const tree = [['R', -1, ''], ['A', 0, '100G'], ['A1', 1, '200T'], ['B', 0, '300C']];
const cons = {}; for (let i = 1; i <= 400; i++) cons[i] = 'A';
cons[100] = 'G'; cons[200] = 'T';
assert.strictEqual(ctx.mtHaplogroup({ consensus: cons, polys: ['100G', '200T'], covered: 400 }, tree, {}, ref).name, 'A1');
// Y: derived calls on R1b and below choose the deepest consistent branch
const ysnps = [[10, 'C', 'T', 'R', 'M207', 'C'], [20, 'G', 'A', 'R1b', 'M343', 'G'], [30, 'A', 'G', 'R1b1', 'P297', 'A'], [40, 'C', 'T', 'R1a', 'M198', 'C']];
const yobs = new Map([['Y:10', { T: 1 }], ['Y:20', { A: 2 }], ['Y:30', { G: 1 }], ['Y:40', { C: 1 }]]);
assert.strictEqual(ctx.yHaplogroup(yobs, ysnps).name, 'R1b1');
// chromosome dosage: X at half density, Y present -> male
const depth = []; for (let c = 0; c < 22; c++) depth.push({ seq: c, start: 1, span: 1e6, n: 1000 });
depth.push({ seq: 22, start: 1, span: 1e6, n: 500 }, { seq: 23, start: 1, span: 1e6, n: 400 });
const names = [...Array(22)].map((_, i) => String(i + 1)).concat(['X', 'Y']);
assert.strictEqual(ctx.chromosomeDosage(depth, names).sex, 'male');

console.log(`ok — ${SNP_DB.length} markers, ${Object.keys(CATEGORIES).length} categories`);
