// Self-check: node test.js  — validates db.js and the scoring, parsing and comparison code in app.js.
const fs = require('fs'), vm = require('vm'), assert = require('assert');
const ctx = vm.createContext({});
vm.runInContext(fs.readFileSync(__dirname + '/db.js', 'utf8'), ctx);
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

console.log(`ok — ${SNP_DB.length} markers, ${Object.keys(CATEGORIES).length} categories`);
