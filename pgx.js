// Pharmacogenomics: gene calls and drug guidance, summarised from CPIC guidelines (cpicpgx.org) and, where CPIC
// has none, the Dutch Pharmacogenetics Working Group (DPWG). Alleles are on the GRCh37 forward strand.
// Calls from consumer data cover only the listed variants: rarer alleles, gene deletions and duplications
// (important for CYP2D6) are not detected, so "Normal" means "none of the tested variants found".
//
// Gene fields
//   variants  { id, alt, name, value }  value = activity of that allele (activity genes) or 0 = no function
//   kind      'activity'  sum allele activities (1 = normal allele), map with bands
//             'count'     number of variant alleles -> phenotype[count]
//             'cyp2c19'   CPIC CYP2C19 rules (*2/*3 no function, *17 increased)
//   bands     [min activity, phenotype] checked from the top
//   note      shown under the gene

const PGX_GENES = [
  { gene: 'CYP2C19', kind: 'cyp2c19', ref: 'CPIC 2022 (clopidogrel), 2023 (SSRIs), 2020 (PPIs), 2016 (voriconazole)',
    variants: [
      { id: 'rs4244285', alt: 'A', name: '*2', value: 0, proxy: { id: 'rs12769205', alt: 'G' } },
      { id: 'rs4986893', alt: 'A', name: '*3', value: 0 },
      { id: 'rs12248560', alt: 'T', name: '*17', value: 'inc' },
    ] },
  { gene: 'CYP2C9', kind: 'activity', ref: 'CPIC 2020 (NSAIDs), 2020 (phenytoin), 2017 (warfarin)',
    variants: [{ id: 'rs1799853', alt: 'T', name: '*2', value: 0.5 }, { id: 'rs1057910', alt: 'C', name: '*3', value: 0 }],
    bands: [[2, 'Normal metabolizer'], [1, 'Intermediate metabolizer'], [0, 'Poor metabolizer']] },
  { gene: 'CYP2D6', kind: 'activity', ref: 'CPIC 2019 (codeine, tramadol), 2018 (tamoxifen), 2023 (SSRIs), 2016 (TCAs)',
    variants: [
      { id: 'rs3892097', alt: 'T', name: '*4', value: 0 },
      { id: 'rs1065852', alt: 'A', name: '*10', value: 0.25, notWith: 'rs3892097' },
      { id: 'rs28371725', alt: 'T', name: '*41', value: 0.5 },
    ],
    bands: [[1.25, 'Normal metabolizer'], [0.25, 'Intermediate metabolizer'], [0, 'Poor metabolizer']],
    note: 'Gene deletions (*5) and duplications cannot be seen in this data, so ultrarapid metabolizers and some poor metabolizers are missed.' },
  { gene: 'DPYD', kind: 'activity', ref: 'CPIC 2017, updated 2023 (fluorouracil, capecitabine)',
    variants: [
      { id: 'rs3918290', alt: 'T', name: '*2A', value: 0 },
      { id: 'rs55886062', alt: 'C', name: '*13', value: 0 },
      { id: 'rs67376798', alt: 'A', name: 'c.2846A>T', value: 0.5 },
      { id: 'rs75017182', alt: 'C', name: 'HapB3', value: 0.5 },
    ],
    bands: [[2, 'Normal metabolizer'], [1, 'Intermediate metabolizer'], [0, 'Poor metabolizer']] },
  { gene: 'TPMT', kind: 'count', ref: 'CPIC 2018 (thiopurines)',
    variants: [
      { id: 'rs1800462', alt: 'G', name: '*2' },
      { id: 'rs1800460', alt: 'T', name: '*3B' },
      { id: 'rs1142345', alt: 'C', name: '*3C' },
    ],
    cis: ['rs1800460', 'rs1142345'], // *3A carries both: one copy of each counts as one allele
    phenotype: ['Normal metabolizer', 'Intermediate metabolizer', 'Poor metabolizer'] },
  { gene: 'NUDT15', kind: 'count', ref: 'CPIC 2018 (thiopurines)',
    variants: [{ id: 'rs116855232', alt: 'T', name: '*3' }],
    phenotype: ['Normal metabolizer', 'Intermediate metabolizer', 'Poor metabolizer'] },
  { gene: 'SLCO1B1', kind: 'count', ref: 'CPIC 2022 (statins)',
    variants: [{ id: 'rs4149056', alt: 'C', name: 'c.521T>C (*5/*15)' }],
    phenotype: ['Normal function', 'Decreased function', 'Poor function'] },
  { gene: 'ABCG2', kind: 'count', ref: 'CPIC 2022 (rosuvastatin)',
    variants: [{ id: 'rs2231142', alt: 'T', name: 'Q141K' }],
    phenotype: ['Normal function', 'Decreased function', 'Poor function'] },
  { gene: 'VKORC1', kind: 'count', ref: 'CPIC 2017 (warfarin)',
    variants: [{ id: 'rs9923231', alt: 'T', name: '-1639G>A', proxy: { id: 'rs9934438', alt: 'A' } }],
    phenotype: ['Normal warfarin sensitivity', 'Increased warfarin sensitivity', 'High warfarin sensitivity'] },
  { gene: 'CYP4F2', kind: 'count', ref: 'CPIC 2017 (warfarin)',
    variants: [{ id: 'rs2108622', alt: 'T', name: '*3' }],
    phenotype: ['Normal', 'Slightly higher warfarin need', 'Higher warfarin need'] },
  { gene: 'CYP3A5', kind: 'expresser', ref: 'CPIC 2015 (tacrolimus)',
    variants: [{ id: 'rs776746', alt: 'C', name: '*3' }, { id: 'rs10264272', alt: 'T', name: '*6' }],
    phenotype: ['Normal metabolizer (expresser)', 'Intermediate metabolizer (expresser)', 'Poor metabolizer (non-expresser)'],
    note: 'Non-expresser is the usual European result.' },
  { gene: 'UGT1A1', kind: 'count', ref: 'CPIC 2016 (atazanavir); DPWG (irinotecan)',
    variants: [{ id: 'rs887829', alt: 'T', name: '*80 (tags *28)' }],
    phenotype: ['Normal metabolizer', 'Intermediate metabolizer', 'Poor metabolizer (Gilbert)'] },
  { gene: 'CYP2B6', kind: 'count', ref: 'CPIC 2019 (efavirenz)',
    variants: [{ id: 'rs3745274', alt: 'T', name: '516G>T (*6/*9)' }],
    phenotype: ['Normal metabolizer', 'Intermediate metabolizer', 'Poor metabolizer'] },
  { gene: 'G6PD', kind: 'count', ref: 'CPIC 2022 (G6PD)', xlinked: true,
    variants: [{ id: 'rs1050828', alt: 'T', name: 'A- (202G>A)' }, { id: 'rs5030868', alt: 'A', name: 'Mediterranean' }],
    phenotype: ['Normal', 'Variable (heterozygous female)', 'Deficient'] },
  { gene: 'HLA-B*57:01', kind: 'count', ref: 'CPIC 2014 (abacavir)',
    variants: [{ id: 'rs2395029', alt: 'G', name: 'HCP5 tag' }],
    phenotype: ['Not carried', 'Carried', 'Carried'], note: 'Tag SNP; confirm with HLA typing before relying on it.' },
  { gene: 'HLA-A*31:01', kind: 'count', ref: 'CPIC 2017 (carbamazepine)',
    variants: [{ id: 'rs1061235', alt: 'T', name: 'tag' }],
    phenotype: ['Not carried', 'Carried', 'Carried'], note: 'Tag SNP, validated in Europeans; confirm with HLA typing.' },
  { gene: 'MT-RNR1', kind: 'mito', ref: 'CPIC 2021 (aminoglycosides)',
    variants: [{ pos: 1555, alt: 'G', name: 'm.1555A>G' }, { pos: 1494, alt: 'T', name: 'm.1494C>T' }],
    phenotype: ['Normal risk', 'Increased risk of aminoglycoside-induced hearing loss'],
    note: 'Mitochondrial: read from a CRAM file only.' },
];

// Drug guidance: for each drug, the genes it depends on and a rule returning [action, advice].
// action: 'standard' | 'adjust' | 'avoid'
const PGX_DRUGS = [
  { drug: 'Clopidogrel', area: 'Heart', genes: ['CYP2C19'], rule: p => /Poor|Intermediate/.test(p.CYP2C19)
    ? ['avoid', 'Reduced activation: prefer prasugrel or ticagrelor if there is no contraindication.']
    : ['standard', 'Standard dose.'] },
  { drug: 'Citalopram, escitalopram', area: 'Depression', genes: ['CYP2C19'], rule: p => /Poor/.test(p.CYP2C19)
    ? ['adjust', 'Consider a 50% lower starting dose and slower increases, or another antidepressant.']
    : /Ultrarapid|Rapid/.test(p.CYP2C19) ? ['adjust', 'May clear the drug too fast: consider an antidepressant not mainly broken down by CYP2C19.']
    : ['standard', 'Standard dose.'] },
  { drug: 'Sertraline', area: 'Depression', genes: ['CYP2C19'], rule: p => /Poor/.test(p.CYP2C19)
    ? ['adjust', 'Consider a lower starting dose and slower increases.'] : ['standard', 'Standard dose.'] },
  { drug: 'Omeprazole, lansoprazole, pantoprazole', area: 'Stomach', genes: ['CYP2C19'], rule: p => /Ultrarapid/.test(p.CYP2C19)
    ? ['adjust', 'Increase the starting dose by 100%.'] : /Rapid|Normal/.test(p.CYP2C19) ? ['standard', 'Standard dose; for H. pylori or erosive oesophagitis consider a 50–100% higher dose.']
    : ['standard', 'Standard starting dose; for long-term use (over 12 weeks) consider a 50% reduction.'] },
  { drug: 'Voriconazole', area: 'Fungal infection', genes: ['CYP2C19'], rule: p => /Ultrarapid|Rapid|Poor/.test(p.CYP2C19)
    ? ['avoid', 'Choose another antifungal (e.g. isavuconazole) or dose by blood levels.'] : ['standard', 'Standard dose.'] },
  { drug: 'Ibuprofen, celecoxib, flurbiprofen', area: 'Pain', genes: ['CYP2C9'], rule: p => p.CYP2C9_AS <= 0.5
    ? ['adjust', 'Start at 25–50% of the lowest dose, or choose another painkiller.'] : p.CYP2C9_AS === 1
    ? ['adjust', 'Start at the lowest dose and increase carefully.'] : ['standard', 'Standard dose.'] },
  { drug: 'Meloxicam', area: 'Pain', genes: ['CYP2C9'], rule: p => p.CYP2C9_AS <= 0.5
    ? ['avoid', 'Choose another painkiller.'] : p.CYP2C9_AS === 1 ? ['adjust', 'Start at 50% of the lowest dose.'] : ['standard', 'Standard dose.'] },
  { drug: 'Phenytoin', area: 'Epilepsy', genes: ['CYP2C9'], rule: p => p.CYP2C9_AS <= 0.5
    ? ['adjust', 'Reduce the maintenance dose by 50% and monitor levels.'] : p.CYP2C9_AS <= 1.5
    ? ['adjust', 'Reduce the maintenance dose by 25% and monitor levels.'] : ['standard', 'Standard dose.'] },
  { drug: 'Warfarin', area: 'Blood thinner', genes: ['CYP2C9', 'VKORC1', 'CYP4F2'], rule: p => {
      const lower = [];
      if (p.CYP2C9_AS < 2) lower.push('CYP2C9');
      if (/Increased|High/.test(p.VKORC1)) lower.push('VKORC1');
      return lower.length ? ['adjust', `Lower than average dose likely (${lower.join(' and ')}); use a genotype-guided dosing algorithm (e.g. warfarindosing.org).${/higher/i.test(p.CYP4F2) ? ' CYP4F2*3 raises the need by about 5–10%.' : ''}`]
        : /higher/i.test(p.CYP4F2) ? ['adjust', 'About 5–10% higher dose than the algorithm without CYP4F2.'] : ['standard', 'Average dose expected; dose by INR as usual.']; } },
  { drug: 'Codeine, tramadol', area: 'Pain', genes: ['CYP2D6'], rule: p => /Poor/.test(p.CYP2D6)
    ? ['avoid', 'Little or no pain relief: use a non-opioid or an opioid not activated by CYP2D6 (e.g. morphine).'] : /Intermediate/.test(p.CYP2D6)
    ? ['standard', 'Standard dose; if pain relief is poor, switch rather than increase.'] : ['standard', 'Standard dose (ultrarapid metabolism cannot be ruled out from this data).'] },
  { drug: 'Tamoxifen', area: 'Breast cancer', genes: ['CYP2D6'], rule: p => /Poor|Intermediate/.test(p.CYP2D6)
    ? ['adjust', 'Lower endoxifen levels: consider an aromatase inhibitor (post-menopause) and avoid CYP2D6 inhibitors.'] : ['standard', 'Standard dose.'] },
  { drug: 'Amitriptyline, nortriptyline', area: 'Depression, pain', genes: ['CYP2D6'], rule: p => /Poor/.test(p.CYP2D6)
    ? ['avoid', 'Choose another drug; if needed, start at 50% of the dose and monitor levels.'] : /Intermediate/.test(p.CYP2D6)
    ? ['adjust', 'Consider a 25% lower starting dose.'] : ['standard', 'Standard dose.'] },
  { drug: 'Paroxetine, fluvoxamine', area: 'Depression', genes: ['CYP2D6'], rule: p => /Poor/.test(p.CYP2D6)
    ? ['adjust', 'Consider a 50% lower starting dose and slower increases.'] : ['standard', 'Standard dose.'] },
  { drug: 'Atomoxetine', area: 'ADHD', genes: ['CYP2D6'], rule: p => /Poor/.test(p.CYP2D6)
    ? ['adjust', 'Start at the usual dose but increase only if response is poor after 2 weeks; monitor side effects.'] : ['standard', 'Standard dose.'] },
  { drug: 'Metoprolol', area: 'Heart', genes: ['CYP2D6'], rule: p => /Poor/.test(p.CYP2D6)
    ? ['adjust', 'Use a lower dose or a beta-blocker not broken down by CYP2D6 (e.g. bisoprolol) (DPWG).'] : ['standard', 'Standard dose.'] },
  { drug: 'Fluorouracil, capecitabine', area: 'Cancer', genes: ['DPYD'], rule: p => p.DPYD_AS <= 0.5
    ? ['avoid', 'Avoid: severe or fatal toxicity. If no alternative, a strongly reduced dose with drug-level monitoring.'] : p.DPYD_AS < 2
    ? ['adjust', 'Start at 50% of the dose, then increase if tolerated.'] : ['standard', 'Standard dose.'] },
  { drug: 'Azathioprine, mercaptopurine, thioguanine', area: 'Immune, leukaemia', genes: ['TPMT', 'NUDT15'], rule: p => /Poor/.test(p.TPMT + p.NUDT15)
    ? ['avoid', 'Drastically reduced dose (about 10% given 3×/week) for cancer; choose another drug for other conditions.'] : /Intermediate/.test(p.TPMT + p.NUDT15)
    ? ['adjust', 'Start at 30–80% of the normal dose and adjust by blood counts.'] : ['standard', 'Standard dose.'] },
  { drug: 'Simvastatin', area: 'Cholesterol', genes: ['SLCO1B1'], rule: p => /Poor/.test(p.SLCO1B1)
    ? ['avoid', 'High myopathy risk: choose another statin.'] : /Decreased/.test(p.SLCO1B1)
    ? ['adjust', 'Keep the dose at 20 mg/day or less, or choose another statin.'] : ['standard', 'Standard dose.'] },
  { drug: 'Atorvastatin', area: 'Cholesterol', genes: ['SLCO1B1'], rule: p => /Poor/.test(p.SLCO1B1)
    ? ['adjust', 'Start at 20 mg/day or less; consider rosuvastatin or pravastatin for higher needs.'] : ['standard', 'Standard dose.'] },
  { drug: 'Rosuvastatin', area: 'Cholesterol', genes: ['SLCO1B1', 'ABCG2'], rule: p => /Poor/.test(p.SLCO1B1) || /Poor/.test(p.ABCG2)
    ? ['adjust', 'Start at 10 mg/day or less (higher exposure).'] : ['standard', 'Standard dose.'] },
  { drug: 'Tacrolimus', area: 'Transplant', genes: ['CYP3A5'], rule: p => /expresser\)/.test(p.CYP3A5) && !/non-expresser/.test(p.CYP3A5)
    ? ['adjust', 'Increase the starting dose 1.5–2× and dose by blood levels.'] : ['standard', 'Standard starting dose.'] },
  { drug: 'Atazanavir', area: 'HIV', genes: ['UGT1A1'], rule: p => /Poor/.test(p.UGT1A1)
    ? ['adjust', 'High chance of jaundice (harmless but may stop treatment): consider another drug.'] : ['standard', 'Standard dose.'] },
  { drug: 'Irinotecan', area: 'Cancer', genes: ['UGT1A1'], rule: p => /Poor/.test(p.UGT1A1)
    ? ['adjust', 'Start at 70% of the dose (DPWG) and monitor neutrophils.'] : ['standard', 'Standard dose.'] },
  { drug: 'Efavirenz', area: 'HIV', genes: ['CYP2B6'], rule: p => /Poor/.test(p.CYP2B6)
    ? ['adjust', 'Start at 400 or 200 mg/day.'] : /Intermediate/.test(p.CYP2B6) ? ['adjust', 'Consider 400 mg/day.'] : ['standard', 'Standard dose.'] },
  { drug: 'Rasburicase, primaquine, dapsone', area: 'Various', genes: ['G6PD'], rule: p => /Deficient/.test(p.G6PD)
    ? ['avoid', 'Risk of haemolysis: avoid (rasburicase contraindicated); check G6PD activity.'] : /Variable/.test(p.G6PD)
    ? ['adjust', 'Check G6PD enzyme activity before use.'] : ['standard', 'No G6PD variant found.'] },
  { drug: 'Abacavir', area: 'HIV', genes: ['HLA-B*57:01'], rule: p => /^Carried/.test(p['HLA-B*57:01'])
    ? ['avoid', 'Hypersensitivity risk: do not use (confirm with HLA typing).'] : ['standard', 'No HLA-B*57:01 tag found; standard pre-treatment HLA test still applies.'] },
  { drug: 'Carbamazepine, oxcarbazepine', area: 'Epilepsy', genes: ['HLA-A*31:01'], rule: p => /^Carried/.test(p['HLA-A*31:01'])
    ? ['avoid', 'Higher risk of severe skin reactions (DRESS, SJS): consider another drug.'] : ['standard', 'No HLA-A*31:01 tag found (HLA-B*15:02, relevant in Asian ancestry, is not tested).'] },
  { drug: 'Gentamicin, amikacin, tobramycin', area: 'Infection', genes: ['MT-RNR1'], rule: p => /Increased/.test(p['MT-RNR1'] || '')
    ? ['avoid', 'Avoid aminoglycosides unless the infection is life-threatening and no alternative exists.'] : ['standard', 'No MT-RNR1 risk variant found.'] },
];
