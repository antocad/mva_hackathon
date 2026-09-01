# Track 2 Drug Repurposing Evaluation: Ataluren for *BUB1B* Nonsense Variant

**Proband Candidate Variant:** `BUB1B` c.2210T>G (p.Leu737*)  
**Drug Candidate:** Ataluren (PTC124 / Translarna)  
**Structural Domain:** Pre-Kinase Linker (aa 737) preceding Catalytic Kinase Domain (aa 750–1045)  
**Evaluation Engine:** PyRosetta FastRelax 15Å MoveMap (`ref2015`) + 4-Link Causal Efficacy Framework  

---

## 1. Executive Summary & Verdict

| Feature | Assessment / Value | Impact on Ataluren Pitch |
| :--- | :--- | :--- |
| **Mutation Class** | Premature Stop Codon (p.Leu737*) | **Ideal Drug Target Class** |
| **Stop Codon Identity** | **UGA** (from codon `TTA` $\rightarrow$ `TGA`) | **Favorable** (UGA is the leakiest stop codon) |
| **Readthrough Species** | Leu737Trp (~70%), Leu737Cys (~15%), Leu737Arg (~15%) | Heterogeneous full-length protein mixture |
| **Structural Stability ($\Delta\Delta G$)** | **-2.51 REU** (15Å MoveMap Relax) | **Foldable** (Trp accommodated by 15Å local flexibility) |
| **Transcript NMD Status** | Exon 17 of 23 (>50 nt from final junction) | **High NMD Risk** (Requires NMD inhibitor / rescue) |
| **Overall Verdict** | **Plausible Partial Rescue Candidate** | **Recommended for Experimental Minigene Validation** |

---

## 2. Metric Definitions & Rosetta Energy Units (REU)

### What is REU?
- **Rosetta Energy Unit (REU):** The standard potential energy metric used in the Rosetta and PyRosetta modeling suites.
- **Physical Scale & Calibration:** Rosetta's full-atom scoring function (`ref2015`; Alford et al., *JCTC* 2017) is parameterized such that **$1\text{ REU} \approx 1\text{ kcal/mol}$** of free energy change.
- **Thermodynamic Interpretation:**
  - $\Delta G < 0\text{ REU}$: Favored, stable structural state. Lower (more negative) values indicate greater thermodynamic stability.
  - $\Delta\Delta G = \Delta G_{\text{mutant}} - \Delta G_{\text{wildtype}}$: Free energy change resulting from a mutation.
    - **$\Delta\Delta G < 1.0\text{ REU}$ ($< 1.0\text{ kcal/mol}$):** Within thermal fluctuation noise ($k_B T \approx 0.6\text{ kcal/mol}$ at 300K). Considered **neutral / wild-type-like** (Kellogg et al., *Proteins* 2011).
    - **$1.0\text{ REU} \le \Delta\Delta G < 5.0\text{ REU}$:** **Mild to moderate thermodynamic perturbation** ($1\text{--}5\text{ kcal/mol}$). The protein remains natively folded under physiological conditions within cellular chaperone buffering capacity (Park et al., *JCTC* 2016).
    - **$\Delta\Delta G \ge 5.0\text{ REU}$:** **Severe destabilization / steric clash**, likely driving protein misfolding or degradation.

### Relaxation Metrics
- **Unrelaxed Energy:** Raw score of the unrelaxed mutation, capturing immediate steric overlap prior to side-chain rotamer optimization.
- **15Å Neighborhood Relax $\Delta G$ / $\Delta\Delta G$:** Side-chains and backbone within $15.0\text{ \AA}$ of position 737 flex, capturing all local packing and H-bonding adaptations while keeping distant regions rigid.
- **Weighted Population Ensemble $\Delta\Delta G$:** UGA stop-codon readthrough produces a heterogeneous protein mixture dominated by Tryptophan (~70%), Cysteine (~15%), and Arginine (~15%) (Blanchet et al., *Nat Commun* 2014; Roy et al., *PNAS* 2016). The ensemble energy is:
  $$\Delta\Delta G_{\text{ensemble}} = 0.70 \cdot \Delta\Delta G_{\text{Trp}} + 0.15 \cdot \Delta\Delta G_{\text{Cys}} + 0.15 \cdot \Delta\Delta G_{\text{Arg}}$$

---

## 3. PyRosetta Thermodynamic Stability & Energy Decomposition (15Å Neighborhood Relax)

Relaxation of side-chains and backbones within a 15Å sphere around position 737 demonstrates how local structural flexibility accommodates the Tryptophan substitution inserted during UGA readthrough:

| Variant | Residue | Unrelaxed Energy (REU) | 15Å Relax $\Delta G$ (REU) | 15Å Relax $\Delta\Delta G$ (REU) |
| :--- | :---: | :---: | :---: | :---: |
| **Wild-Type** | Leu737 | -105.22 | -191.23 | **0.00** |
| **Major Readthrough** | Trp737 | 219.18 | -194.80 | **-3.57** |
| **Minor Readthrough** | Cys737 | -106.79 | -192.04 | **-0.81** |
| **Minor Readthrough** | Arg737 | -106.66 | -190.52 | **+0.71** |

*Weighted Population Ensemble $\Delta\Delta G$ (70% Trp, 15% Cys, 15% Arg):* **-2.51 REU**

---

## 4. Four-Link Efficacy Decision Chain

### 1. Transcript Survival & NMD Evasion
- **Claim:** Does mRNA survive NMD degradation?
- **Verdict:** `STRONG RISK (NMD Substrate)`
- **Literature Source:** *Nagy & Maquat, TIBS 1998 (50-55 nt rule); Popp & Maquat, Cell 2016*

### 2. Ataluren-Induced Readthrough at UGA
- **Claim:** Does Ataluren induce translation readthrough at this specific stop codon?
- **Verdict:** `FAVORABLE (Optimal Stop Codon Identity)`
- **Literature Source:** *Welch et al., Nature 2007; Manuvakhova et al., RNA 2000; PMC9272376*

### 3. Pharmacokinetics & Tissue Distribution
- **Claim:** Does Ataluren reach relevant tissues at therapeutic concentration?
- **Verdict:** `PLAUSIBLE BUT UNTESTED IN MVA`
- **Literature Source:** *Hirawat et al., J Clin Pharmacol 2007; Kerem et al., Lancet 2008*

### 4. Rescued Protein Folding & Kinase Domain Retention
- **Claim:** Does the readthrough protein (L737W/C/R) fold and retain function?
- **Verdict:** `THERMODYNAMICALLY STABLE & FOLDABLE`
- **Literature Source:** *Alford et al., JCTC 2017 (ref2015); Kellogg et al., Proteins 2011*

