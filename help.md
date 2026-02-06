How it all works
The pipeline takes your mouse Micro-C contact matrix (.mcool file) and extracts where TAD boundaries are, how strong they are, and whether they're real — in three stages:

Stage 1: Detect boundaries (two complementary methods)
Insulation Score — slides a diamond-shaped window along the diagonal of the contact map. Inside a TAD, the window captures lots of contacts (high insulation). At a boundary, contacts drop sharply because the window straddles two domains. The output is a 1D track where valleys = boundaries.

Directionality Index (DI) — for each genomic bin, compares "does this bin talk more to its left neighbors or right neighbors?" Inside a TAD, bins near the left edge talk rightward (positive DI), bins near the right edge talk leftward (negative DI). A sharp sign flip from negative to positive = boundary. This catches boundaries that insulation might miss when the valley is shallow.

Stage 2: Score and call TADs
Boundary Prominence — not all insulation valleys are equal. Some are deep canyons (strong CTCF-anchored boundaries), others are shallow dips (noise or transient sub-structures). The code uses scipy.signal.find_peaks to measure the topographic prominence of each valley — how deep it is relative to the surrounding landscape. Then classifies:

Class	Prominence	Meaning
Strong	>= 0.5	Constitutive boundary, likely CTCF/cohesin-anchored
Weak	0.2 – 0.5	Possible tissue-specific or transient domain edge
Sub-threshold	< 0.2	Probably noise
TAD Interval Calling — takes the boundary positions and says "the genomic segment between boundary A and boundary B is one TAD." Outputs a BED-like table of TAD intervals with start/end coordinates. Filters out anything suspiciously large (>3Mb) since that's likely a gap, not a real domain.

Multi-Scale Detection — runs insulation at many window sizes (25kb through 500kb). A boundary that appears at all scales is constitutive (a major structural barrier). One that only appears at small windows is a nested sub-TAD boundary. One that only appears at large windows is a meta-TAD boundary.

Stage 3: Validate with pileup
Boundary Pileup — the "does this actually look like a real boundary?" test. It takes every detected boundary, cuts out a square sub-matrix centered on it from the contact map, normalizes by expected distance decay, and averages them all together. If your boundaries are real, the averaged heatmap should show:


 High  |  Low
-------+-------  
 Low   |  High
The top-left and bottom-right quadrants (intra-domain) should be enriched. The off-diagonal quadrants (inter-domain) should be depleted. This is the classic "TAD corner" signal.

What each output plot shows
Plot	What to look for
_directionality_index.png	Red/blue oscillations. Sharp sign changes align with domain edges.
_boundary_strength.png	Bar chart ranked by prominence. A clear drop-off separates real boundaries from noise.
_multiscale_insulation.png	Vertical blue streaks = boundaries that persist across scales (strong). Streaks only at the top or bottom = scale-specific boundaries.
_tad_overlay.png	Hi-C heatmap with black L-shaped corners marking each called TAD. Do the corners match the visible domain triangles in the map?
_boundary_pileup.png	The checkerboard pattern described above. Strong quadrant contrast = your boundaries are biologically real.
_combined_boundary.png	Everything together — heatmap + insulation + DI + boundary ticks — so you can visually confirm that all signals agree at the same positions.
The key question is always: do the insulation valleys, the DI sign flips, and the visible domain triangles in the contact map all line up at the same genomic positions? If yes, those are high-confidence TAD boundaries.