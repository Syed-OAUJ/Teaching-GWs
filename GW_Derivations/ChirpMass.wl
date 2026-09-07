(* ::Package:: *)

(* ::Package:: *)

(* ::Title:: *)
(*Chirp Mass Scaling from the Quadrupole Formula*)

(* ::Subtitle:: *)
(*Symbolic verification, numerical inspiral integration, and visualization companion to the LaTeX derivation.*)

(* ::Text:: *)
(*Open this file in Mathematica (it will auto-convert to a notebook with Section/Text/Input cells), or run it top-to-bottom with wolframscript. Evaluate section by section.*)


(* ::Section:: *)
(*1. Symbolic verification: quadrupole moment and its third time derivative*)


(* ::Text:: *)
(*We reproduce, symbolically, the hand derivation: a reduced mass \[Mu] on a circular orbit of radius a and angular frequency \[Omega], and its trace-free quadrupole moment Q_ij(t). We then verify that Overscript[Q, ...]_ij Overscript[Q, ...]^ij is exactly constant and equal to 32 \[Mu]^2 a^4 \[Omega]^6, as claimed in the notes.*)


ClearAll[xt, yt, zt, Qmat, sym\[Mu], syma, sym\[Omega], t];

xt[t_] := syma Cos[sym\[Omega] t];
yt[t_] := syma Sin[sym\[Omega] t];
zt[t_] := 0;

r2 = xt[t]^2 + yt[t]^2 + zt[t]^2 // Simplify;

posVec = {xt[t], yt[t], zt[t]};
Qmat = Table[
   sym\[Mu] (posVec[[i]] posVec[[j]] - (1/3) KroneckerDelta[i, j] r2),
   {i, 1, 3}, {j, 1, 3}
   ] // TrigReduce // Simplify;

Print["Quadrupole tensor Q_ij(t):"];
Print[MatrixForm[Qmat]];
Print["Trace (should be 0): ", Simplify[Tr[Qmat]]];

(* Third time derivative of each component *)
Q3 = D[Qmat, {t, 3}] // TrigReduce // Simplify;
Print["Third time derivative dddot(Q)_ij(t):"];
Print[MatrixForm[Q3]];

(* Contract: sum_ij (Q3_ij)^2 *)
contraction = Sum[Q3[[i, j]]^2, {i, 1, 3}, {j, 1, 3}] // TrigReduce // Simplify;
Print["Contraction dddotQ_ij dddotQ^ij  =  ", contraction];
Print["Matches hand-derived 32 mu^2 a^4 omega^6 ?  ",
  Simplify[contraction - 32 sym\[Mu]^2 syma^4 sym\[Omega]^6] === 0];


(* ::Section:: *)
(*2. Symbolic checks: chirp-mass algebraic identities*)


(* ::Text:: *)
(*Verify the two identities used to collapse (mu, M) into the single combination \[ScriptCapitalM] = mu^(3/5) M^(2/5):*)
(*  (i)  \[ScriptCapitalM]^(10/3) = mu^2 M^(4/3)          [used to get P_GW(omega) in terms of \[ScriptCapitalM]]*)
(*  (ii) \[ScriptCapitalM]^(5/3)  = mu   M^(2/3)          [used to get omegadot in terms of \[ScriptCapitalM]]*)


ClearAll[symM, symMc];
symMc = sym\[Mu]^(3/5) symM^(2/5);

check1 = Simplify[symMc^(10/3) - sym\[Mu]^2 symM^(4/3),
   Assumptions -> {sym\[Mu] > 0, symM > 0}];
check2 = Simplify[symMc^(5/3) - sym\[Mu] symM^(2/3),
   Assumptions -> {sym\[Mu] > 0, symM > 0}];

Print["Identity (i)  Mc^(10/3) - mu^2 M^(4/3) = 0 ?  ", check1 === 0];
Print["Identity (ii) Mc^(5/3)  - mu   M^(2/3) = 0 ?  ", check2 === 0];


(* ::Section:: *)
(*3. Symbolic derivation of the chirp equation fdot(f)*)


(* ::Text:: *)
(*Redo the energy-balance step symbolically: P_GW(omega) from Section 1's contraction plus Kepler's law, orbital energy E(omega), then solve dE/dt = -P_GW for omegadot.*)


ClearAll[G, c, PGW, Eorb, om, \[Omega]dotSol];

(* P_GW(a,omega) from the quadrupole contraction *)
PGWaOmega = (G/(5 c^5)) contraction /. {syma -> syma, sym\[Omega] -> sym\[Omega]} // Simplify;
Print["P_GW(a,\[Omega]) = ", PGWaOmega];

(* Eliminate a via Kepler: omega^2 = G M / a^3  =>  a = (G M/omega^2)^(1/3) *)
PGWomega = PGWaOmega /. syma -> (G symM/sym\[Omega]^2)^(1/3) // PowerExpand // Simplify;
Print["P_GW(\[Omega]) after Kepler substitution = ", PGWomega];

(* Orbital energy E(omega) = -(mu/2)(GM)^(2/3) omega^(2/3) *)
EorbOmega = -(sym\[Mu]/2) (G symM)^(2/3) sym\[Omega]^(2/3);

(* Energy balance: dE/domega * omegadot = - P_GW(omega)  =>  solve for omegadot *)
dEdomega = D[EorbOmega, sym\[Omega]];
\[Omega]dotSol = Simplify[-PGWomega/dEdomega,
   Assumptions -> {G > 0, symM > 0, sym\[Mu] > 0, sym\[Omega] > 0}];
Print["\[Omega]dot(\[Omega]) = ", \[Omega]dotSol];

(* Re-express using the chirp mass identity Mc^(10/3) = mu^2 M^(4/3) *)
\[Omega]dotInMc = Simplify[
   \[Omega]dotSol /. sym\[Mu]^2 symM^(4/3) -> symMc^(10/3),
   Assumptions -> {G > 0, symM > 0, sym\[Mu] > 0, sym\[Omega] > 0}];
Print["\[Omega]dot in terms of Mc (should reduce to (96/5) G^(5/3)/c^5 Mc^(5/3) \[Omega]^(11/3)):"];
Print[\[Omega]dotInMc];


(* ::Section:: *)
(*4. Numerical setup: a GW150914-like binary*)


ClearAll["Global`*"];

(* Physical constants, SI units *)
G = 6.67430*10^-11;      (* m^3 kg^-1 s^-2 *)
c = 2.99792458*10^8;      (* m/s *)
Msun = 1.98892*10^30;     (* kg *)
Mpc = 3.0857*10^22;       (* m *)

(* Component masses -- GW150914-like *)
m1 = 36 Msun;
m2 = 29 Msun;
M = m1 + m2;
\[Mu]num = m1 m2/M;
Mc = \[Mu]num^(3/5) M^(2/5);

Print["m1 = ", m1/Msun, " Msun,  m2 = ", m2/Msun, " Msun"];
Print["Total mass M = ", M/Msun, " Msun"];
Print["Chirp mass \[ScriptCapitalM] = ", Mc/Msun, " Msun"];

(* Start the inspiral when the GW frequency enters the ground-based band *)
f0 = 20;                    (* Hz *)
\[Omega]0 = \[Pi] f0;
a0 = (G M/\[Omega]0^2)^(1/3);
Print["Initial separation a0 = ", a0/1000, " km  =  ", a0/(G M/c^2), " (in units of GM/c^2)"];

(* Analytic (Peters) merger-time estimate, Eq. (t_merge) from the notes *)
tMergeEst = (5/256) c^5 a0^4/(G^3 \[Mu]num M^2);
Print["Analytic estimate of merger time t_merge = ", tMergeEst, " s"];

(* Schwarzschild-like ISCO cutoff, purely as a numerical stopping radius *)
rISCO = 6 G M/c^2;
Print["ISCO cutoff radius = ", rISCO/1000, " km"];


(* ::Section:: *)
(*5. Numerically integrate the inspiral: da/dt and the orbital phase*)


(* ::Text:: *)
(*Coupled ODEs: da/dt from the Peters formula, and orbital phase \[Phi]orb'[t] = \[Omega]orb(t) = Sqrt[G M / a(t)^3], integrated together so we get both a(t) and the orbital phase (hence the GW phase = 2 \[Phi]orb) directly as InterpolatingFunctions.*)


tMaxTry = 1.02 tMergeEst;

solList = NDSolve[
   {
    a'[\[Tau]] == -(64/5) (G^3 \[Mu]num M^2)/(c^5 a[\[Tau]]^3),
    \[Phi]orb'[\[Tau]] == Sqrt[G M/a[\[Tau]]^3],
    a[0] == a0,
    \[Phi]orb[0] == 0,
    WhenEvent[a[\[Tau]] <= rISCO, "StopIntegration"]
   },
   {a, \[Phi]orb}, {\[Tau], 0, tMaxTry},
   MaxStepFraction -> 10^-4
   ];

aFunc = a /. First[solList];
\[Phi]orbFunc = \[Phi]orb /. First[solList];

tNumMerge = aFunc["Domain"][[1, 2]];
Print["Numerically integrated inspiral duration (down to ISCO cutoff) = ", tNumMerge, " s"];
Print["Ratio to analytic t_merge estimate: ", tNumMerge/tMergeEst];

fFunc[\[Tau]_] := Sqrt[G M/aFunc[\[Tau]]^3]/\[Pi];          (* GW frequency f = omega_orb/pi *)
\[CapitalPhi]GWFunc[\[Tau]_] := 2 \[Phi]orbFunc[\[Tau]];      (* GW phase = 2 x orbital phase *)


(* ::Section:: *)
(*6. Plot: separation a(t) and GW frequency f(t), numeric vs. analytic*)


plotA = Plot[aFunc[\[Tau]]/1000, {\[Tau], 0, tNumMerge},
   PlotRange -> All, Frame -> True,
   FrameLabel -> {"t  [s]", "separation a(t)  [km]"},
   PlotLabel -> "Inspiral separation (Peters 1964 ODE)",
   PlotStyle -> Directive[RGBColor[0.20, 0.35, 0.65], Thickness[0.004]],
   ImageSize -> 500];

(* Analytic f(t) from inverting tau(f) = 5/256 pi^(-8/3) (G Mc/c^3)^(-5/3) f^(-8/3) *)
fAnalytic[\[Tau]_] := (1/\[Pi]) ((5/(256 (tMergeEst - \[Tau]))) (G Mc/c^3)^(-5/3))^(3/8);

plotF = Plot[
   {fFunc[\[Tau]], fAnalytic[\[Tau]]},
   {\[Tau], 0, 0.995 tNumMerge},
   PlotRange -> All, Frame -> True,
   FrameLabel -> {"t  [s]", "GW frequency f(t)  [Hz]"},
   PlotLabel -> "Chirp: numeric ODE (blue) vs. analytic fdot formula (dashed red)",
   PlotStyle -> {Directive[RGBColor[0.20, 0.35, 0.65], Thickness[0.004]],
     Directive[RGBColor[0.75, 0.15, 0.15], Dashed, Thickness[0.003]]},
   PlotLegends -> {"a(t) ODE + Kepler", "analytic \!\(\*OverscriptBox[\(f\), \(.\)]\)(f) solution"},
   ImageSize -> 500];

GraphicsRow[{plotA, plotF}, ImageSize -> 1000]


(* ::Section:: *)
(*7. Leading-order (restricted) chirp waveform h+(t)*)


(* ::Text:: *)
(*Face-on, leading quadrupole-order plus-polarization strain:*)
(*  h+(t) = (4/D) (G Mc/c^2)^(5/3) (pi f(t)/c)^(2/3) Cos[Phi_GW(t)]*)
(*This is the standard "restricted" Newtonian waveform used for illustration; it omits amplitude PN corrections and inclination/polarization dependence.*)


D0 = 400 Mpc;  (* luminosity distance, GW150914-like *)

hplus[\[Tau]_] := (4/D0) (G Mc/c^2)^(5/3) (\[Pi] fFunc[\[Tau]]/c)^(2/3) Cos[\[CapitalPhi]GWFunc[\[Tau]]];

(* Zoom on the last ~0.3 s before the ISCO cutoff, where the chirp is visually obvious *)
tZoomStart = Max[0, tNumMerge - 0.30];

plotWaveform = Plot[hplus[\[Tau]], {\[Tau], tZoomStart, 0.999 tNumMerge},
   PlotPoints -> 3000, MaxRecursion -> 4,
   Frame -> True, FrameLabel -> {"t  [s]", "h+(t)"},
   PlotLabel -> "Restricted leading-order chirp waveform (last ~0.3 s)",
   PlotStyle -> Directive[RGBColor[0.1, 0.1, 0.1], Thickness[0.0015]],
   ImageSize -> 700]


(* ::Section:: *)
(*8. Animation: shrinking binary orbit*)


(* ::Text:: *)
(*Positions of each body about the center of mass: r1(t) = (m2/M) a(t) {cos,sin}[phi_orb(t)], r2(t) = -(m1/M) a(t) {cos,sin}[phi_orb(t)].*)


nFrames = 150;
frameTimes = Subdivide[0, 0.995 tNumMerge, nFrames - 1];

orbitFrame[\[Tau]_] := Module[{ang, sep, r1, r2, trail},
   ang = \[Phi]orbFunc[\[Tau]];
   sep = aFunc[\[Tau]];
   r1 = (m2/M) sep {Cos[ang], Sin[ang]};
   r2 = -(m1/M) sep {Cos[ang], Sin[ang]};
   Graphics[
    {
     RGBColor[0.20, 0.35, 0.65], PointSize[0.05], Point[r1/1000],
     RGBColor[0.75, 0.15, 0.15], PointSize[0.04], Point[r2/1000],
     Black, PointSize[0.015], Point[{0, 0}]
    },
    PlotRange -> {{-450, 450}, {-450, 450}},
    Axes -> False, Frame -> True,
    FrameLabel -> {"x  [km]", "y  [km]"},
    PlotLabel -> "Binary inspiral,  t = " <> ToString[NumberForm[\[Tau], {5, 4}]] <> " s",
    ImageSize -> 400
    ]
   ];

orbitFrames = orbitFrame /@ frameTimes;

Export[FileNameJoin[{NotebookDirectory[], "binary_orbit_inspiral.gif"}],
  orbitFrames, "DisplayDurations" -> ConstantArray[0.04, nFrames]];

Print["Exported binary_orbit_inspiral.gif (", nFrames, " frames) to the notebook directory."];

ListAnimate[orbitFrames, AnimationRate -> 25]


(* ::Section:: *)
(*9. Sanity check against the real GW150914 event*)


(* ::Text:: *)
(*LIGO's measured chirp mass for GW150914 was \[ScriptCapitalM] \[TildeTilde] 28 Msun, with a signal duration of order 0.2 s from ~35 Hz to merger. Compare: *)


tauCheck = tMergeEst - 0 // N;
f35 = 35;
tauFrom35Hz = (5/256) \[Pi]^(-8/3) (G Mc/c^3)^(-5/3) f35^(-8/3);

Print["This toy binary's chirp mass \[ScriptCapitalM] = ", Mc/Msun, " Msun  (cf. GW150914 \[TildeTilde] 28 Msun)"];
Print["Time from f = 35 Hz to merger (this system) = ", tauFrom35Hz, " s  (cf. GW150914 \[TildeTilde] 0.2 s)"];



