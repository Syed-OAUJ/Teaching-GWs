(* ::Package:: *)

(* ::Title:: *)
(*Newtonian vs. 1PN Chirp: Side-by-Side Comparison*)


(* ::Subtitle:: *)
(*Companion to chirp_mass_1PN_correction.tex (Part II). Self-contained: does not require running ChirpMassDerivation.wl first.*)


(* ::Text:: *)
(*Open in Mathematica (auto-converts to a notebook), or run with wolframscript. Evaluate section by section. This script has NOT been executed here -- run and debug locally.*)


(* ::Section:: *)
(*1. Constants and the binary (same GW150914-like system as Part I)*)


ClearAll["Global`*"];

G = 6.67430*10^-11;      (* m^3 kg^-1 s^-2 *)
c = 2.99792458*10^8;      (* m/s *)
Msun = 1.98892*10^30;     (* kg *)
Mpc = 3.0857*10^22;       (* m *)

m1 = 36 Msun;
m2 = 29 Msun;
M = m1 + m2;
\[Mu]num = m1 m2/M;
\[Eta] = \[Mu]num/M;
Mc = \[Mu]num^(3/5) M^(2/5);

Print["m1 = ", m1/Msun, " Msun,  m2 = ", m2/Msun, " Msun"];
Print["M = ", M/Msun, " Msun,   \[Eta] = ", \[Eta], "   (\[Eta]=1/4 for equal masses)"];
Print["Chirp mass \[ScriptCapitalM] = ", Mc/Msun, " Msun"];

(* PN expansion parameter x(f) = (pi G M f / c^3)^(2/3) *)
xOfF[f_] := (\[Pi] G M f/c^3)^(2/3);

(* the 1PN bracket appearing in both fdot and tau *)
pn1Coeff = 743/336 + (11/4) \[Eta];
Print["1PN coefficient (743/336 + 11/4 \[Eta]) = ", pn1Coeff // N];


(* ::Section:: *)
(*2. fdot(f): Newtonian vs. 1PN, and the size of the correction*)


fdotN[f_] := (96/5) \[Pi]^(8/3) (G Mc/c^3)^(5/3) f^(11/3);
fdot1PN[f_] := fdotN[f] (1 - pn1Coeff xOfF[f]);

(* Frequency band of interest: 20 Hz (detector low end) up to a rough
   "1PN validity ceiling" where x is still modest, say x=0.15 (~ 200-300 Hz here) *)
fLow = 20;
xCeil = 0.15;
fHigh = fLow;
While[xOfF[fHigh] < xCeil, fHigh *= 1.01];
Print["Frequency range shown: ", fLow, " Hz to ", fHigh // Round, " Hz  (x from ",
  xOfF[fLow] // N, " to ", xOfF[fHigh] // N, ")"];

plotFdot = LogPlot[{fdotN[f], fdot1PN[f]}, {f, fLow, fHigh},
   PlotRange -> All, Frame -> True,
   FrameLabel -> {"f  [Hz]", "\!\(\*OverscriptBox[\(f\), \(.\)]\)  [Hz/s]"},
   PlotLabel -> "Frequency sweep rate: Newtonian (blue) vs. 1PN (dashed red)",
   PlotStyle -> {Directive[RGBColor[0.20, 0.35, 0.65], Thickness[0.004]],
     Directive[RGBColor[0.75, 0.15, 0.15], Dashed, Thickness[0.003]]},
   PlotLegends -> {"Newtonian (0PN)", "1PN-corrected"},
   ImageSize -> 500];

plotFdotRatio = Plot[100 (fdot1PN[f]/fdotN[f] - 1), {f, fLow, fHigh},
   Frame -> True,
   FrameLabel -> {"f  [Hz]", "(\!\(\*SubscriptBox[\(f\), \(1PN\)]\)/\!\(\*SubscriptBox[\(f\), \(N\)]\) - 1)  [%]"},
   PlotLabel -> "Relative 1PN correction to \!\(\*OverscriptBox[\(f\), \(.\)]\)",
   PlotStyle -> Directive[RGBColor[0.15, 0.55, 0.25], Thickness[0.004]],
   ImageSize -> 500];

GraphicsRow[{plotFdot, plotFdotRatio}, ImageSize -> 1000]


(* ::Section:: *)
(*3. Time-to-coalescence tau(f): Newtonian vs. 1PN*)


tauN[f_] := (5/256) (G Mc/c^3)^(-5/3) (\[Pi] f)^(-8/3);
tau1PN[f_] := tauN[f] (1 + (4/3) pn1Coeff xOfF[f]);

plotTau = LogLogPlot[{tauN[f], tau1PN[f]}, {f, fLow, fHigh},
   Frame -> True,
   FrameLabel -> {"f  [Hz]", "\[Tau](f)  [s]"},
   PlotLabel -> "Time to coalescence: Newtonian (blue) vs. 1PN (dashed red)",
   PlotStyle -> {Directive[RGBColor[0.20, 0.35, 0.65], Thickness[0.004]],
     Directive[RGBColor[0.75, 0.15, 0.15], Dashed, Thickness[0.003]]},
   PlotLegends -> {"Newtonian (0PN)", "1PN-corrected"},
   ImageSize -> 500];

plotTauDiff = Plot[tau1PN[f] - tauN[f], {f, fLow, fHigh},
   Frame -> True,
   FrameLabel -> {"f  [Hz]", "\[Tau]\!\(\*SubscriptBox[\(\), \(1PN\)]\) - \[Tau]\!\(\*SubscriptBox[\(\), \(N\)]\)  [s]"},
   PlotLabel -> "Absolute 1PN correction to the time-to-merger estimate",
   PlotStyle -> Directive[RGBColor[0.55, 0.25, 0.55], Thickness[0.004]],
   ImageSize -> 500];

GraphicsRow[{plotTau, plotTauDiff}, ImageSize -> 1000]


(* ::Section:: *)
(*4. Integrate x(t) with and without the 1PN term*)


(* dx/dt_N   = (64 c^3 eta)/(5 G M) x^5
   dx/dt_1PN = (64 c^3 eta)/(5 G M) x^5 [1 - pn1Coeff x]           *)

f0 = 20;
x0 = xOfF[f0];
xISCO = 1/6;  (* standard PN "Schwarzschild ISCO" cutoff, x = 1/6 *)

tMergeEstN = tauN[f0];
tMaxTry = 1.02 tMergeEstN;

solN = NDSolve[
   {
    xN'[\[Tau]] == (64/5) (c^3 \[Eta])/(G M) xN[\[Tau]]^5,
    xN[0] == x0,
    WhenEvent[xN[\[Tau]] >= xISCO, "StopIntegration"]
   },
   xN, {\[Tau], 0, tMaxTry}, MaxStepFraction -> 10^-4];

sol1PN = NDSolve[
   {
    x1'[\[Tau]] == (64/5) (c^3 \[Eta])/(G M) x1[\[Tau]]^5 (1 - pn1Coeff x1[\[Tau]]),
    x1[0] == x0,
    WhenEvent[x1[\[Tau]] >= xISCO, "StopIntegration"]
   },
   x1, {\[Tau], 0, tMaxTry}, MaxStepFraction -> 10^-4];

xNFunc = xN /. First[solN];
x1PNFunc = x1 /. First[sol1PN];

tEndN = xNFunc["Domain"][[1, 2]];
tEnd1PN = x1PNFunc["Domain"][[1, 2]];
Print["Time to reach x=1/6 (ISCO):  Newtonian = ", tEndN, " s,   1PN = ", tEnd1PN, " s"];
Print["Difference = ", tEnd1PN - tEndN, " s  (", 100 (tEnd1PN - tEndN)/tEndN // N, " %)"];

fNFunc[\[Tau]_] := (c^3/(\[Pi] G M)) xNFunc[\[Tau]]^(3/2);
f1PNFunc[\[Tau]_] := (c^3/(\[Pi] G M)) x1PNFunc[\[Tau]]^(3/2);

plotFvsT = Plot[{fNFunc[\[Tau]], f1PNFunc[\[Tau]]}, {\[Tau], 0, Min[tEndN, tEnd1PN]},
   Frame -> True,
   FrameLabel -> {"t  [s]", "f(t)  [Hz]"},
   PlotLabel -> "GW frequency vs. time: Newtonian (blue) vs. 1PN (dashed red)",
   PlotStyle -> {Directive[RGBColor[0.20, 0.35, 0.65], Thickness[0.004]],
     Directive[RGBColor[0.75, 0.15, 0.15], Dashed, Thickness[0.003]]},
   PlotLegends -> {"Newtonian (0PN)", "1PN-corrected"},
   ImageSize -> 700]


(* ::Section:: *)
(*5. Accumulated dephasing between the two models*)


(* ::Text:: *)
(*This is the physically important comparison for matched filtering: even a small fractional difference in f(t) accumulates, over many cycles, into an O(1) radian phase difference -- which is enough to ruin the overlap between a Newtonian-only template and the true (1PN) signal. Re-integrate with the orbital phase as a second ODE variable for each model.*)


solNphase = NDSolve[
   {
    xN2'[\[Tau]] == (64/5) (c^3 \[Eta])/(G M) xN2[\[Tau]]^5,
    \[Phi]N'[\[Tau]] == (c^3/(G M)) xN2[\[Tau]]^(3/2),
    xN2[0] == x0, \[Phi]N[0] == 0,
    WhenEvent[xN2[\[Tau]] >= xISCO, "StopIntegration"]
   },
   {xN2, \[Phi]N}, {\[Tau], 0, tMaxTry}, MaxStepFraction -> 10^-4];

sol1PNphase = NDSolve[
   {
    x12'[\[Tau]] == (64/5) (c^3 \[Eta])/(G M) x12[\[Tau]]^5 (1 - pn1Coeff x12[\[Tau]]),
    \[Phi]1'[\[Tau]] == (c^3/(G M)) x12[\[Tau]]^(3/2),
    x12[0] == x0, \[Phi]1[0] == 0,
    WhenEvent[x12[\[Tau]] >= xISCO, "StopIntegration"]
   },
   {x12, \[Phi]1}, {\[Tau], 0, tMaxTry}, MaxStepFraction -> 10^-4];

\[Phi]NFunc = \[Phi]N /. First[solNphase];
\[Phi]1Func = \[Phi]1 /. First[sol1PNphase];

tCommon = Min[\[Phi]NFunc["Domain"][[1, 2]], \[Phi]1Func["Domain"][[1, 2]]];

(* GW phase = 2 x orbital phase; dephasing in GW cycles = DeltaPhi/(2 pi) *)
dephasing[\[Tau]_] := 2 (\[Phi]1Func[\[Tau]] - \[Phi]NFunc[\[Tau]]);

plotDephasing = Plot[dephasing[\[Tau]], {\[Tau], 0, 0.999 tCommon},
   Frame -> True,
   FrameLabel -> {"t  [s]", "\[CapitalDelta]\[CapitalPhi]\!\(\*SubscriptBox[\(\), \(GW\)]\)(t)  [rad]"},
   PlotLabel -> "Accumulated 1PN dephasing relative to a Newtonian-only template",
   PlotStyle -> Directive[RGBColor[0.75, 0.15, 0.15], Thickness[0.004]],
   ImageSize -> 700];

Print["Total accumulated dephasing over the shown band \[TildeTilde] ",
  dephasing[0.999 tCommon], " rad  =  ",
  dephasing[0.999 tCommon]/(2 \[Pi]) // N, " GW cycles"];
Print["(For reference: a dephasing of order ~1 radian is already enough to significantly ",
  "degrade the matched-filter overlap between a Newtonian-only template and the true signal ",
  "-- this is the practical reason PN order matters for template accuracy.)"];

plotDephasing


(* ::Section:: *)
(*6. Waveform overlay near merger: Newtonian-only vs. 1PN-corrected frequency evolution*)


D0 = 400 Mpc;

hplusN[\[Tau]_] := (4/D0) (G Mc/c^2)^(5/3) (\[Pi] fNFunc[\[Tau]]/c)^(2/3) Cos[2 \[Phi]NFunc[\[Tau]]];
hplus1PN[\[Tau]_] := (4/D0) (G Mc/c^2)^(5/3) (\[Pi] f1PNFunc[\[Tau]]/c)^(2/3) Cos[2 \[Phi]1Func[\[Tau]]];

tZoomStart = Max[0, tCommon - 0.15];

plotWaveformCompare = Plot[{hplusN[\[Tau]], hplus1PN[\[Tau]]},
   {\[Tau], tZoomStart, 0.999 tCommon},
   PlotPoints -> 4000, MaxRecursion -> 4,
   Frame -> True, FrameLabel -> {"t  [s]", "h+(t)"},
   PlotLabel -> "Restricted waveform, last ~0.15 s: Newtonian-only (blue) vs. 1PN (dashed red)",
   PlotStyle -> {Directive[RGBColor[0.20, 0.35, 0.65], Thickness[0.0015]],
     Directive[RGBColor[0.75, 0.15, 0.15], Dashed, Thickness[0.0015]]},
   PlotLegends -> {"Newtonian (0PN)", "1PN-corrected"},
   ImageSize -> 800];

plotWaveformCompare


(* ::Section:: *)
(*7. Where does the PN parameter x itself sit? (validity check)*)


Print["x at f = 20 Hz:  ", xOfF[20] // N];
Print["x at f = 100 Hz: ", xOfF[100] // N];
Print["x at ISCO (x=1/6 by convention): ", 1/6 // N,
  "  <-> f_ISCO = ", (c^3/(\[Pi] G M)) (1/6)^(3/2) // N, " Hz"];
Print["Interpretation: x is a few percent at 20-100 Hz (PN expansion is trustworthy there), ",
  "but grows to O(1/6) by the time the PN parameter reaches its formal ISCO value -- ",
  "the regime where the PN series itself is no longer trustworthy and full IMR/NR-tuned ",
  "waveform models are required."];



