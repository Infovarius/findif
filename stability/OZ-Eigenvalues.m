(*******************************************************************
This file was generated automatically by the Mathematica front end.
It contains Initialization cells from a Notebook file, which
typically will have the same name as this file except ending in
".nb" instead of ".m".

This file is intended to be loaded into the Mathematica kernel using
the package loading commands Get or Needs.  Doing so is equivalent
to using the Evaluate Initialization Cells menu command in the front
end.

DO NOT EDIT THIS FILE.  This entire file is regenerated
automatically each time the parent Notebook file is saved in the
Mathematica front end.  Any changes you make to this file will be
overwritten.
***********************************************************************)

<<ACPackages`

\!\(\(H = 2;\)\[IndentingNewLine]
  \(v[z_] = z \((H - z)\) 4/H\^2;\)\)

\!\(\(eqOZ[\[Omega]_, k_, R_] := \((v[z] - \[Omega]\/k)\) \((\(\(f'\)'\)[z] - \(k\^2\) f[z])\) - \(\(v'\)'\)[z] f[z] - \(\[ImaginaryI]\/\(k\ R\)\) \((D[f[z], {z, 4}] - 2 \( k\^2\) \(\(f'\)'\)[z] + \(k\^4\) f[z])\);\)\)











\!\(\(FiniteMatrix4Eigen[n_, sample_:  5, subst___] := Block[{h = H\/n, krseq, krs, unk, unkmain, unkfict, i, sm, s, order, lhs, rhs, \[Omega]}, \[IndentingNewLine]sm = Max[sample, 5]; \[IndentingNewLine]unk = Table[f\_\(i + 1/2\), {i, 0 - \(sm - 1\)\/2, n + \(sm - 1\)\/2 - 1}]; \[IndentingNewLine]unkmain = Table[f\_i, {i, 1/2, n - 1/2}]; \[IndentingNewLine]unkfict = Complement[unk, unkmain]; \[IndentingNewLine]dersam[order_, s_] := Table[f\_j, {j, i - \(s - 1\)\/2, i + \(s - 1\)\/2}] . NDCoefficientList[order, s]/h\^order; \[IndentingNewLine]rhs = \(\(-D[eqOZ[\[Omega], k, R], \[Omega]]\) /. {z \[Rule] i\ h, f[z_] \[Rule] f\_i, \(f''\)[z_] \[Rule] dersam[2, sm], \(f''''\)[z_] \[Rule] dersam[4, sm]}\) /. {subst}; \[IndentingNewLine]lhs = \(Simplify[\((eqOZ[\[Omega], k, R] - \[Omega]\ D[eqOZ[\[Omega], k, R], \[Omega]])\)] /. {z \[Rule] i\ h, f[z_] \[Rule] f\_i, \(f''\)[z_] \[Rule] dersam[2, sm], \(f''''\)[z_] \[Rule] dersam[4, sm]}\) /. {subst}; \[IndentingNewLine]virfict = First[Solve[Thread[{dersam[0, sm - 1] /. i \[Rule] 0, dersam[1, sm - 1] /. i \[Rule] 0, dersam[1, sm - 1] /. i \[Rule] n, dersam[0, sm - 1] /. i \[Rule] n} \[Equal] 0], unkfict]]; \[IndentingNewLine]fict = unkfict /. virfict; rhs = \(Coefficient[#, unkmain] &\) /@ \((Table[rhs, {i, 1/2, n - 1/2}] /. virfict)\); \[IndentingNewLine]lhs = \(Coefficient[#, unkmain] &\) /@ \((Table[lhs, {i, 1/2, n - 1/2}] /. virfict)\); \[IndentingNewLine] (*rhs = \(Coefficient[#, unk] &\) /@ Join[Table[0, {i, 2}, {j, n + sm - 1}], Table[rhs, {i, 1/2, n - 1/2}]/\[Omega], Table[0, {i, 2}, {j, n + sm - 1}]]; \[IndentingNewLine]lhs = \(Coefficient[#, unk] &\) /@ Join[{dersam[0, sm - 1] /. i \[Rule] 0, dersam[1, sm - 1] /. i \[Rule] 0}, Table[lhs, {i, 1/2, n - 1/2}], {dersam[1, sm - 1] /. i \[Rule] n, dersam[0, sm - 1] /. i \[Rule] n}];*) \[IndentingNewLine] (*\(Return[rhs . Inverse[lhs]];\)*) \[IndentingNewLine]Return[{lhs, rhs}];\[IndentingNewLine]];\)\)



GetIncrements[n_,k1_,R1_]:=
    Eigenvalues[
      Inverse[Last[#]].First[#]&[
        N[FiniteMatrix4Eigen[n,5,k\[Rule]k1,R\[Rule]R1]]]];

GetEigenWave[n1_,k1_,R1_,mode_:1]:=
    Block[{eigval,eigvec,nom,syst},
      syst=Inverse[Last[#]].First[#]&[
          N[FiniteMatrix4Eigen[n1,5,k\[Rule]k1,R\[Rule]R1]]];
      {eigval,eigvec}=
        Select[Transpose[Eigensystem[syst]],NumericQ[First[#]]&]//Transpose;
      nom=
        Flatten[Position[eigval,
              Sort[eigval,Im[#1]<Im[#2]&]\[LeftDoubleBracket]
                mode\[RightDoubleBracket]]]//First;
      Print["Main eigen value is ",eigval[[nom]]];
      Return[eigvec[[nom]]];];







GetMainIncrement[n_,k1_,R1_,
      num_:1]:=(*GetMainIncrement[n,k1,R1,num]=*)(Print[{n,k1,R1}];
      Take[Sort[Im/@GetIncrements[n,k1,R1]],{num}]//First);



<<"MainIncrement.m";

































































GetRecr[n2_,k2_,Rmin_,Rmax_]:=
  r/.FindRoot[Hold[GetMainIncrement[n2,k2,r,1]],{r,Rmin,Rmax}]















































































































































































































































