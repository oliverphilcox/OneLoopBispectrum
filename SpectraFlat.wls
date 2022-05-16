#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*In this script, the flattened (k1 = k2 + k3) one-loop bispectrum templates are computed using a given input power spectrum. This is designed to run on a HPC cluster and computes the templates for all biases simultaneously, but for a single angular component. *)
(**)
(*The script can be run via "./SpectraFlat.wls N", where N is the angular index (between 1 and 7). If the index is not specified, the script will run in testing mode, using a small value of the FFTLog index for speed. If the script is *not* run in testing mode, it will compute templates only if they do not already exist in the output directory*)


(* Read mu index from command line, or take a fixed value for testing *)
Check[muIndex=ToExpression[$ScriptCommandLine[[2]]];test=False;,muIndex=3;test=True;];
Print[ToString[StringForm["Angular index: ``",muIndex]]];


(* ::Subsection:: *)
(*Set-Up*)


(* ::Subsubsection::Closed:: *)
(*Options*)


(* ::Code::Initialization::Plain:: *)
(* Output string [prepended to file names] *)
string = "ptc"; 

(* ----------- Directories -----------  *)
If[test,
(* Directory containing power spectrum *)
pkdir = "/Users/ophilcox/Desktop/OneLoopBispectrum/pk";
(* Directory containing kernels *)
indir = "/Users/ophilcox/Desktop/OneLoopBispectrum/kernels";
(* Output directory [for testing] *)
outdir = "/Users/ophilcox/Desktop/OneLoopBispectrum/output";
,
(* Directory containing power spectrum *)
pkdir = "/home/ophilcox/OneLoopBispectrum/pk";
(* Directory containing kernels *)
indir = "/home/ophilcox/OneLoopBispectrum/kernels";
(* Output directory [should be large] *)
outdir = "/scratch/ophilcox/OneLoopBispectrum/output";
];

(* ----------- FFTLog Settings ------------- *)
(* Choices of FFTLog bias *)
(* NB: this allows two values of the B321-I bias *)
bias222 = -0.6;
bias321Iv1 = -0.3;
bias321Iv2 = -0.8;
bias321II = -0.6;
bias411 = -0.3;

(* Precision parameters *)
eps=10.^-8;
kmin=1 10^-5;
kmax=10.0;

(* Number of FFT frequencies *)
(* The default uses reduce = 1 *)
If[test, reduce=16, reduce=1];
Nmax222=64/reduce;
Nmax321I=96/reduce;
Nmax321II=128/reduce;
Nmax411=128/reduce;
Print["Nmax-222 = ",Nmax222];
Print["Nmax-321-I = ",Nmax321I];
Print["Nmax-321-II = ",Nmax321II];
Print["Nmax-411 = ",Nmax411];

(* Specification of the k-array. This has NmaxPlot logarithmically spaced values of k in [kminPlot, kmaxPlot] *)
NmaxPlot=64;
kminPlot=0.001;
kmaxPlot=0.5;

Clear[x,y,k];


(* ::Subsubsection::Closed:: *)
(*Load inputs*)


(* ::Code::Initialization::Plain:: *)
SetDirectory[outdir];
(* Check if output file exists *)
If[test,
includeB211=True;
includeB222=True;
includeB321Iv1=True;
includeB321Iv2=True;
includeB321II=True;
includeB411=True;
,
includeB211=Not[FileExistsQ[ToString[StringForm["``-flat_b211_m``.wdx",string,muIndex]]]];
includeB222=Not[FileExistsQ[ToString[StringForm["``-flat_b222_nu``_m``.wdx",string,bias222,muIndex]]]];
includeB321Iv1=Not[FileExistsQ[ToString[StringForm["``-flat_b321I_nu``_m``.wdx",string,bias321Iv1,muIndex]]]];
includeB321Iv2=Not[FileExistsQ[ToString[StringForm["``-flat_b321I_nu``_m``.wdx",string,bias321Iv2,muIndex]]]];
includeB321II=Not[FileExistsQ[ToString[StringForm["``-flat_b321II_nu``_m``.wdx",string,bias321II,muIndex]]]];
includeB411=Not[FileExistsQ[ToString[StringForm["``-flat_b411_nu``_m``.wdx",string,bias411,muIndex]]]];
];

(* Define angular indices *)
powsAng={};
For[i=0,i<=12,i++,
If[OddQ[i],Continue[]];
powsAng=Join[powsAng,{{i,0}}]];
\[Mu]index = powsAng[[muIndex,1]];
\[Chi]index = powsAng[[muIndex,2]];
Print[ToString[StringForm["Using \[Mu] index `` and \[Chi] index ``.",\[Mu]index,\[Chi]index]]];

(* Get angular coefficients *)
getCoeff[dat_]:=Module[{datLoc},
datLoc=dat;
If[\[Mu]index==0,datLoc=datLoc/.{\[Mu]->0},datLoc=Coefficient[datLoc,\[Mu]^\[Mu]index]];
datLoc];

Print["Computing B211: ",includeB211];
Print["Computing B222: ",includeB222];
Print["Computing B321I (v1): ",includeB321Iv1];
Print["Computing B321I (v2): ",includeB321Iv2];
Print["Computing B321II: ",includeB321II];
Print["Computing B411: ",includeB411];


(* ::Subsubsection::Closed:: *)
(*Define triangles*)


(* ::Code::Initialization::Plain:: *)
(* Define k binning *)
kBins[Nmax_,kmin_,kmax_]:=Module[{\[CapitalDelta],result},
\[CapitalDelta]=1/(Nmax-1) Log[kmax/kmin]//N;
result=Table[kmin Exp[(i-1) \[CapitalDelta]],{i,1,Nmax}];
result
];
knPlot=kBins[NmaxPlot,kminPlot,kmaxPlot];

(* Define array of x, y *)
xArr = Range[0.0,1.,0.025]^2;
yArr = Range[0.0,1.,0.025]^2;
xy1={};
For[i=1,i<=Length[xArr],i++,
For[j=1,j<=Length[yArr],j++,
x1=xArr[[i]];
y1=yArr[[j]];
(* NB omit flattened triangles *)
If[(Abs[Sqrt[x1]-Sqrt[y1]]<=1)&&(Sqrt[x1]+Sqrt[y1]==1)&&(y1>=x1),xy1=Join[xy1,{{x1,y1}}]]]];
ntriag=Length[xy1]
Print["Ntriag=",ntriag];
(* Replace first value to avoid errors *)
xy1[[1]]={N[10^-4],1};


(* ::Subsubsection::Closed:: *)
(*Kernel Definitions*)


(* ::Text:: *)
(*These are needed for B321-II computation*)


(* ::Code::Initialization::Plain:: *)
perm12=Permutations[{q1,q2}];
perm123=Permutations[{q1,q2,q3}];
perm1234=Permutations[{q1,q2,q3,q4}];

Clear[Fn];
Fn[n_,v_]:=If[n==1,1,Sum[Gn[m,v[[1;;m]]]/((2n+3)(n-1)) ((2n+1)al[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Fn[n-m,v[[m+1;;n]]]+2be[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Gn[n-m,v[[m+1;;n]]]),{m,1,n-1}]];
Gn[n_,v_]:=If[n==1,1,Sum[Gn[m,v[[1;;m]]]/((2n+3)(n-1)) (3al[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Fn[n-m,v[[m+1;;n]]]+2n be[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Gn[n-m,v[[m+1;;n]]]),{m,1,n-1}]];
F2s[q1_,q2_]=Simplify[Sum[1/Length[perm12]Fn[2,{perm12[[i,1]],perm12[[i,2]]}],{i,1,Length[perm12]}]];
F3s[q1_,q2_,q3_]=Simplify[Sum[1/Length[perm123]Fn[3,{perm123[[i,1]],perm123[[i,2]],perm123[[i,3]]}],{i,1,Length[perm123]}]];
F4s[q1_,q2_,q3_,q4_]=Sum[1/Length[perm1234]Fn[4,{perm1234[[i,1]],perm1234[[i,2]],perm1234[[i,3]],perm1234[[i,4]]}],{i,1,Length[perm1234]}];
(* alpha, beta kernels *)
alfm[k1_,k2_]:=1+1/2 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2)/mag[k1]^2;
befm[k1_,k2_]:=(mag[k1+k2]^2 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2))/(4mag[k1]^2 mag[k2]^2);

(* Angular operators *)
ang[qa_,qb_]:=(mag[qa+qb]^2-mag[qa]^2-mag[qb]^2)/(2mag[qa]mag[qb]) (* angle between qa and qb *)
angK[qa_,qb_]:=ang[qa,qb]^2-1 (* K kernel appearing in G2 *)
(* NB: ang1 will be replaced by ang later, it just allows removal of some nasty terms *)
angK1[qa_,qb_]:=ang1[qa,qb]^2-1 (* K kernel appearing in G2, unevaluated *)

(* RSD kernels *)
K1=b1;
Z1s[q1_]:=K1+ f \[Eta][q1]^2;
K2[q1_,q2_]:=(b2/2+g2 angK[q1,q2])+b1 Fn[2,{q1,q2}]
Z2[q1_,q2_]:=K2[q1,q2]+f \[Eta][q1+q2]^2 Gn[2,{q1,q2}]+(f \[Eta][q1+q2]mag[q1+q2])/2 (\[Eta][q1]/mag[q1]+\[Eta][q2]/mag[q2])K1+(f \[Eta][q1+q2]mag[q1+q2])^2/2 \[Eta][q1]/mag[q1] \[Eta][q2]/mag[q2]
Z2s[q1_,q2_]=Simplify[Sum[1/Length[perm12]Z2[perm12[[i,1]],perm12[[i,2]]],{i,1,Length[perm12]}]];
K3[q1_,q2_,q3_]:=(b3 /6+g2x angK[q1,q2]+g3 angL[q1,q2,q3]+g21 angK[q1,q2+q3]angK[q2,q3] )+(b2 Fn[2,{q1,q2}]+2g2 angK[q1,q2+q3]Gn[2,{q2,q3}] ) + (b1 Fn[3,{q1,q2,q3}])
Z3[q1_,q2_,q3_]:=K3[q1,q2,q3]+f \[Eta][q1+q2+q3]^2 Gn[3,{q1,q2,q3}]+(f mag[q1+q2+q3] \[Eta][q1+q2+q3])(\[Eta][q1+q2]/mag[q1+q2] Gn[2,{q1,q2}]+\[Eta][q3]/mag[q3] K2[q1,q2])+(f mag[q1+q2+q3] \[Eta][q1+q2+q3])^2/2 (2 \[Eta][q1+q2]/mag[q1+q2] \[Eta][q3]/mag[q3] Gn[2,{q1,q2}]+\[Eta][q1]/mag[q1] \[Eta][q2]/mag[q2] K1)+(f mag[q1+q2+q3] \[Eta][q1+q2+q3])^3/6 \[Eta][q1]/mag[q1] \[Eta][q2]/mag[q2] \[Eta][q3]/mag[q3]
Z3s[q1_,q2_,q3_]=Simplify[Sum[1/Length[perm123]Z3[perm123[[i,1]],perm123[[i,2]],perm123[[i,3]]],{i,1,Length[perm123]}]];


(* ::Subsubsection::Closed:: *)
(*Angular Substitutions*)


(* ::Code::Initialization::Plain:: *)
SetDirectory[indir];
muTable=Import["muTensorsFlat_n6.wdx"];
muSubs={};
For[n=6,n>=1,n--,
muSubs=Join[muSubs,{\[Eta][q]^n-> Sum[ang1[k1,q]^i1  muTable[[n,i1+1]],{i1,0,n}]}]];
muSubsExt={\[Eta][k1]->\[Mu]1,\[Eta][k2]->\[Mu]2,\[Eta][k3]->\[Mu]3,\[Eta][-k1-k2]->\[Mu]3,\[Eta][k1+k2]->-\[Mu]3,\[Eta][k2+k3]->-\[Mu]1,\[Eta][k1+k3]->-\[Mu]2,\[Eta][k1-q]:> k1/k1mq \[Mu]1-q/k1mq \[Eta][q],\[Eta][k1+q]:> k1/k1pq \[Mu]1+q/k1pq \[Eta][q],\[Eta][-k1-q]:> -(k1/k1pq)\[Mu]1-q/k1pq \[Eta][q],\[Eta][k2-q]:> k2/k2mq \[Mu]2-q/k2mq \[Eta][q],\[Eta][-k2+q]:> -(k2/k2mq)\[Mu]2+q/k2mq \[Eta][q],\[Eta][k2+q]:> k2/k2pq \[Mu]2+q/k2pq \[Eta][q],\[Eta][k3-q]:> k3/k3mq \[Mu]3-q/k3mq \[Eta][q],\[Eta][k3+q]:> k3/k3pq \[Mu]3+q/k3pq \[Eta][q],\[Eta][-k1-k2+q]:> k3/k3pq \[Mu]3+q/k3pq \[Eta][q],\[Eta][k1+k2+q]:>- (k3/k3mq)\[Mu]3+q/k3mq \[Eta][q],\[Eta][k1+k2-q]:>- (k3/k3pq)\[Mu]3-q/k3pq \[Eta][q]};


(* ::Subsubsection::Closed:: *)
(*Load Tables*)


(* ::Text:: *)
(*Here we load the FFTLog kernels computed using the Kernels.nb file*)


(* ::Code::Initialization::Plain:: *)
(* B211 *)
b1s={b1,f};
b2s={b2,g2,b1,f, b1 f, f^2};

(* Compute B211 and correct mu coefficient *)
If[includeB211,
biasListB211=Union[Flatten[Table[b2s[[i]]b1s[[j]]b1s[[k]],{i,1,Length[b2s]},{j,1,Length[b1s]},{k,1,Length[b1s]}]]];
Print[ToString[StringForm["`` contributions to B211",Length[biasListB211]]]];
B211tmp[k1_,k2_,k3_,\[Mu]1_,\[Mu]2_,\[Mu]3_]=(2(Z1s[k1]Z1s[k2]Z2s[k1,k2]Plin[k1]Plin[k2]+Z1s[k1]Z1s[k3]Z2s[k1,k3]Plin[k1]Plin[k3]+Z1s[k2]Z1s[k3]Z2s[k2,k3]Plin[k2]Plin[k3])/.{al->alfm,be->befm})/.muSubsExt/.{mag[k1]->k1,mag[k2]->k2,mag[k3]->k3,mag[k1+k2]->k3,mag[k1+k3]->k2,mag[k2+k3]->k1}//Simplify;
B211[k1_,k2_,k3_]=getCoeff[B211tmp[k1,k2,k3,\[Mu],-\[Mu],-\[Mu]]];];


(* ::Code::Initialization::Plain:: *)
(* B222 *)
(* Load table and assert flattened angular structure *)
If[includeB222,
biasListB222init=Union[Flatten[Table[b2s[[i]]b2s[[j]]b2s[[k]],{i,1,Length[b2s]},{j,1,Length[b2s]},{k,1,Length[b2s]}]]];
SetDirectory[indir];
Tab222derivsinit=getCoeff[Import["tab222derivs-flat.wdx"]/.{\[Mu]1->\[Mu],\[Mu]2->-\[Mu],\[Mu]3->-\[Mu]}];
Tab222freq=Import["tab222freq-flat.wdx"];

(* Remove any unneeded bias parameters *)
biasListB222cut = biasListB222init;
Tab222derivscut = Tab222derivsinit;
For[i=1,i<=Length[biasListB222init],i++,s=Sum[Abs[Tab222derivsinit[[i,j]]],{j,1,Length[Tab222derivsinit[[1]]]}];If[s==0,biasListB222cut[[i]]=0;]];
biasListB222 = DeleteCases[biasListB222cut,0];
Print[ToString[StringForm["`` contributions to B222",Length[biasListB222]]]];
Tab222derivs = DeleteCases[Tab222derivscut,ConstantArray[0,Length[Tab222freq]]];
];


(* ::Code::Initialization::Plain:: *)
(* B321-I *)
b3sI={b3,g2x, g3, g21, b2,g2,b1,f^3,b1 f^2,f^2,b2 f, g2 f, b1 f,f};
If[Or[includeB321Iv1,includeB321Iv2],
SetDirectory[indir];
Print["Loading B321-I"];
biasListB321Iinit=Union[Flatten[Table[b3sI[[i]]b2s[[j]]b1s[[k]],{i,1,Length[b3sI]},{j,1,Length[b2s]},{k,1,Length[b1s]}]]];
Tab321derivs=Import["tab321derivs-flat.wdx"];
Tab321freq=Import["tab321freq-flat.wdx"];
Tab321derivs123=getCoeff[Tab321derivs/.{x->x,y->y,\[Mu]1->\[Mu]1,\[Mu]2->\[Mu]2,\[Mu]3->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab321derivs132=getCoeff[Tab321derivs/.{x->y,y->x,\[Mu]1->\[Mu]1,\[Mu]3->\[Mu]2,\[Mu]2->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab321derivs213=getCoeff[Tab321derivs/.{x->x/y,y->1/y,\[Mu]1->\[Mu]2,\[Mu]2->\[Mu]1,\[Mu]3->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab321derivs231=getCoeff[Tab321derivs/.{x->1/y,y->x/y,\[Mu]1->\[Mu]2,\[Mu]3->\[Mu]1,\[Mu]2->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab321derivs321=getCoeff[Tab321derivs/.{x->1/x,y->y/x,\[Mu]2->\[Mu]2,\[Mu]3->\[Mu]1,\[Mu]1->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab321derivs312=getCoeff[Tab321derivs/.{x->y/x,y->1/x,\[Mu]2->\[Mu]1,\[Mu]3->\[Mu]2,\[Mu]1->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Remove[Tab321derivs];

(* Remove any unneeded bias parameters *)
biasListB321Icut = biasListB321Iinit;
For[i=1,i<=Length[biasListB321Iinit],i++,s=Sum[Abs[Tab321derivs123[[i,j]]]+Abs[Tab321derivs132[[i,j]]]+Abs[Tab321derivs213[[i,j]]]+Abs[Tab321derivs213[[i,j]]]+Abs[Tab321derivs321[[i,j]]]+Abs[Tab321derivs312[[i,j]]],{j,1,Length[Tab321derivs123[[1]]]}];If[s==0,biasListB321Icut[[i]]=0;]];
biasListB321I = DeleteCases[biasListB321Icut,0];
Print[ToString[StringForm["`` contributions to B321-I",Length[biasListB321I]]]];
Tab321derivs123 = DeleteCases[Table[If[biasListB321Icut[[i]]==0,{},Tab321derivs123[[i,All]],Tab321derivs123[[i,All]]],{i,1,Length[biasListB321Iinit]}],{}];
Tab321derivs132 = DeleteCases[Table[If[biasListB321Icut[[i]]==0,{},Tab321derivs132[[i,All]],Tab321derivs132[[i,All]]],{i,1,Length[biasListB321Iinit]}],{}];
Tab321derivs213 = DeleteCases[Table[If[biasListB321Icut[[i]]==0,{},Tab321derivs213[[i,All]],Tab321derivs213[[i,All]]],{i,1,Length[biasListB321Iinit]}],{}];
Tab321derivs231 = DeleteCases[Table[If[biasListB321Icut[[i]]==0,{},Tab321derivs231[[i,All]],Tab321derivs231[[i,All]]],{i,1,Length[biasListB321Iinit]}],{}];
Tab321derivs312 = DeleteCases[Table[If[biasListB321Icut[[i]]==0,{},Tab321derivs312[[i,All]],Tab321derivs312[[i,All]]],{i,1,Length[biasListB321Iinit]}],{}];
Tab321derivs321 = DeleteCases[Table[If[biasListB321Icut[[i]]==0,{},Tab321derivs321[[i,All]],Tab321derivs321[[i,All]]],{i,1,Length[biasListB321Iinit]}],{}];
];


(* ::Code::Initialization::Plain:: *)
(* B321-II *)
b3sII={b1,g2,g21,f, f^2,b1 f^2,f^3}; 
If[includeB321II,
biasListB321II=Union[Flatten[Table[b3sII[[i]]b2s[[j]]b1s[[k]],{i,1,Length[b3sII]},{j,1,Length[b2s]},{k,1,Length[b1s]}]]];
Print[ToString[StringForm["`` contributions to B321-II",Length[biasListB321II]]]];
];


(* ::Code::Initialization::Plain:: *)
(* B411 *)
b4s={g21x, g211, g22, g31, g21,g3,b3, g2x, g2,b2,b1,f b3, f g2x, f g3, f g21,f b2, f g2, f b1, f,  f^2 b2, f^2 g2, f^2 b1, f^2,f^3 b1, f^3, f^4};
If[includeB411,
SetDirectory[indir];
biasListB411init=Union[Flatten[Table[b4s[[i]]b1s[[j]]b1s[[k]],{i,1,Length[b4s]},{j,1,Length[b1s]},{k,1,Length[b1s]}]]];
Tab411derivs=Import["tab411derivs-flat.wdx"];
Tab411freq=Import["tab411freq-flat.wdx"];

Tab411derivs123=getCoeff[Tab411derivs/.{x->x,y->y,\[Mu]1->\[Mu]1,\[Mu]2->\[Mu]2,\[Mu]3->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab411derivs231=getCoeff[Tab411derivs/.{x->1/y,y->x/y,\[Mu]1->\[Mu]2,\[Mu]3->\[Mu]1,\[Mu]2->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
Tab411derivs132=getCoeff[Tab411derivs/.{x->y,y->x,\[Mu]1->\[Mu]1,\[Mu]3->\[Mu]2,\[Mu]2->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];

(* Remove any unneeded bias parameters *)
biasListB411cut = biasListB411init;
For[i=1,i<=Length[biasListB411init],i++,s=Sum[Abs[Tab411derivs123[[i,j]]]+Abs[Tab411derivs231[[i,j]]]+Abs[Tab411derivs132[[i,j]]],{j,1,Length[Tab411derivs[[1]]]}];If[s==0,biasListB411cut[[i]]=0;]];
biasListB411 = DeleteCases[biasListB411cut,0];
Tab411derivs123 = DeleteCases[Table[If[biasListB411cut[[i]]==0,{},Tab411derivs123[[i,All]],Tab411derivs123[[i,All]]],{i,1,Length[biasListB411init]}],{}];
Tab411derivs231 = DeleteCases[Table[If[biasListB411cut[[i]]==0,{},Tab411derivs231[[i,All]],Tab411derivs231[[i,All]]],{i,1,Length[biasListB411init]}],{}];
Tab411derivs132 = DeleteCases[Table[If[biasListB411cut[[i]]==0,{},Tab411derivs132[[i,All]],Tab411derivs132[[i,All]]],{i,1,Length[biasListB411init]}],{}];
Print[ToString[StringForm["`` contributions to B411",Length[biasListB411]]]];

Remove[Tab411derivs];
];


(* ::Subsubsection::Closed:: *)
(*Linear Power Spectrum*)


(* ::Text:: *)
(*Here we load in the power spectrum file - this should be an IR-resummed spectrum.*)
(*The input script may need to be modified depending on the Pk set-up in use.*)


(* ::Code::Initialization::Plain:: *)
Print["Loading power spectrum"];
SetDirectory[pkdir];
InputTab=Import["ptc_pk_IR.dat"][[5;;All]]; (* This is IR resummed PT challenge P(k) *)
plindat=Table[{InputTab[[i,1]],-InputTab[[i,3]]/InputTab[[i,1]]^2},{i,1,Length[InputTab]}];

anorm=plindat[[1,2]];
kb=plindat[[1,1]];
pextralow=Table[{10^(-10+( i-1) 0.2), anorm (10^(-10+ (i-1) 0.2)/kb)^0.965},{i,1,25}];
pextrahighk=Table[{18+10^(( i) 0.1), 0.00015((18+10^(( i) 0.1))/(18+10^(( 20) 0.1)))^-3},{i,20,30}];
plindat2=Join[pextralow,plindat,pextrahighk];
plinus=Interpolation[plindat2];
plin=Interpolation[Log[plindat2]];
klow=10^-8; khigh=100.;
pcut[k_]= Exp[plin[Log[k]]]Exp[-(k/khigh)^2];
Plin[k_]=pcut[k];


(* ::Subsubsection::Closed:: *)
(*UV Limits*)


(* ::Text:: *)
(*This loads the UV limits, which are non-trivial only for B321-II and B411.*)


(* ::Code::Initialization::Plain:: *)
SetDirectory[indir];
Print["Loading UV limits"];
If[includeB222,
B222UVkernel = 0;
(* B222/k^2 UV limit [NB: zero here!] *)
lim1=1/(2\[Pi]^2) NIntegrate[Plin[q]^3,{q,0,Infinity},WorkingPrecision->10,MaxRecursion->1000,PrecisionGoal->5];
B222UV[x_,y_]=getCoeff[(B222UVkernel lim1)/.{k1->1,k2->y^(1/2),k3->x^(1/2)}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
];
If[Or[includeB321Iv1,includeB321Iv2],
B321UVkernel = 0;
(* B321-I/k^2 UV limit [NB: zero here!] *)
lim2=1/(2\[Pi]^2) NIntegrate[Plin[q]^2,{q,0,Infinity},WorkingPrecision->10,MaxRecursion->1000,PrecisionGoal->5];
If[B321UVkernel==0,B321UV[x_,y_]=getCoeff[B321UVkernel lim2/.{k1->1,k2->y^(1/2),k3->x^(1/2)}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}],Print["Need to add B321 UV symmetries!"];Abort[];];
];
If[includeB321II,
(* Explicit B321-II form, given P13 *)
z2[k1_,k2_,k3_,\[Mu]1_,\[Mu]2_,\[Mu]3_]:=(Z2s[k1,k2]/.{al->alfm,be->befm}/.{\[Eta][k1]->\[Mu]1,\[Eta][k2]->\[Mu]2,\[Eta][k1+k2]->-\[Mu]3})/.{mag[k1]->k1,mag[k2]->k2,mag[k1+k2]->k3}//Simplify;
z1[kT_,\[Mu]T_]:=Z1s[kT]/.{\[Eta][kT]->\[Mu]T};
B321II[k1_,k2_,k3_]:=(z2[k1,k2,k3,\[Mu]1,\[Mu]2,\[Mu]3]z1[k1,\[Mu]1]Plin[k1]P13[k2,\[Mu]2]+z2[k1,k3,k2,\[Mu]1,\[Mu]3,\[Mu]2]z1[k1,\[Mu]1]Plin[k1]P13[k3,\[Mu]3]+z2[k2,k1,k3,\[Mu]2,\[Mu]1,\[Mu]3]z1[k2,\[Mu]2]Plin[k2]P13[k1,\[Mu]1]+z2[k2,k3,k1,\[Mu]2,\[Mu]3,\[Mu]1]z1[k2,\[Mu]2]Plin[k2]P13[k3,\[Mu]3]+z2[k3,k1,k2,\[Mu]3,\[Mu]1,\[Mu]2]z1[k3,\[Mu]3]Plin[k3]P13[k1,\[Mu]1]+z2[k3,k2,k1,\[Mu]3,\[Mu]2,\[Mu]1]z1[k3,\[Mu]3]Plin[k3]P13[k2,\[Mu]2]);(* P13 FFTlog kernel and table [NB: does not include a Z1 factor!] *)
kerP13=6(Z3s[kP,q,-q]/.{g2x->0,g3->0,b3->0})/.muSubsExt/.{ang1->ang,al->alfm,be->befm,\[Eta][0]->0}/.{mag[kP]->kP,mag[q]->q,mag[-q]->q,mag[kP-q]->kPmq,mag[0]->0,mag[kP+q]->kPpq,\[Eta][kP]->\[Mu]P,\[Eta][-q]->-\[Eta][q],\[Eta][kP-q]->kP/kPmq \[Mu]P-q/kPmq \[Eta][q],\[Eta][kP+q]->kP/kPmq \[Mu]P+q/kPmq \[Eta][q]}//Expand;
kerP13=((((kerP13/.{\[Eta][q]^2-> 1/2 (1-\[Mu]P^2)+ang[q,kP]^2/2 (3 \[Mu]P^2-1),\[Eta][q]->ang[q,kP]\[Mu]P})/.{mag[kP+q]->kPpq,mag[kP]->kP,mag[q]->q})//Expand)/.{kPmq^2 kPpq^2-> kPmq^2 (2 kP^2+2q^2-kPmq^2),kPmq^4 kPpq^2-> kPmq^4 (2 kP^2+2q^2-kPmq^2),kPmq^2 kPpq^4-> kPpq^4 (2 kP^2+2q^2-kPpq^2),kPpq^2/kPmq^2-> (2 kP^2+2q^2-kPmq^2)/kPmq^2,kPpq^4/kPmq^2-> (2 kP^2+2q^2-kPmq^2)^2/kPmq^2,kPmq^2/kPpq^2-> (2 kP^2+2q^2-kPpq^2)/kPpq^2}//Expand)/.{kPpq->kPmq};
J2[k_,\[Nu]1_,\[Nu]2_]:=1/(2\[Pi])^3 (Gamma[3/2-\[Nu]1]Gamma[3/2-\[Nu]2]Gamma[\[Nu]1+\[Nu]2-3/2])/(Gamma[\[Nu]1]Gamma[\[Nu]2]Gamma[3-\[Nu]1-\[Nu]2]) \[Pi]^(3/2) k^(3-2\[Nu]1-2\[Nu]2);
Tab13=Table[{-(1/2)q D[kerP13[[i]] ,q]/kerP13[[i]],-(1/2)kPmq D[kerP13[[i]] ,kPmq]/kerP13[[i]],kerP13[[i]]/.{q->1,kPmq->1,kP->1}},{i,1,Length[kerP13]}];
P13ab[n1_]=Sum[Tab13[[i,3]]J2[1,Tab13[[i,1]]+n1,Tab13[[i,2]]],{i,1,Length[Tab13]}]//FullSimplify;
(* P13 UV limits (assuming -1 < \[Nu] < 1 *)
\[Alpha]UV13=1/315 (-61 b1-576 g2+672 g21-3 f (-9+48 f+35 b1 (2+f)) \[Mu]P^2-3 f^2 (46+35 f) \[Mu]P^4) 1/(2 \[Pi]^2) NIntegrate[Plin[q],{q,0,Infinity}];
];
If[includeB411,
B411UVkernel = Import["b411uv-flat.wdx"];
(* B411/k^2 UV limit *)
\[Sigma]v=1/(6\[Pi]^2) NIntegrate[Plin[q],{q,0,Infinity},WorkingPrecision->10,MaxRecursion->1000,PrecisionGoal->5];
B411UV123[x_,y_]=getCoeff[(B411UVkernel(3 \[Sigma]v)/.{k1->1,k2->y^(1/2),k3->x^(1/2)})/.{x->x,y->y,\[Mu]1->\[Mu]1,\[Mu]2->\[Mu]2,\[Mu]3->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
B411UV231[x_,y_]=getCoeff[(B411UVkernel(3 \[Sigma]v)/.{k1->1,k2->y^(1/2),k3->x^(1/2)})/.{x->1/y,y->x/y,\[Mu]1->\[Mu]2,\[Mu]3->\[Mu]1,\[Mu]2->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
B411UV132[x_,y_]=getCoeff[(B411UVkernel(3 \[Sigma]v)/.{k1->1,k2->y^(1/2),k3->x^(1/2)})/.{x->y,y->x,\[Mu]1->\[Mu]1,\[Mu]3->\[Mu]2,\[Mu]2->\[Mu]3}/.{\[Mu]2->-\[Mu],\[Mu]3->-\[Mu],\[Mu]1->\[Mu]}];
];


(* ::Subsubsection::Closed:: *)
(*FFTLog Functions*)


(* ::Text:: *)
(*Define the FFTLog transformations, as well as the Gamma matrices and Hypergeometric functions. This is based on Marko Simonovic's FFTLog code, and vectorizes over all x values.*)


(* ::Code::Initialization::Plain:: *)
Print["Compiling FFTlog functions"];
CoeffPow[b_,Nmax_,kmin_,kmax_,Pk_]:=Module[{\[CapitalDelta],n,m,kn,\[Eta]m,\[Eta]n,Pn,cn,cnsym,result},
\[CapitalDelta]=1/(Nmax-1) Log[kmax/kmin]//N;
kn=Table[kmin Exp[(i-1) \[CapitalDelta]],{i,1,Nmax}];
\[Eta]m=Table[b+(2\[Pi] I)/(Nmax \[CapitalDelta]) (j-Nmax/2-1),{j,1,Nmax+1}];

Pn=Table[{Pk[kn[[i]]]Exp[-b(i-1)\[CapitalDelta]]},{i,1,Nmax}];
cn=Fourier[Pn,FourierParameters->{-1,-1}];
cnsym=Table[If[(j-Nmax/2)<1,kmin^-\[Eta]m[[j]] Conjugate[cn[[-j+Nmax/2+2,1]]],kmin^-\[Eta]m[[j]] cn[[j-Nmax/2,1]]],{j,1,Nmax+1}];
cnsym[[1]]=cnsym[[1]]/2;
cnsym[[Length[cnsym]]]=cnsym[[Length[cnsym]]]/2;
result=Table[{cnsym[[i]],\[Eta]m[[i]]},{i,1,Length[\[Eta]m]}];
result
];
GammaC=Compile[{{z,_Complex}},
Block[{q0,q1,q2,q3,q4,q5,q6,p1,p2,result},
q0=75122.6331530+0.0I;
	q1=80916.6278952+0.0I;
	q2=36308.2951477+0.0I;
	q3=8687.24529705+0.0I;
	q4=1168.92649479+0.0I;
	q5=83.8676043424+0.0I;
	q6=2.50662827511+0.0I;
If[Re[z]>=0,p1=(q0+q1 z+q2 z^2+q3 z^3+q4 z^4+q5 z^5+q6 z^6)/(z(z+1)(z+2)(z+3)(z+4)(z+5)(z+6));result=p1 (z+5.5)^(z+0.5) Exp[-z-5.5],p1=(q0+q1(1-z)+q2 (1-z)^2+q3 (1-z)^3+q4 (1-z)^4+q5 (1-z)^5+q6 (1-z)^6)/((1-z)(2-z)(3-z)(4-z)(5-z)(6-z)(7-z));p2=p1 (1-z+5.5)^(1-z+0.5) Exp[-1+z-5.5];result=\[Pi]/(Sin[\[Pi] z]p2)];
result
],
{{Sin[_],_Complex},{Exp[_],_Complex},{Re[_],_Real}}, 
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];
Hyp2F1C=Compile[{{a,_Complex},{b,_Complex},{c,_Complex},{z,_Complex}},
Block[{p,s,sold,eps,n},
s=0.0+0.0 I;
p=1.0+0.0I;
eps=1.0;
n=0.0;
While[eps>10.^-10,
sold=s;
s=s+p;
p=p ((a+n)(b+n))/((c+n)(n+1)) z;
eps=Abs[(s-sold)/s];
n++];
s
],
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];

JCx = Compile[{{\[Nu]1,_Complex},{\[Nu]2,_Complex},{\[Nu]3,_Complex},{xtab,_Complex,1},{y,_Complex},{epscr,_Real}},
Block[
{p1,p2,s1,s2,s1old,s2old,eps,n,result,unittab,jmax},
jmax=Length[xtab]-1;
unittab=Table[1.+I 0.,{j,0,jmax}];

s1=(0.+I 0.)xtab;
s2=(0.+I 0.)xtab;

p1=unittab (GammaC[2-2 \[Nu]2] GammaC[2-2 \[Nu]3] Sin[\[Pi] \[Nu]2] Sin[\[Pi] \[Nu]3])/GammaC[5-2 \[Nu]2-2 \[Nu]3];
p2=unittab (GammaC[2-2 \[Nu]1] GammaC[2 (-2+\[Nu]1+\[Nu]2+\[Nu]3)] Sin[\[Pi] \[Nu]1] Sin[\[Pi] (\[Nu]1+\[Nu]2+\[Nu]3)])/GammaC[-1+2 \[Nu]2+2 \[Nu]3];

eps=1.;
n=0;
While[eps>epscr,
s1old=s1;
s2old=s2;
s1=s1+p1 Hyp2F1C[\[Nu]1+n,1.5-\[Nu]2+n,3-\[Nu]2-\[Nu]3+2 n,1-y];
s2=s2+p2 Hyp2F1C[\[Nu]2+n,1.5-\[Nu]1+n,\[Nu]2+\[Nu]3+2 n,1-y];
p1=p1*((n+\[Nu]1) (3-\[Nu]1-\[Nu]2-\[Nu]3+n) (1.5-\[Nu]3+n) (1.5-\[Nu]2+n))/((1+n) (2.5-\[Nu]2-\[Nu]3+n) (4-\[Nu]2-\[Nu]3+2 n) (3-\[Nu]2-\[Nu]3+2 n)) xtab;
p2=p2*((1.5-\[Nu]1+n) (\[Nu]1+\[Nu]2+\[Nu]3-1.5+n) (n+\[Nu]2) (n+\[Nu]3))/((1+n) (\[Nu]2+\[Nu]3-0.5+n) (1+\[Nu]2+\[Nu]3+2 n) (\[Nu]2+\[Nu]3+2 n)) xtab;

eps=Abs[(s1[[jmax+1]]-s2[[jmax+1]]-s1old[[jmax+1]]+s2old[[jmax+1]])/(s1[[jmax+1]]-s2[[jmax+1]])];
n++];

result=Sec[\[Pi] (\[Nu]2+\[Nu]3)]/(2 \[Pi]^2) (xtab^(1.5-\[Nu]2-\[Nu]3)s1- y^(1.5-\[Nu]1-\[Nu]3)s2);
result
],
{{Hyp2F1C[_],_Complex},{GammaC[_],_Complex},{Sec[_],_Complex},{Abs[_],_Real},{Sin[_],_Complex}}
,
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];


(* ::Subsection:: *)
(*Main Computation*)


(* ::Subsubsection::Closed:: *)
(*Compute k, x, y grid*)


(* ::Text:: *)
(*Output an array containing all k, x, y values used in the computation.*)


(* ::Code::Initialization::Plain:: *)
Print["Computing k, x, y bins"];
time=AbsoluteTiming[
	xyout=xy1;
	][[1]];
Print["Runtime: ",time];

SetDirectory[outdir];
kcoords=Join[{knPlot},{xyout[[All,1]]},{xyout[[All,2]]}];
Export[ToString[StringForm["``-flat_rsd_k_x_y.wdx",string]],kcoords];


(* ::Subsubsection::Closed:: *)
(*Compute B211 spectra*)


(* ::Text:: *)
(*Compute the tree-level bispectrum contributions. These are included for reference only: it is much faster to compute them directly in Python!*)


computeB211[x1_,y1_]:=Module[{Btree},
Btree=Table[B211[knPlot[[i]],Sqrt[y1]knPlot[[i]],Sqrt[x1]knPlot[[i]]],{i,1,Length[knPlot]}];
B211plot = Table[Coefficient[Btree[[All]],biasListB211[[j]]]/.{f->0,b1->0,b2->0,g2->0},{j,1,Length[biasListB211]}];
Transpose[B211plot]
];


(* ::Code::Initialization::Plain:: *)
If[includeB211,
Print["Computing B211 spectra"];
time=AbsoluteTiming[
	B211spec=Table[computeB211[xyout[[i,1]],xyout[[i,2]]],{i,1,Length[xyout]}]
][[1]];
Print["Runtime: ",time];
];

SetDirectory[outdir];
If[includeB211,Export[ToString[StringForm["``-flat_b211_m``.wdx",string,muIndex]],B211spec]];


(* ::Subsubsection::Closed:: *)
(*Compute B222 spectra*)


(* ::Text:: *)
(*Here, we compute the B222 loop integrals via FFTLog. This is usually the most expensive part of the calculation (for large Subscript[N, max]), since it is cubic in the number of FFTLog points.*)


(* ::Code::Initialization::Plain:: *)
(* Define UV derivatives *)
B222UVderivs={};
For[i=1,i<=Length[biasListB222],i++,
biasFunc=Coefficient[B222UV[x,y],biasListB222[[i]]]/.{f->0,b1->0,b2->0,g2->0};
B222UVderivs=Join[B222UVderivs,{biasFunc}];
];

(* Function to compute B222 *)
computeB222[x1vec_,y1_]:=Module[{Bn1n2n3step1,b222,coefTab,coefTab1,cnc,Nmax,bias,cn,mu3,tab222,b222plot},
(* Compute FFTLog coefficients *)
Nmax=Nmax222;
bias=bias222;
cn=CoeffPow[bias,Nmax,kmin,kmax,Plin];
(* Compute FFTLog table with relevant angles*)
tab222=Tab222derivs/.{y->y1};
coefTab=Table[tab222[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB222]},{k,1,Length[Tab222freq]}];
(* Sum over coefficients *)
Bn1n2n3step1=ParallelTable[Sum[coefTab[[All,All,i]]JCx[Tab222freq[[i,1]]/.{\[Nu]1->cn[[i1,2]],\[Nu]2->cn[[i2,2]],\[Nu]3->cn[[i3,2]]},Tab222freq[[i,2]]/.{\[Nu]1->cn[[i1,2]],\[Nu]2->cn[[i2,2]],\[Nu]3->cn[[i3,2]]},Tab222freq[[i,3]]/.{\[Nu]1->cn[[i1,2]],\[Nu]2->cn[[i2,2]],\[Nu]3->cn[[i3,2]]},x1vec,y1,eps],{i,1,Length[Tab222freq]}],{i1,1,Nmax/2+1},{i2,1,Nmax+1},{i3,1,Nmax+1}];
b222=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,i3,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,Nmax+2-i3,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{i3,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB222]}];
b222plot=ConstantArray[0,{NmaxPlot,Length[x1vec],Length[biasListB222]}];
For[i=1,i<=Length[b222plot],i++,
cnc=knPlot[[i]]^(1+cn[[All,2]]) cn[[All,1]];
For[n=1,n<=Length[x1vec],n++,
For[j=1,j<=Length[biasListB222],j++,
b222plot[[i,n,j]]=(Re[((b222[[All,All,All,n,j]].cnc).cnc).cnc+(knPlot[[i]])^2 B222UVderivs[[j]]/.{x->x1vec[[n]],y->y1}])]]];
b222plotOut=ConstantArray[0,{Length[b222plot],Length[b222plot[[1]]],Length[biasListB222init]}];
j=1;
For[i=1,i<=Length[biasListB222init],i++,If[biasListB222cut[[i]]==0,Continue[],Continue[],b222plotOut[[All,All,i]]=b222plot[[All,All,j]];j+=1]];
b222plotOut];


If[includeB222,
Print["Computing B222 spectra"];
time=AbsoluteTiming[
	B222spec={};
	For[m=1,m<=Length[xyout],m++,
		(* Define relevant x values (NB: should be ordered!) *)
		b222=computeB222[{xyout[[m,1]]},xyout[[m,2]]];
	B222spec=Join[B222spec,{b222[[All,1,All]]}]];
	][[1]];
Print["Runtime: ",time];
];

SetDirectory[outdir];
If[includeB222,Export[ToString[StringForm["``-flat_b222_nu``_m``.wdx",string,bias222,muIndex]],B222spec]];


(* ::Subsubsection::Closed:: *)
(* Compute B321-I spectra*)


(* ::Text:: *)
(*Repeat for B321-I. This is quadratic in the number of FFTLog points, but also has a sixfold sum over permutations. We compute spectra for two choices of FFTLog frequency here.*)


(* ::Code::Initialization::Plain:: *)
(* Define UV derivatives *)
B321UVderivs={};
For[i=1,i<=Length[biasListB321I],i++,
biasFunc=Coefficient[B321UV[x,y],biasListB321I[[i]]]/.{f->0,b1->0,b2->0,g2->0,b3->0,g21->0,g2x->0,g3->0};
B321UVderivs=Join[B321UVderivs,{biasFunc}];
]

(* Function to compute B321-I *)
computeB321I[x1vec_,y1_,bias_]:=Module[{Bn1n2n3step1,b321,b321p1,b321p2,b321p3,b321p4,b321p5,b321p6,b321plot,coefTab,coefTabExpan,coefTabSafe,cnc,Nmax,cn,nu23,nu13,nu12,rho12,rho23,rho13,tab321},
(* Compute FFTlog coefficients *)
Nmax=Nmax321I;
cn=CoeffPow[bias,Nmax,kmin,kmax,Plin];

(* Create output array *)
b321plot=ConstantArray[0,{NmaxPlot,Length[x1vec],Length[biasListB321I]}];

(* Permutation 1: (x,y)\[Rule](x,y) *)
coefTab2=Tab321derivs123/.{y->y1};
coefTab=Table[coefTab2[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB321I]},{k,1,Length[Tab321freq]}];
Bn1n2n3step1=ParallelTable[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[Tab321freq]\)]\(coefTab[\([All, All, i]\)]\ JCx[eps + Tab321freq[\([i, 1]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 2]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, x1vec, y1, eps]\)\),{i1,1,Nmax/2+1},{i2,1,Nmax+1}];
b321p1=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB321I]}];
For[i=1,i<=Length[b321plot],i++,cnc=knPlot[[i]]^(1.5` +cn[[All,2]]) cn[[All,1]];
For[n=1,n<=Length[x1vec],n++,
For[j=1,j<=Length[biasListB321I],j++,b321plot[[i,n,j]]+=Plin[knPlot[[i]]]Re[(b321p1[[All,All,n,j]].cnc.cnc+knPlot[[i]]^2B321UVderivs[[j]]/.{x->x1vec[[n]],y->y1})]]]];

(* Permutation 2: (x,y)\[Rule](y,x) *)
(* Note that we have already imposed permutations in the derivatives here!*)
coefTab2=Tab321derivs132/.{y->y1};
coefTab=Table[coefTab2[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB321I]},{k,1,Length[Tab321freq]}];
Bn1n2n3step1=ParallelTable[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[Tab321freq]\)]\(coefTab[\([All, All, i]\)]\ JCx[eps + Tab321freq[\([i, 2]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 1]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, x1vec, y1, eps]\)\),{i1,1,Nmax/2+1},{i2,1,Nmax+1}];
b321p2=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB321I]}];
For[i=1,i<=Length[b321plot],i++,cnc=knPlot[[i]]^(1.5` +cn[[All,2]]) cn[[All,1]];
For[n=1,n<=Length[x1vec],n++,
For[j=1,j<=Length[biasListB321I],j++,b321plot[[i,n,j]]+=Plin[knPlot[[i]]]Re[(b321p2[[All,All,n,j]].cnc.cnc+knPlot[[i]]^2B321UVderivs[[j]]/.{x->y1,y->x1vec[[n]]})]]]];

(* Permutation 3: (x,y)\[Rule](x/y,1/y) *)
coefTab2=Tab321derivs213/.{y->y1};
coefTab=Table[coefTab2[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB321I]},{k,1,Length[Tab321freq]}];
Bn1n2n3step1=ParallelTable[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[Tab321freq]\)]\(coefTab[\([All, All, i]\)]\ 
\*SuperscriptBox[\(y1\), \(\(-1.5`\) + \((Tab321freq[\([i, 1]\)] + Tab321freq[\([i, 2]\)] + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]})\)\)] JCx[eps + Tab321freq[\([i, 1]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 2]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, x1vec, y1, eps]\)\),{i1,1,Nmax/2+1},{i2,1,Nmax+1}];
b321p3=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB321I]}];
For[i=1,i<=Length[b321plot],i++,cnc=(Sqrt[y1] knPlot[[i]])^(1.5` +cn[[All,2]]) cn[[All,1]];
For[n=1,n<=Length[x1vec],n++,
For[j=1,j<=Length[biasListB321I],j++,b321plot[[i,n,j]]+=Plin[Sqrt[y1] knPlot[[i]]] Re[(b321p3[[All,All,n,j]].cnc.cnc+(Sqrt[y1] knPlot[[i]])^2 B321UVderivs[[j]]/.{x->1/y1,y->x1vec[[n]]/y1})]]]];

(* Permutation 4: (x,y)\[Rule](1/y,x/y) *)
coefTab2=Tab321derivs231/.{y->y1};
coefTab=Table[coefTab2[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB321I]},{k,1,Length[Tab321freq]}];
Bn1n2n3step1=ParallelTable[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[Tab321freq]\)]\(coefTab[\([All, All, i]\)]\ 
\*SuperscriptBox[\(y1\), \(\(-1.5`\) + \((Tab321freq[\([i, 1]\)] + Tab321freq[\([i, 2]\)] + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]})\)\)]\ JCx[eps + Tab321freq[\([i, 2]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 1]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, x1vec, y1, eps]\)\),{i1,1,Nmax/2+1},{i2,1,Nmax+1}];
b321p4=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB321I]}];
For[i=1,i<=Length[b321plot],i++,cnc=(Sqrt[y1] knPlot[[i]])^(1.5` +cn[[All,2]]) cn[[All,1]];
For[n=1,n<=Length[x1vec],n++,
For[j=1,j<=Length[biasListB321I],j++,b321plot[[i,n,j]]+=Plin[Sqrt[y1] knPlot[[i]]] Re[(b321p4[[All,All,n,j]].cnc.cnc+(Sqrt[y1] knPlot[[i]])^2 B321UVderivs[[j]]/.{x->1/y1,y->x1vec[[n]]/y1})]]]];

(* Permutation 5: (x,y)\[Rule](1/x,y/x) *)
coefTab2=Tab321derivs321/.{y->y1};
coefTab=Table[coefTab2[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB321I]},{k,1,Length[Tab321freq]}];
Bn1n2n3step1=ParallelTable[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[Tab321freq]\)]\(coefTab[\([All, All, i]\)]\ 
\*SuperscriptBox[\(x1vec\), \(\(-1.5`\) + \((Tab321freq[\([i, 1]\)] + Tab321freq[\([i, 2]\)] + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]})\)\)]\ JCx[eps + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 2]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 1]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, x1vec, y1, eps]\)\),{i1,1,Nmax/2+1},{i2,1,Nmax+1}];
b321p5=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB321I]}];
For[i=1,i<=Length[b321plot],i++,
For[n=1,n<=Length[x1vec],n++,cnc=(Sqrt[x1vec[[n]]] knPlot[[i]])^(1.5` +cn[[All,2]]) cn[[All,1]];
For[j=1,j<=Length[biasListB321I],j++,b321plot[[i,n,j]]+=Plin[Sqrt[x1vec[[n]]] knPlot[[i]]] Re[(b321p5[[All,All,n,j]].cnc.cnc+(Sqrt[x1vec[[n]]] knPlot[[i]])^2 B321UVderivs[[j]]/.{x->1/x1vec[[n]],y->y1/x1vec[[n]]})]]]];

(* Permutation 6: (x,y)\[Rule](y/x,1/x) *)
coefTab2=Tab321derivs312/.{y->y1};
coefTab=Table[coefTab2[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB321I]},{k,1,Length[Tab321freq]}];
Bn1n2n3step1=ParallelTable[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[Tab321freq]\)]\(coefTab[\([All, All, i]\)]\ 
\*SuperscriptBox[\(x1vec\), \(\(-1.5`\) + \((Tab321freq[\([i, 1]\)] + Tab321freq[\([i, 2]\)] + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]})\)\)]\ JCx[eps + Tab321freq[\([i, 3]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 1]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, eps + Tab321freq[\([i, 2]\)] /. {\[Nu]1 -> cn[\([i1, 2]\)], \[Nu]2 -> cn[\([i2, 2]\)]}, x1vec, y1, eps]\)\),{i1,1,Nmax/2+1},{i2,1,Nmax+1}];
b321p6=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,i2,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,Nmax+2-i2,n,j]]]],{i1,1,Nmax+1},{i2,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB321I]}];
For[i=1,i<=Length[b321plot],i++,
For[n=1,n<=Length[x1vec],n++,cnc=(Sqrt[x1vec[[n]]] knPlot[[i]])^(1.5` +cn[[All,2]]) cn[[All,1]];
For[j=1,j<=Length[biasListB321I],j++,b321plot[[i,n,j]]+=Plin[Sqrt[x1vec[[n]]] knPlot[[i]]] Re[(b321p6[[All,All,n,j]].cnc.cnc+(Sqrt[x1vec[[n]]] knPlot[[i]])^2 B321UVderivs[[j]]/.{x->y1/x1vec[[n]],y->1/x1vec[[n]]})]]]];

b321plotOut=ConstantArray[0,{Length[b321plot],Length[b321plot[[1]]],Length[biasListB321Iinit]}];
j=1;
For[i=1,i<=Length[biasListB321Iinit],i++,If[biasListB321Icut[[i]]==0,Continue[],Continue[],b321plotOut[[All,All,i]]=b321plot[[All,All,j]];j+=1]];

b321plotOut]


(* ::Code::Initialization::Plain:: *)
If[includeB321Iv1,
Print["Computing B321-I spectra (version 1)"];
time=AbsoluteTiming[
	B321Ispec={};
	For[m=1,m<=Length[xyout],m++,
	b321I=computeB321I[{xyout[[m,1]]},xyout[[m,2]],bias321Iv1];
	B321Ispec=Join[B321Ispec,{b321I[[All,1,All]]}]]
	][[1]];
Print["Runtime: ",time];

SetDirectory[outdir];
Export[ToString[StringForm["``-flat_b321I_nu``_m``.wdx",string,bias321Iv1,muIndex]],B321Ispec];
];

If[includeB321Iv2,
Print["Computing B321-I spectra (version 2)"];
time=AbsoluteTiming[
	B321Ispec={};
	For[m=1,m<=Length[xyout],m++,
	b321I=computeB321I[{xyout[[m,1]]},xyout[[m,2]],bias321Iv2];
	B321Ispec=Join[B321Ispec,{b321I[[All,1,All]]}]]
	][[1]];
Print["Runtime: ",time];

SetDirectory[outdir];
Export[ToString[StringForm["``-flat_b321I_nu``_m``.wdx",string,bias321Iv2,muIndex]],B321Ispec];
];


(* ::Subsubsection::Closed:: *)
(*Compute B321-II spectra*)


(* ::Text:: *)
(*Next we compute the B321-II spectra. These are obtained by first computing P13-like function, then inserting this into the functional form of B321-II (post-symmetrization). This is quite expensive to compute, but depends only linearly on the number of FFTLog points.*)


(* ::Code::Initialization::Plain:: *)
(* Compute FFTLog coefficients *)
cn=CoeffPow[bias321II,Nmax321II,kmin,kmax,Plin];
\[CapitalDelta]=1/(Nmax321II-1) Log[kmax/kmin]//N;
kn=Table[kmin Exp[(i-1) \[CapitalDelta]],{i,1,Nmax321II}];

(* Compute P13 *)
P13biases=b3sII;
p13s={};
For[n=1,n<=Length[P13biases],n++,
p13matrix=ParallelTable[Coefficient[P13ab[-0.5cn[[j,2]]],P13biases[[n]]]/.{b1->0,g2->0,g21->0,f->0},{j,1,Nmax321II+1}];
p13=ParallelTable[{kn[[i]],kn[[i]]^3 PPlin[kn[[i]]]ComplexExpand[Re[p13matrix.(kn[[i]]^cn[[All,2]] cn[[All,1]])]]+kn[[i]]^2 PPlin[kn[[i]]]Coefficient[\[Alpha]UV13,P13biases[[n]]]/.{b1->0,b2->0,g2->0,g21->0,f->0}},{i,1,Nmax321II}];
p13s=Join[p13s,{p13}]];
CoeffSafe[x_,fac_,pow_]:=If[pow==0,x/. {fac->0},Coefficient[x,fac^pow]];
P13loc[kT_,\[Mu]T_]:=Sum[P13biases[[i]]Sum[\[Mu]T^n Interpolation[Transpose[Join[{p13s[[i]][[All,1]],CoeffSafe[p13s[[i]][[All,2]],\[Mu]P,n]}]],"ExtrapolationHandler" -> {(0.0 &), "WarningMessage" -> False}][kT],{n,0,12}],{i,1,Length[P13biases]}];
b321IIxy=getCoeff[((B321II[k,Sqrt[y]k,Sqrt[x]k]/.{P13->P13loc}))/.{\[Mu]1->\[Mu],\[Mu]2->-\[Mu],\[Mu]3->-\[Mu]}]/.{PPlin->Plin};

(* Function to compute B321-II *)
computeB321II[x1vec_,y1_]:=Module[{Bn1n2n3step1,P13loc,p13,p13s,P13biases,b321IIplot,cnc,Nmax,bias,cn},
(* Compute output [NB: this is already symmetrized] *)
b321IIplot=ParallelTable[Coefficient[b321IIxy/.{k->knPlot[[i]],y->Sqrt[y1],x->Sqrt[x1vec[[j]]]},biasListB321II[[n]]]/.{b1->0,b2->0,g2->0,g21->0,f->0},{i,1,Length[knPlot]},{j,1,Length[x1vec]},{n,1,Length[biasListB321II]}];
Transpose[b321IIplot]];


(* ::Code::Initialization::Plain:: *)
If[includeB321II,
Print["Computing B321-II spectra"];
time=AbsoluteTiming[
	B321IIspec={};
	For[m=1,m<=Length[xyout],m++,
	b321II=computeB321II[{xyout[[m,1]]},xyout[[m,2]]];
	B321IIspec=Join[B321IIspec,{b321II[[1,All,All]]}]]
	][[1]];
Print["Runtime: ",time];
];

SetDirectory[outdir];
If[includeB321II,Export[ToString[StringForm["``-flat_b321II_nu``_m``.wdx",string,bias321II,muIndex]],B321IIspec]];


(* ::Subsubsection::Closed:: *)
(*Compute B411 spectra*)


(* ::Text:: *)
(*Finally, we compute B411. This is analogous to the other loop-integrals, but linear in the number of FFTLog points, and includes a sum over three permutations.*)


computeB411[x1vec_,y1_]:=
Module[{Bn1n2n3step1,tab411,b411,b411p1,b411p2,b411p3,b411plot,cnc,Nmax,bias,cn,biasZeros, coefTab},
(* Compute FFTlog coefficients *)
Nmax=Nmax411;
bias=bias411;
cn=CoeffPow[bias,Nmax,kmin,kmax,Plin];

biasZeros={g21x->0,g211->0,g22->0,g31->0,g21->0,g3->0,b3->0,g2x->0,g2->0,b2->0,b1->0,f->0};

(* Define output array *)
b411plot=ConstantArray[0,{NmaxPlot,Length[x1vec],Length[biasListB411]}];

(* Permutation 1: (x,y)\[Rule](x,y) *)
tab411=Tab411derivs123/.{y->y1};
coefTab=Table[tab411[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB411]},{k,1,Length[Tab411freq]}];
Bn1n2n3step1=ParallelTable[Sum[coefTab[[All,All,i]]JCx[eps+Tab411freq[[i,1]]/.{\[Nu]->cn[[i1,2]]},eps+Tab411freq[[i,2]]/.{\[Nu]->cn[[i1,2]]},eps+Tab411freq[[i,3]]/.{\[Nu]->cn[[i1,2]]},x1vec,y1,eps],{i,1,Length[Tab411freq]}],{i1,1,Nmax/2+1}];
b411p1=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,n,j]]]],{i1,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB411]}];
For[i=1,i<=Length[b411plot],i++,
For[n=1,n<=Length[x1vec],n++,
cnc=knPlot[[i]]^(3+cn[[All,2]]) cn[[All,1]];
For[j=1,j<=Length[biasListB411],j++,
b411plot[[i,n,j]]+=Plin[knPlot[[i]]]Plin[Sqrt[y1]knPlot[[i]]](Re[b411p1[[All,n,j]].cnc]+Coefficient[knPlot[[i]]^2B411UV123[x1vec[[n]],y1],biasListB411[[j]]]/.biasZeros)]]];

(* Permutation 2: (x,y)\[Rule](1/y,x/y) *)
(* Note that permutations have already been counted! *)
tab411=Tab411derivs231/.{y->y1};
coefTab=Table[tab411[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB411]},{k,1,Length[Tab411freq]}];
Bn1n2n3step1=Table[Sum[coefTab[[All,All,i]]y1^(-1.5+(Tab411freq[[i,1]]+Tab411freq[[i,2]]+Tab411freq[[i,3]]/.{\[Nu]->cn[[i1,2]]})) JCx[eps+Tab411freq[[i,2]]/.{\[Nu]->cn[[i1,2]]},eps+Tab411freq[[i,3]]/.{\[Nu]->cn[[i1,2]]},eps+Tab411freq[[i,1]]/.{\[Nu]->cn[[i1,2]]},x1vec,y1,eps],{i,1,Length[Tab411freq]}],{i1,1,Nmax/2+1}];
b411p2=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,n,j]]]],{i1,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB411]}];
For[i=1,i<=Length[b411plot],i++,
For[n=1,n<=Length[x1vec],n++,
cnc=(Sqrt[y1]knPlot[[i]])^(3+cn[[All,2]]) cn[[All,1]];
For[j=1,j<=Length[biasListB411],j++,
b411plot[[i,n,j]]+=Plin[Sqrt[x1vec[[n]]]knPlot[[i]]]Plin[Sqrt[y1]knPlot[[i]]](Re[b411p2[[All,n,j]].cnc]+Coefficient[(Sqrt[y1]knPlot[[i]])^2 B411UV231[x1vec[[n]],y1],biasListB411[[j]]]/.biasZeros)]]];

(* Permutation 3: (x,y)\[Rule](y,x) *)
tab411=Tab411derivs132/.{y->y1};
coefTab=Table[tab411[[j,k]]/.{x->x1vec[[i]]},{i,1,Length[x1vec]},{j,1,Length[biasListB411]},{k,1,Length[Tab411freq]}];
Bn1n2n3step1=Table[Sum[coefTab[[All,All,i]]JCx[eps+Tab411freq[[i,2]]/.{\[Nu]->cn[[i1,2]]},eps+Tab411freq[[i,1]]/.{\[Nu]->cn[[i1,2]]},eps+Tab411freq[[i,3]]/.{\[Nu]->cn[[i1,2]]},x1vec,y1,eps],{i,1,Length[Tab411freq]}],{i1,1,Nmax/2+1}];
b411p3=Table[If[i1<=Nmax/2+1,Bn1n2n3step1[[i1,n,j]],Conjugate[Bn1n2n3step1[[Nmax+2-i1,n,j]]]],{i1,1,Nmax+1},{n,1,Length[x1vec]},{j,1,Length[biasListB411]}];
For[i=1,i<=Length[b411plot],i++,
For[n=1,n<=Length[x1vec],n++,
cnc=knPlot[[i]]^(3+cn[[All,2]]) cn[[All,1]];
For[j=1,j<=Length[biasListB411],j++,
b411plot[[i,n,j]]+=Plin[knPlot[[i]]]Plin[Sqrt[x1vec[[n]]]knPlot[[i]]](Re[b411p3[[All,n,j]].cnc]+Coefficient[knPlot[[i]]^2 B411UV132[x1vec[[n]],y1],biasListB411[[j]]]/.biasZeros)]]];

b411plotOut=ConstantArray[0,{Length[b411plot],Length[b411plot[[1]]],Length[biasListB411init]}];
j=1;
For[i=1,i<=Length[biasListB411init],i++,If[biasListB411cut[[i]]==0,Continue[],Continue[],b411plotOut[[All,All,i]]=b411plot[[All,All,j]];j+=1]];

b411plotOut]


(* ::Code::Initialization::Plain:: *)
If[includeB411,
Print["Computing B411 spectra"];
time=AbsoluteTiming[
	B411spec={};
	For[m=1,m<=Length[xyout],m++,
	b411=computeB411[{xyout[[m,1]]},xyout[[m,2]]];
	B411spec=Join[B411spec,{b411[[All,1,All]]}]]
	][[1]];
Print["Runtime: ",time];
];

SetDirectory[outdir];
If[includeB411,Export[ToString[StringForm["``_b411-flat_nu``_m``.wdx",string,bias411,muIndex]],B411spec]];



