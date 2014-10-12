(* ::Package:: *)

(* Utility functions *)

SubtractList[eall_,es_]:=Delete[eall,First[Position[eall,#]]&/@es]


(* Calculate irreducible LS eigenspace quantum numbers *)

quantML[quantL_]:=Flatten[{#,#}&/@Range[quantL,-quantL,-1]]
quantMS[quantL_]:=Flatten[KroneckerProduct[ConstantArray[1,2quantL+1],{1/2,-1/2}]]

(* Calculate all Lz-Sz quantum numbers *)
CalcQuantLzSz[quantL_,n_]:=Module[{
		fm=FermiMap[2(2quantL+1),n],
		qmL=quantML[quantL],
		qmS=quantMS[quantL]},
	{Total[qmL[[#]]],Total[qmS[[#]]]}&/@(FermiToCoords/@fm)]

QuantMRange[l_,s_]:=Flatten[Outer[{##}&,Range[l,-l,-1],Range[s,-s,-1]],1]

(* Irreducible LS eigenspace quantum numbers *)
CalcIrredLSQuant[quantLzSz_]:=Module[{irred={},qm},
	(* lexicograpical ordering *)
	qm=Sort[quantLzSz,First[#1]>First[#2]\[Or](First[#1]==First[#2]\[And]Last[#1]>Last[#2])&];
	While[qm!={},
		irred=Append[irred,First[qm]];
		qm=SubtractList[qm,QuantMRange@@First[qm]]];
	irred]


(* Angular momentum and spin (explicit matrix representation) *)

AngularMomentumUp[quantL_]:=If[quantL==0,SparseArray[{},{2,2}],
	KroneckerProduct[SparseArray[DiagonalMatrix[Sqrt[quantL (quantL+1)-# (#+1)]&/@Range[quantL-1,-quantL,-1],1]],IdentityMatrix[2]]]

AngularMomentumDown[quantL_]:=If[quantL==0,SparseArray[{},{2,2}],
	KroneckerProduct[SparseArray[DiagonalMatrix[Sqrt[quantL (quantL+1)-# (#-1)]&/@Range[quantL,-quantL+1,-1],-1]],IdentityMatrix[2]]]

AngularMomentumZ[quantL_]:=
	KroneckerProduct[SparseArray[DiagonalMatrix[Range[quantL,-quantL,-1]]],IdentityMatrix[2]]

AngularMomentum[quantL_]:=Module[{Lu,Ld,Lz},
		Lu=AngularMomentumUp[quantL];
		Ld=AngularMomentumDown[quantL];
		Lz=AngularMomentumZ[quantL];
		{(Lu+Ld)/2,-I(Lu-Ld)/2,Lz}]

SpinUp[quantL_]  :=KroneckerProduct[SparseArray[IdentityMatrix[2quantL+1]],{{0,1},{0,0}}]
SpinDown[quantL_]:=KroneckerProduct[SparseArray[IdentityMatrix[2quantL+1]],{{0,0},{1,0}}]
SpinZ[quantL_]   :=KroneckerProduct[SparseArray[IdentityMatrix[2quantL+1]],{{1/2,0},{0,-1/2}}]

Spin[quantL_]:=Module[{Su,Sd,Sz},
		Su=SpinUp[quantL];
		Sd=SpinDown[quantL];
		Sz=SpinZ[quantL];
		{(Su+Sd)/2,-I(Su-Sd)/2,Sz}]

(* angular momentum and spin symmetry quantum numbers *)
Symmetry[V_]:=(Dimensions[V,2]-1)/2

SymmetryToString[l_,s_,parity_]:=Module[{Lnames={"S","P","D","F","G","H","I","J","K","L","M","N","O"}},
	"^"<>ToString[2s+1]<>Lnames[[l+1]]<>If[parity==-1,"^o",""]]


(* Calculate irreducible LS eigenspaces *)

(* Partition the fermionic base indices into Lz-Sz eigenspaces *)
CalcLzSzBaseIX[quantLzSz_]:=Module[{lsmax,baseIX,i},
	lsmax={Max[First/@quantLzSz],Max[Last/@quantLzSz]};
	baseIX=ConstantArray[{},2lsmax+1];
	MapIndexed[(i=lsmax-#1+1;
		baseIX[[First[i],Last[i]]]=Append[baseIX[[First[i],Last[i]]],First[#2]])&,quantLzSz];
	baseIX]

(* Extract a normalized vector v from the range of a projection matrix P.
	In case v is not unique (rank (P)>1), try to find a v with simple coefficients,
	i.e., the numerators and denominators should be as small as possible. *)
ExtractProjVector[P_]:=Module[{c,Z,v1,v2},
	If[Length[P]<16,
		(* try 'NullSpace' explicitly, only works for small dimensions *)
		c=Sqrt[LCM@@Denominator[Flatten[P]^2]];
		Z=Normalize/@FullSimplify[NullSpace[c (P-IdentityMatrix[Length[P]])]];
		(* select vector with smallest maximal prime factor in denominator *)
		Z[[First[Ordering[Max[First/@FactorInteger[LCM@@Denominator[#^2]]]&/@Z,1]]]],
		(* else resort to QR decomposition, which works typically faster than NullSpace! *)
		Z=First[QRDecomposition[P]];
		v1=First[Z]; v2=FullSimplify[Normalize[First[RowReduce[Z]]]];
		If[LCM@@Denominator[v1^2]<LCM@@Denominator[v2^2],v1,v2]]]

SubtractProjVector[P_,v_]:=FullSimplify[P-KroneckerProduct[v,v]]

SubtractIrredLS[projLzSz_,V_]:=Module[{i},
	Map[(
		i=Symmetry[V]-First[#]+1;
		If[First[i]>0\[And]Last[i]>0,{First[#],SubtractProjVector[Last[#],V[[First[i],Last[i]]]]},#])&,
	projLzSz]]

CreateLadderMap[Ld_,baseIX_]:=Map[Ld[[Last[#],First[#]]]&,Transpose[{baseIX,RotateLeft[baseIX]},{3,1,2}],{2}]

ApplyLadderMap[l_,Ld_,\[Psi]_]:=Module[{\[Phi]=\[Psi]},Join[{\[Psi]},Outer[Function[m,\[Phi]=FullSimplify[Ld[[l-m+1]].\[Phi]/Sqrt[l (l+1)-m (m-1)]]],Range[l,-l+1,-1]]]]

CalcIrredLS[quantL_,n_]:=Module[{quantLzSz,irredSym,baseIX,projLzSz,Ld,Sd,LdMap,SdMap,\[Psi],V,i},
	quantLzSz=CalcQuantLzSz[quantL,n];
	irredSym=CalcIrredLSQuant[quantLzSz];
	projLzSz={#,IdentityMatrix[Count[quantLzSz,#]]}&/@DeleteDuplicates[irredSym];
	Ld=p2N[AngularMomentumDown[quantL],2(2quantL+1),1,n];
	Sd=p2N[SpinDown[quantL],2(2quantL+1),1,n];
	baseIX=CalcLzSzBaseIX[quantLzSz];
	LdMap=CreateLadderMap[Ld,baseIX];
	SdMap=CreateLadderMap[Sd,Transpose[baseIX]];
	Function[ls,
		i=First[Flatten[Position[projLzSz,{ls,_},{1}]]];
		\[Psi]=ExtractProjVector[projLzSz[[i,2]]];
		i=Symmetry[baseIX]-ls+1;
		V=ApplyLadderMap[First[ls],LdMap[[First[i];;-First[i],Last[i]]],\[Psi]];
		V=MapIndexed[ApplyLadderMap[Last[ls],SdMap[[Last[i];;-Last[i],First[i]+First[#2]-1]],#]&,V];
		projLzSz=SubtractIrredLS[projLzSz,V];
		SparseArray[Map[SparseArray[First[#]->Last[#],{Length[Ld]}]&,
			Transpose[{baseIX[[First[i];;-First[i],Last[i];;-Last[i]]],V},{3,1,2}],{2}]]]
	/@irredSym]


(* Check computed irreducible LS eigenspaces *)

CheckIrredLS[quantL_,n_,irredLS_]:=Module[{L,S,LL,SS,l,s},
	L=p2N[#,2(2quantL+1),1,n]&/@AngularMomentum[quantL];
	S=p2N[#,2(2quantL+1),1,n]&/@Spin[quantL];
	LL=Total[#.#&/@L];
	SS=Total[#.#&/@S];
	Total[Function[V,
		{l,s}=Symmetry[V];
		Norm[MapIndexed[
			Norm[PurgeSparseArray[L[[3]].#-(l-First[#2]+1)#]]+
			Norm[PurgeSparseArray[S[[3]].#-(s-Last[#2]+1)#]]+
			Norm[PurgeSparseArray[LL.#-l (l+1)#]]+
			Norm[PurgeSparseArray[SS.#-s (s+1)#]]+
			FullSimplify[Norm[#]-1]&,
		V,{2}]]
	]/@irredLS]]
