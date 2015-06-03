(* ::Package:: *)

SlaterToTeX[orbnames_,s_]:="\\ket{"<>StringJoin@@Flatten[{#," "}&/@orbnames[[FermiToCoords[s]]]][[1;;-2]]<>"}"

FermiToTeX[orbnames_,fm_,\[Psi]_]:=Module[{f=LCM@@Denominator[\[Psi]^2],ret},
	ret=StringJoin@@MapIndexed[If[First[#2]==1,#1,If[First[Characters[#1]]!="-","+",""]<>#1]&,
		If[First[#]==1,"",If[First[#]==-1,"-",ToString[TeXForm[First[#]]]<>"\\cdot"]]<>SlaterToTeX[orbnames,#[[2]]]&/@
			Select[Transpose[{Sqrt[f]\[Psi],fm}],First[#]!=0&]];
	If[f!=1,ToString[TeXForm[1/Sqrt[f]]]<>"\\left("<>ret<>"\\right)",ret]]


FermiLineBreak[\[Psi]tex_,width_]:=Module[{breakpos,ketstart=First[First[StringPosition[\[Psi]tex,"ket"]]]},
	breakpos=Select[First[#]&/@StringPosition[\[Psi]tex,{"+","-"}],#>ketstart&];
	If[Length[breakpos]<width,\[Psi]tex,
		StringInsert[\[Psi]tex,"\\right.$ $\\vphantom{\\biggl(}$\\\\\n&&&&$\\left.",breakpos[[width;; ;;width]]]]]


(* Irreducible LS eigenspaces *)

SymmetryToTex[l_,s_,parity_]:=Module[{Lnames={"\\mathrm{S}","\\mathrm{P}","\\mathrm{D}","\\mathrm{F}","\\mathrm{G}","\\mathrm{H}","\\mathrm{I}","\\mathrm{J}","\\mathrm{K}","\\mathrm{L}","\\mathrm{M}","\\mathrm{N}","\\mathrm{O}"}},
	"^"<>ToString[2s+1]<>Lnames[[l+1]]<>If[parity==-1,"^{\\mathrm{o}}",""]]

IrredWriteHeader[fid_]:=WriteString[fid,"\n\\begin{table}[!ht]\n\\centering\n"<>"\\begin{tabular}{|c|c|cc|c|}\n"<>
	"\\hline\nconfig&sym&$L_z$&$S_z$&$\\Psi$ $\\vphantom{\\Bigl(}$\\\\\n"<>"\\hline\\hline\n"];

IrredWriteClosure[fid_,count_]:=(
	WriteString[fid,"\\end{tabular}\n"];
	If[count==1,WriteString[fid,"\\caption{Irreducible LS eigenspaces, showing states with maximal $L_z, S_z$ only}\n"],
		WriteString[fid,"\\caption{Irreducible LS eigenspaces (continued)}\n"]];
	WriteString[fid,"\\label{tab:irredLS"<>ToString[count]<>"}\n\\end{table}\n\n"])

IrredLSToTex[quantL_,n_,LS_,fid_,linestart_,tabstart_,width_]:=Module[{fm=FermiMap[2(2quantL+1),n],
	orbnames,linecount=linestart,tabcount=tabstart},
	orbnames=Switch[quantL,
		0,{"\\mathrm{s}","\\conj{\\mathrm{s}}"},
		1,{"\\mathrm{p}_1","\\conj{\\mathrm{p}_1}","\\mathrm{p}_0","\\conj{\\mathrm{p}_0}","\\mathrm{p}_{\\neg1}","\\conj{\\mathrm{p}_{\\neg1}}"},
		2,{"\\mathrm{d}_2","\\conj{\\mathrm{d}_2}","\\mathrm{d}_1","\\conj{\\mathrm{d}_1}","\\mathrm{d}_0","\\conj{\\mathrm{d}_0}","\\mathrm{d}_{\\neg1}","\\conj{\\mathrm{d}_{\\neg1}}","\\mathrm{d}_{\\neg2}","\\conj{\\mathrm{d}_{\\neg2}}"},
		3,{"\\mathrm{f}_3","\\conj{\\mathrm{f}_3}","\\mathrm{f}_2","\\conj{\\mathrm{f}_2}","\\mathrm{f}_1","\\conj{\\mathrm{f}_1}","\\mathrm{f}_0","\\conj{\\mathrm{f}_0}","\\mathrm{f}_{\\neg1}","\\conj{\\mathrm{f}_{\\neg1}}","\\mathrm{f}_{\\neg2}","\\conj{\\mathrm{f}_{\\neg2}}","\\mathrm{f}_{\\neg3}","\\conj{\\mathrm{f}_{\\neg3}}"}];
	MapIndexed[
		(If[linecount>25,
			If[First[#2]>1,WriteString[fid,"\\hline\n"]];
				IrredWriteClosure[fid,tabcount];IrredWriteHeader[fid];
				linecount=0;tabcount++,
			If[First[#2]>1,WriteString[fid,"\\cline{2-5}\n"]]];
		If[First[#2]==1,WriteString[fid,"$\\wedge^{"<>ToString[n]<>"}V_{\\mathrm{"<>{"s","p","d","f","g"}[[quantL+1]]<>"}}$"]];
		WriteString[fid,"&$"<>SymmetryToTex@@Append[Symmetry[#1],(-1)^(n quantL)]<>"$&"];
		WriteString[fid,"$"<>ToString[TeXForm[#]]<>"$&"]&/@Symmetry[#1];
		WriteString[fid,"$"<>StringReplace[FermiLineBreak[FermiToTeX[orbnames,fm,#1[[1,1]]],width],"\\neg"->"\\text{-}"]<>"$ $\\vphantom{\\Bigl(}$\\\\\n"];
		linecount+=Ceiling[Length[Select[#1[[1,1]],#!=0&]]/width])&,
		LS];
	WriteString[fid,"\\hline\n"];
	{linecount,tabcount}]

IrredToTex[irredLS_,filename_,widthtable_]:=Module[{fid,linecount=0,tabcount=1,Lnames={"s","p","d","f","g"}},
	fid=OpenWrite[filename];
	IrredWriteHeader[fid];
	Function[quantL,
		Function[n,Print[Lnames[[quantL+1]]^n];
			{linecount,tabcount}=IrredLSToTex[quantL,n,irredLS[[quantL+1,n]],fid,linecount,tabcount,widthtable[[quantL+1,n]]]]/@
		Range[Length[irredLS[[quantL+1]]]];
		If[quantL<Length[irredLS]-1,WriteString[fid,"\\hline\n"]]]/@
	Range[0,Length[irredLS]-1];
	IrredWriteClosure[fid,tabcount];
	Close[fid];]
