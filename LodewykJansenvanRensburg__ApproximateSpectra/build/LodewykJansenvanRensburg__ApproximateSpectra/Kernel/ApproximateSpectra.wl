(* ::Package:: *)

BeginPackage["LodewykJansenvanRensburg`ApproximateSpectra`"]

(* Declare your package's public symbols here. *)


(* Declare your package's public symbols here. *)
(* Declare your package's public symbols here. *)
charPoly::usage = "charPoly is the characteristic polynomial computed by AnalyzeCharPoly.";
charPolyMinus::usage = "charPolyMinus is the characteristic polynomial of the minus case in AnalyzeCharPoly.";
z::usage = "Variable of charateristic polynomial.";
AnalyzeCharPoly::usage = "AnalyzeCharPoly[A] computes the characteristic polynomial and other spectral data.";

Begin["`Private`"]

(* Define your public and private symbols here. *)

(* 1. Generate Random Haar Unitaries *)
(* Full function to generate random complex Gaussian matrix and apply QR decomposition *)
LodewykJansenvanRensburg`ApproximateSpectra`GenerateHaarUnitaryMatrix[n_, \[Sigma]_: 1] := Module[{A, Q, R},
  (* Generate random complex Gaussian matrix *)
  A = RandomVariate[NormalDistribution[0, \[Sigma]], {n, n}] + 
      I RandomVariate[NormalDistribution[0, \[Sigma]], {n, n}];

  (* QR decomposition to get orthonormal matrix Q *)
  {Q, R} = QRDecomposition[A];

  (* Return the unitary matrix Q *)
  Q
];
(* 2. Generate Self-adjoint Tuples --------------------------------------------- *)
(* n = size, g = tuple length *)
LodewykJansenvanRensburg`ApproximateSpectra`RandomSelfAdjointMatrixTuple[n_, g_] := 
  Table[
    Module[{a},
      a = RandomComplex[{-1 - I, 1 + I}, {n, n}];
      (a + ConjugateTranspose[a])/2
    ],
    {g}
  ];
 
  
   
(* 3. Generate Self-adjoint Tuples --------------------------------------------- *)
Options[LodewykJansenvanRensburg`ApproximateSpectra`EstimateSpectralRadius] = {
  Verbose -> False,
  Extremum -> Automatic (* Options: "Max", "Min", or Automatic *)
};
(*


*)
LodewykJansenvanRensburg`ApproximateSpectra`EstimateSpectralRadius[matrixTuple_, mode_: Automatic, opts : OptionsPattern[]] := Module[
  {Aconj, Aneg, lambda, outerSpec, innerSpec, verboseQ, extremum, computeValue},

  verboseQ = OptionValue[Verbose];
  extremum = OptionValue[Extremum];

  (* Lambda operator: sum of conj(Ai) \[CircleTimes] Ai *)
  lambda[t1_, t2_] := Total[
    MapThread[KroneckerProduct[Conjugate[#1], #2] &, {t1, t2}]
  ];

  Aconj = Conjugate /@ matrixTuple;

	 (* Handle inner spec negation safely *)
  Aneg = If[
    Length[Aconj] >= 2,
     KroneckerProduct[Conjugate[matrixTuple[[1]]],matrixTuple[[1]]]-KroneckerProduct[Conjugate[matrixTuple[[2]]],matrixTuple[[2]]]
  ];

  (* Function to compute either max, min, or spectral radius *)
  computeValue[M_] := Module[{vals},
    vals = N[Abs[Eigenvalues[M]]];
   
    Switch[extremum,
      "Max", Max[Abs[vals]],
      "Min", Min[Abs[vals]],
      _, Max[Abs[vals]] (* Spectral radius *)
    ]
  ];
  (* *)
  Aconj =  KroneckerProduct[Conjugate[matrixTuple[[1]]],matrixTuple[[1]]]+KroneckerProduct[Conjugate[matrixTuple[[2]]],matrixTuple[[2]]];
	(* *)
	
  outerSpec = Sqrt[computeValue[Aconj]];
  innerSpec = Sqrt[computeValue[Aneg] ];

  If[verboseQ,
    Print["Plus: ",MatrixForm[Aconj]];
    Print["Minus: ",MatrixForm[Aneg]];
    Print["Outer Estimate (", extremum, "): ", N[outerSpec]];
    Print["Inner Estimate (", extremum, "): ", N[innerSpec]];
  ];

  Switch[mode,
    "outer", outerSpec,
    "inner", innerSpec,
    Automatic | "both", {innerSpec, outerSpec},
    _, Message[LodewykJansenvanRensburg`ApproximateSpectra`EstimateSpectralRadius::badmode]; $Failed
  ]
];

LodewykJansenvanRensburg`ApproximateSpectra`EstimateSpectralRadius::badmode = 
  "Mode must be one of \"inner\", \"outer\", or omitted to return both.";

(* Returns minimum Eigenvalue*)
LodewykJansenvanRensburg`ApproximateSpectra`SpectralMinimum[A_?MatrixQ] := Module[
  {eigvals},
  eigvals = Eigenvalues[A];
  N[Min[Abs[eigvals]]];
];
(* _____________________________________________________________________________ *)
(* This computes the two norm of A1u1 +...Agug*)
LodewykJansenvanRensburg`ApproximateSpectra`Norm2Vneum[mats_List] := Module[
  {conjProdSum, n, traceVal, radius},
  
  (* Sum of A_i^* A_i *)
  conjProdSum = Total[ConjugateTranspose[#] . # & /@ mats];
  
  (* Matrix size *)
  n = Length[conjProdSum];
  
  (* Normalized trace *)
  traceVal = Tr[conjProdSum]/n;
  
  (* Final result *)
  radius = Re[Sqrt[traceVal]];
  
  radius
];

(* 5. ------------------------- Function that does testing -------------------------- *)
Options[LodewykJansenvanRensburg`ApproximateSpectra`PlotSpectrum] = {PlotStyle -> "2D", NumBins -> 20};
LodewykJansenvanRensburg`ApproximateSpectra`PlotSpectrum[AU_, opts : OptionsPattern[]] := Module[
  {
     eigval1, eigMax,  hist,
    plotStyle, numBins, legendLabels, basePlot
  },
   plotStyle = OptionValue[PlotStyle];
  numBins = OptionValue[NumBins];
  eigval1 = {Re[#], Im[#]} & /@ Eigenvalues[AU];
  eigMax = N[Max[Abs[eigval1]]];
  (* Histogram or scatter depending on style *)
  basePlot = Switch[plotStyle,
    "3D",
    Histogram3D[
      eigval1, numBins, "Probability",
      PlotRange -> {{-eigMax - 1, eigMax + 1}, {-eigMax - 1, eigMax + 1}, All},
      ColorFunction -> (Blend[{White, Blue}, #] &),
      ChartLegends -> Automatic
    ],
    
    "2D",
    ListPlot[
      eigval1,
      PlotStyle -> Blue,
      AxesLabel -> {"Re", "Im"},
      AspectRatio -> 1,
      PlotRange -> {{-eigMax - 1, eigMax + 1}, {-eigMax - 1, eigMax + 1}},
      Frame -> True,
      FrameLabel -> {"Re", "Im"},
      GridLines -> Automatic
    ],
    
    _, (Message[LodewykJansenvanRensburg`ApproximateSpectra`PlotSpectrum::badstyle]; Return[$Failed])
  ];

  (* Points to plot depending on 2D vs 3D *)
  hist = Switch[plotStyle,
    "3D",
    Show[
      basePlot,
      Graphics3D[{Red, PointSize[Large], Point[{0,eigMax , 0}]}],
      ImageSize -> 500
    ],

    "2D",
    Show[
      basePlot,
      Graphics[{
        Red, PointSize[Large], Point[{0, eigMax}]
      }],
      ImageSize -> 500
    ]
  ];
(* Legend as side column *)
  legendLabels = Column[{
    Style["Legend", Bold, 14],
    Style[eigMax , ": Spectral Max. " , Red],
    Style[SafeDivide[1,eigMax] , ": 1/Spectral Max. (inner radius of original)" ]
    
  }, Spacings -> 1.5];

  Grid[{{hist, legendLabels}}, Spacings -> {2, 1}]
];
LodewykJansenvanRensburg`ApproximateSpectra`PlotSpectrum::badstyle = 
  "PlotStyle must be either \"2D\" or \"3D\".";
  
LodewykJansenvanRensburg`ApproximateSpectra`SafeDivide[a_, b_] := If[b == 0, "Undefined", a/b];

(* 5. ------------------------- Function that does testing -------------------------- *)
ClearAll[LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrum]
Options[LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrum] = {PlotStyle -> "3D", NumBins -> 20};

LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrum[A1_, U_, opts : OptionsPattern[]] := Module[
  {
    AU, eigval1, Ac1, lambda, hist, 
    SumMax, SumMin, DiffMax, DiffMin,  norm1,
    plotStyle, numBins, legendLabels, basePlot,sqrtEigValAUPlus,sqrtEigValAUMinus,
    AnnuliRadiiPlus,AnnuliRadiiMinus,AnnuliRadiiJoin,histInverse,nonZeroEig,eigvalRecip
  },

  plotStyle = OptionValue[PlotStyle];
  numBins = OptionValue[NumBins];

  (* Lambda operator *)
  lambda[t1_, t2_] := Total[
    MapThread[KroneckerProduct[Conjugate[#1], #2] &, {t1, t2}]
  ];

  (* Matrix of interest *)
  Ac1 = Conjugate /@ A1;
  AU = lambda[A1, U];

  (* Multiplicity *)	
	sqrtEigValAUPlus = Sqrt[
           N[
               Abs[ 
                   Eigenvalues[ 
                                 KroneckerProduct[ Conjugate[ A1[[1]] ],A1[[1]]] +KroneckerProduct[Conjugate[A1[[2]]],A1[[2]]] 
                               ] 
                  ]
            ] 
        ];
        	sqrtEigValAUMinus = Sqrt[
           N[
               Abs[ 
                   Eigenvalues[ 
                                 KroneckerProduct[Conjugate[A1[[1]]],A1[[1]]]-KroneckerProduct[Conjugate[A1[[2]]],A1[[2]]]
                               ] 
                  ]
            ] 
        ];
  
   (* 10^-8 accuracy*)
   (* Plus *)
    sqrtEigValAUPlus = N[Round[#, 10^-8] & /@ sqrtEigValAUPlus];
	AnnuliRadiiPlus = Select[sqrtEigValAUPlus, Count[sqrtEigValAUPlus, #] == 1 &]; (* Only Keeps multiplicity one eig values*)
   (* Minus *)
    sqrtEigValAUMinus = N[Round[#, 10^-8] & /@ sqrtEigValAUMinus];
	AnnuliRadiiMinus = Select[sqrtEigValAUMinus, Count[sqrtEigValAUMinus, #] == 1 &]; (* Not needed now -  Only Keeps multiplicity one eig values*)

	AnnuliRadiiJoin = Union[sqrtEigValAUPlus,sqrtEigValAUMinus]; (* ---------------------------- Change Join to Union ------------ *)
  (* Spectral radius estimates *)
  SumMax = N[Abs[EstimateSpectralRadius[A1, "outer",  Verbose -> False, Extremum -> "Max"]]];
  DiffMax = Abs[EstimateSpectralRadius[A1, "inner", Verbose -> False, Extremum -> "Max"]];
  SumMin = Abs[EstimateSpectralRadius[A1, "outer", Verbose -> False, Extremum -> "Min"]];
  DiffMin = Abs[EstimateSpectralRadius[A1, "inner", Verbose -> False, Extremum -> "Min"]];
  norm1 = Norm2Vneum[A1];

  (* Eigenvalues and xy projection *)
  eigval1 = Eigenvalues[AU];
  nonZeroEig = Select[eigval1, Abs[#] != 0 &];
  eigval1 = {Re[#], Im[#]} & /@ eigval1;
  eigvalRecip = 1/nonZeroEig;
  eigvalRecip = {Re[#], Im[#]} & /@ eigvalRecip;
  
  
  
    Print["==================================================================================================="];
	Print[" AnalyzeTupleSpectrum: "];
	Print["===================================================================================================\n"];
  Print["Sqrt of eig values of  A1bar \[CircleTimes] A1 + A2bar \[CircleTimes] A2: ",
    Sqrt[N[Abs[Eigenvalues[KroneckerProduct[Conjugate[A1[[1]]],A1[[1]]]+KroneckerProduct[Conjugate[A1[[2]]],A1[[2]]]]]]]
    ];
     Print["Sqrt of eig values of A1bar \[CircleTimes] A1 - A2bar \[CircleTimes] A2 : ",
    Sqrt[N[Abs[Eigenvalues[KroneckerProduct[Conjugate[A1[[1]]],A1[[1]]]-KroneckerProduct[Conjugate[A1[[2]]],A1[[2]]]]]]]
    , "\n"];
	Print["Sqrt of eig values of  A1bar \[CircleTimes] A1 + A2bar \[CircleTimes] A2, and A1bar \[CircleTimes] A1 - A2bar \[CircleTimes] A2: \n", AnnuliRadiiJoin ];
    Print["==================================================================================================="];

  (* Histogram or scatter depending on style *)
  basePlot = Switch[plotStyle,
    "3D",
    Histogram3D[
      eigval1, numBins, "Probability",
      PlotRange -> {{-SumMax - 1, SumMax + 1}, {-SumMax - 1, SumMax + 1}, All},
      ColorFunction -> (Blend[{White, Blue}, #] &),
      ChartLegends -> Automatic
    ],
    
    "2D",
    ListPlot[
      eigval1,
      PlotStyle -> Blue,
      AxesLabel -> {"Re", "Im"},
      AspectRatio -> 1,
      PlotRange -> {{-SumMax - 1, SumMax + 1}, {-SumMax - 1, SumMax + 1}},
      Frame -> True,
      FrameLabel -> {"Re", "Im"},
      GridLines -> Automatic
    ],
    
    _, (Message[LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrum::badstyle]; Return[$Failed])
  ];
  
	(* Hist for Inverse Functions *)
	
	histInverse =  Show[ ListPlot[
      N[eigvalRecip],
      PlotStyle -> Blue,
      AxesLabel -> {"Re", "Im"},
      AspectRatio -> 1,
      PlotRange -> {{-SumMax - 1, SumMax + 1}, {-SumMax - 1, SumMax + 1}},
      Frame -> True,
      FrameLabel -> {"Re", "Im"},
      GridLines -> Automatic
    ],ImageSize -> 500,
      PlotLabel -> Style[
  Row[{
    "Approximate Spectrum of inverse of: ", 
    TraditionalForm[
      HoldForm[
        Subscript[A, 1] \[CircleTimes] Subscript[u, 1] + 
        Subscript[A, 2] \[CircleTimes] Subscript[u, 2]
      ]
    ]
  }],
  16]
    
    ];
 
  (* Points to plot depending on 2D vs 3D *)
  hist = Switch[plotStyle,
    "3D",
    Show[
      basePlot,
      Graphics3D[{Red, PointSize[Large], Point[{0, SumMax, 0}]}],
      Graphics3D[{Cyan, PointSize[Large], Point[{0, DiffMax, 0}]}],
      Graphics3D[{Green, PointSize[Large], Point[{SumMin, 0, 0}]}],
      Graphics3D[{Purple, PointSize[Large], Point[{DiffMin, 0, 0}]}],
      Graphics3D[{Magenta, PointSize[Large], Point[{-norm1, 0, 0}]}],
      ImageSize -> 500
    ],

    "2D",
    Show[
      basePlot,
      Graphics[{
        Red, PointSize[Large], Point[{0, SumMax}],
        Cyan, PointSize[Large], Point[{0, DiffMax}],
        Green, PointSize[Large], Point[{SumMin,0 }],
        Purple, PointSize[Large], Point[{DiffMin, 0}],
        Magenta, PointSize[Large], Point[{-N[norm1], 0}],
        Orange, PointSize[Medium], Point[Transpose[{ConstantArray[0, Length[AnnuliRadiiJoin]], -AnnuliRadiiJoin}]]
      }],
      PlotLabel -> Style[
  Row[{
    "Approximate Spectrum of: ", 
    TraditionalForm[
      HoldForm[
        Subscript["A", 1] \[CircleTimes] Subscript["u", 1] + 
        Subscript["A", 2] \[CircleTimes] Subscript["u", 2]
      ]
    ], " at size ", Dimensions[U[[1]]]
  }],
  16
],
      ImageSize -> 500
    ]
  ];

  (* Legend as side column *)
  legendLabels = Column[{
    Style["Legend", Bold, 14],
    Style[SumMax ": OuterSpec (spectral rad). " , Red],
    Style[DiffMax ": Inner Rad (larger).", Cyan],
    Style[SumMin ": Outer Rad (smaller).",Green],
    Style[DiffMin ": Inner Rad (smaller).", Purple],
    Style[N[norm1] ": Norm2Vneum.", Magenta],
    Style[AnnuliRadiiJoin , Orange]
  }, Spacings -> 1.5];
Print[];
 Grid[
  {
    {hist, legendLabels}(* ,   (* Top row: histogram + legend *)
    {histInverse} *)   (* Bottom row: histogram spanning both columns *) 
  },
  Spacings -> {2, 1}
]
];

LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrum::badstyle = 
  "PlotStyle must be either \"2D\" or \"3D\".";
(* ------------------------------------------------------------------------------------------------------ *)  
LodewykJansenvanRensburg`ApproximateSpectra`GeneratesFullMatrixAlgebraQ[A_, B_] := Module[
  {gens, words, monomials, matList, rank},

  gens = {IdentityMatrix[2], A, B};

  (* Step 1: Generate all "words" of length 1 to 4 from gens *)
  words = Flatten[
    Table[
      Tuples[gens, len],
      {len, 1, 4}
    ],
    1
  ];

  (* Step 2: Evaluate each word as a matrix product *)
  monomials = Dot @@@ words;

  (* Step 3: Remove duplicates up to numerical tolerance *)
  matList = DeleteDuplicates[monomials, MatrixEqualQ];

  (* Step 4: Flatten matrices into vectors and compute matrix rank *)
  rank = MatrixRank[Flatten /@ matList];

  (* Step 5: Return whether the span is 4D (i.e., full M\:2082) *)
  Return[rank == 4]
];

(* Matrix equality up to numerical precision *)
LodewykJansenvanRensburg`ApproximateSpectra`MatrixEqualQ[A_, B_] := Chop[Norm[A - B]] < 10^-10;

(* A way to Compute the Positive Part*)
LodewykJansenvanRensburg`ApproximateSpectra`PolarDecomposition[A_] := Module[{U, S, V, P},
  {U, S, V} = SingularValueDecomposition[A];
  P = V . DiagonalMatrix[S] . ConjugateTranspose[V]; (* P = sqrt(A\[Dagger]A) *)
  {U . ConjugateTranspose[V], P}
];
  

(* 5. ------------------------- Function that does testing -------------------------- *)
ClearAll[LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrumSingleVar]
Options[LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrumSingleVar] = {PlotStyle -> "3D", NumBins -> 20};

LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrumSingleVar[A1_, U_, opts : OptionsPattern[]] := Module[
  {
    AU, eigval1, Ac1, lambda, hist, 
    SumMax, SumMin, DiffMax, DiffMin,  norm1,
    plotStyle, numBins, legendLabels, basePlot,AnnuliRadii,
    sqrtEigValAU
  },

  plotStyle = OptionValue[PlotStyle];
  numBins = OptionValue[NumBins];

  (* Lambda operator *)
  lambda[t1_, t2_] := Total[
    MapThread[KroneckerProduct[Conjugate[#1], #2] &, {t1, t2}]
  ];

  (* Matrix of interest *)
  Ac1 = Conjugate /@ A1;
  AU = lambda[A1, U];

  (* Spectral radius estimates *)
  SumMax = Abs[EstimateSpectralRadius[A1, "outer",  Verbose -> False, Extremum -> "Max"]];
  DiffMax = Abs[EstimateSpectralRadius[A1, "inner", Verbose -> False, Extremum -> "Max"]];
  SumMin = Abs[EstimateSpectralRadius[A1, "outer", Verbose -> False, Extremum -> "Min"]];
  DiffMin = Abs[EstimateSpectralRadius[A1, "inner", Verbose -> False, Extremum -> "Min"]];
  norm1 = Norm2Vneum[A1];

  (* Eigenvalues and xy projection *)
  eigval1 = {Re[#], Im[#]} & /@ Eigenvalues[AU];
  sqrtEigValAU = Sqrt[
           N[
               Abs[ 
                   Eigenvalues[ 
                                  KroneckerProduct[ Conjugate[A1[[1]]],A1[[1]] ] 
                               ] 
                  ]
            ] 
        ];
  
   (* 10^-8 accuracy*)
    sqrtEigValAU = N[Round[#, 10^-8] & /@ sqrtEigValAU];
	AnnuliRadii = Select[sqrtEigValAU, Count[sqrtEigValAU, #] == 1 &]; (* Only Keeps multiplicity one eig values*)

    Print["==================================================================================================="];
	Print[" Multiplicity Info Tells us When Formulas Work: "];
	Print["===================================================================================================\n"];
  Print["Sqrt of eig values of  Conj(A_1)A_1: ", sqrtEigValAU];
  Print["Sqrt of eig values of  Conj(A_1)A_1 with multiplicity 1: ", AnnuliRadii];
    Print["==================================================================================================="];


  (* Histogram or scatter depending on style *)
  basePlot = Switch[plotStyle,
    "3D",
    Histogram3D[
      eigval1, numBins, "Probability",
      PlotRange -> {{-SumMax - 1, SumMax + 1}, {-SumMax - 1, SumMax + 1}, All},
      ColorFunction -> (Blend[{White, Blue}, #] &),
      ChartLegends -> Automatic
    ],
    
    "2D",
    ListPlot[
      eigval1,
      PlotStyle -> Blue,
      AxesLabel -> {"Re", "Im"},
      AspectRatio -> 1,
      PlotRange -> {{-SumMax - 1, SumMax + 1}, {-SumMax - 1, SumMax + 1}},
      Frame -> True,
      FrameLabel -> {"Re", "Im"},
      GridLines -> Automatic
    ],
    
    _, (Message[LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrumSingleVar::badstyle]; Return[$Failed])
  ];

  (* Points to plot depending on 2D vs 3D *)
  hist = Switch[plotStyle,
    "3D",
    Show[
      basePlot,
      Graphics3D[{Red, PointSize[Large], Point[{0, SumMax, 0}]}],
      Graphics3D[{Cyan, PointSize[Large], Point[{0, DiffMax, 0}]}],
      Graphics3D[{Green, PointSize[Large], Point[{SumMin, 0, 0}]}],
      Graphics3D[{Purple, PointSize[Large], Point[{DiffMin, 0, 0}]}],
      Graphics3D[{Magenta, PointSize[Large], Point[{-norm1, 0, 0}]}],
      ImageSize -> 500
    ],

    "2D",
    Show[
      basePlot,
      Graphics[{
        Red, PointSize[Large], Point[{0, SumMax}],
        Cyan, PointSize[Large], Point[{0, DiffMax}],
        Green, PointSize[Large], Point[{SumMin,0 }],
        Purple, PointSize[Large], Point[{DiffMin, 0}],
        Magenta, PointSize[Large], Point[{-N[norm1], 0}],
        Orange, PointSize[Medium], Point[Transpose[{ConstantArray[0, Length[AnnuliRadii]], -AnnuliRadii}]]
      }],
      ImageSize -> 500
    ]
  ];

  (* Legend as side column *)
  legendLabels = Column[{
    Style["Legend", Bold, 14],
    Style[SumMax ": OuterSpec (spectral rad). " , Red],
    Style[DiffMax ": Inner Rad (larger).", Cyan],
    Style[SumMin ": Outer Rad (smaller).",Green],
    Style[DiffMin ": Inner Rad (smaller).", Purple],
    Style[N[norm1] ": Norm2Vneum.", Magenta],
    Style[AnnuliRadii, Orange]
  }, Spacings -> 1.5];

  Grid[{{hist, legendLabels}}, Spacings -> {2, 1}]
];

LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeTupleSpectrumSingleVar::badstyle = 
  "PlotStyle must be either \"2D\" or \"3D\".";

(* Compute the moment matrix of order n for a tuple of matrices *)
Options[LodewykJansenvanRensburg`ApproximateSpectra`MomentMatrix] = {
  ConjugateOrder -> "StarFirst" (* Options: "Max", "Min", or Automatic *)
};

(* Pass in ConjugateOrder-> "StarFirst" if you want type Aw^*Aw and 
ConjugateOrder-> "StarSecond" if you want type Aw^*Aw .*)

(* Examples: 
	MomentMatrix[A,1, ConjugateOrder->"StarFirst"] 
	MomentMatrix[A,1, ConjugateOrder->"StarSecond"]
*)
LodewykJansenvanRensburg`ApproximateSpectra`MomentMatrix[A_List, n_Integer, opts : OptionsPattern[]] := Module[
  {d, words, momentMat,conjugateOrder},
  
  
  conjugateOrder = OptionValue[ConjugateOrder];
  Print[conjugateOrder, " selected in MomentMatrix"];
  d = Length[A];                         (* Number of matrices in tuple *)
  words = Tuples[Range[d], n];          (* All words of length n *)
 
  If[ conjugateOrder == "StarSecond" ,
  momentMat = Total[
    Table[
      Module[{Aw},
        Aw = Fold[Dot, A[[w[[1]]]], A[[w[[2 ;;]]]]];  (* Construct A_w *)
         Aw . ConjugateTranspose[Aw]                   (* Compute A_w A_w^* *)
      ],
      {w, words}
    ]
  ];
  ,  momentMat = Total[
    Table[
      Module[{Aw},
        Aw = Fold[Dot, A[[w[[1]]]], A[[w[[2 ;;]]]]];  (* Construct A_w *)
         ConjugateTranspose[Aw] . Aw                   (* Compute A_w^* A_w *)
      ],
      {w, words}
    ]
  ];
  ];
  
  momentMat
]

 (* 
(* Compute the moment matrix of order n for a tuple of matrices *)
LodewykJansenvanRensburg`ApproximateSpectra`MomentMatrix[A_List, n_Integer] := Module[
  {d, words, momentMat},
  
  d = Length[A];                         (* Number of matrices in tuple *)
  words = Tuples[Range[d], n];          (* All words of length n *)
  
  momentMat = Total[
    Table[
      Module[{Aw},
        Aw = Fold[Dot, A[[w[[1]]]], A[[w[[2 ;;]]]]];  (* Construct A_w *)
        Aw . ConjugateTranspose[Aw]                   (* Compute A_w A_w^* *)
      ],
      {w, words}
    ]
  ];
  
  momentMat
]
*)
(*
  FindTracePreservingTuple[maxTries_: 10, intRange_: {-5, 5}, verbose_: False]

  Description:
    Attempts to construct a pair of 2\[Times]2 real matrices {A1, A2} that satisfy
    the trace-preserving condition:
    
      A1\[Dagger] A1 + A2\[Dagger] A2 == I  (up to norm constraints)

    The function randomly samples values for entries {a, c, e, f, h}, 
    then solves for {b, d, g} such that:
    
      a^2 + c^2 + e^2 + g^2 == b^2 + d^2 + f^2 + h^2,
      a*b + c*d + e*f + g*h == 0

    If a valid solution is found within the allowed number of attempts, 
    it returns the matrix tuple A = {A1, A2}. Otherwise, it returns $Failed.

  Parameters:
    maxTries   - Maximum number of random attempts (default: 10)
    intRange   - Range from which random integers are drawn (default: {-5, 5})
    verbose    - If True, prints details of attempts and solutions (default: False)

  Returns:
    A list of two 2\[Times]2 matrices {A1, A2}, or $Failed if no solution is found.

  Examples:
    (* Basic usage *)
    A = FindTracePreservingTuple[]

    (* Use wider sampling range *)
    A = FindTracePreservingTuple[20, {-10, 10}]

    (* Verbose output for debugging *)
    A = FindTracePreservingTuple[15, {-3, 3}, True]

    (* Display the result nicely *)
    MatrixForm /@ A
*)

LodewykJansenvanRensburg`ApproximateSpectra`FindTracePreservingTuple[maxTries_: 10, intRange_: {-5, 5}, verbose_: False] := Module[
  {tries = 0, solFound = False, Soladg, a, b, c, d, e, f, g, h, A},

  While[tries < maxTries && !solFound,
    
    (* Random values for fixed variables *)
    Clear[a, b, c, d, e, f, g, h];
    a = RandomInteger[intRange];
    c = RandomInteger[intRange];
    e = RandomInteger[intRange];
    f = RandomInteger[intRange];
    h = RandomInteger[intRange];
    
    If[verbose, Print["Attempt ", tries + 1, ": Random values {a, c, e, f, h} = {", a, ", ", c, ", ", e, ", ", f, ", ", h, "}"]];
    
    (* Solve for b, d, g *)
    Soladg = Quiet@NSolve[
      {
        a^2 + c^2 + e^2 + g^2 == b^2 + d^2 + f^2 + h^2,
        a*b + c*d + e*f + g*h == 0,
        b \[Element] Reals, d \[Element] Reals, g \[Element] Reals
      },
      {b, d, g}
    ];
    
    tries++;

    If[Length[Soladg] > 0,
      b = b /. Soladg[[1, 1]];
      d = d /. Soladg[[1, 2]];
      g = g /. Soladg[[1, 3]];
      
      A = {{{a, b}, {c, d}}, {{e, f}, {g, h}}};
      solFound = True;

      If[verbose,
        Print["\:2705 Solution found on attempt ", tries];
        Print["Matrix A1 = ", MatrixForm[A[[1]]]];
        Print["Matrix A2 = ", MatrixForm[A[[2]]]];
        Print["Trace Preserving Check: A1^*A1 + A2^*A2 = ", 
          MatrixForm[
            ConjugateTranspose[A[[1]]] . A[[1]] + ConjugateTranspose[A[[2]]] . A[[2]]
          ]
        ];
      ];
    ];
  ];

  If[solFound,
    A,
    If[verbose, Print["\:274c No Solution Found in ", maxTries, " attempts."]];
    $Failed
  ]
]   

(* Just a way to print a matrix tuple*)         
LodewykJansenvanRensburg`ApproximateSpectra`PrintTuple[A_List] := Module[{},
  Do[
    Print["A", i, " = ", MatrixForm[A[[i]]]],
    {i, Length[A]}
  ];
]

(* 
  Purpose:
    Computes and plots the variance of the spectral radius of the matrix 
    sum A\:2081\[CircleTimes]U\:2081 + A\:2082\[CircleTimes]U\:2082 where the U\:1d62 are Haar-random unitaries of size d x d.
    The variance is computed empirically for each d from 1 to maxSize using sampleSize realizations.
    Fits the resulting curve to a power law decay model a/d^b + c.
  
  Parameters:
    A: A list of fixed k\[Times]k matrices (typically of length 2), forming the Kronecker coefficients.
    sampleSize: Number of random samples used at each matrix size d.
    maxSize: The maximum size of the d\[Times]d unitaries to test.

  Returns:
    A plot of variance vs d with an overlaid fitted curve and fit exponent displayed.
*)

LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusVarianceFitPlot[A_, sampleSize_,minSize_, maxSize_] := Module[
  {
    g = Length[A], (* Number of matrices in the tuple *)
    HaarUnitary,
    spectralRadiusVariance,
    statsData, varianceData,
    fit, fitEquation, exponentB,
    varPlot, fitPlot
  },

  (* Generate Haar-random unitary of size d *)
  HaarUnitary[d_] := RandomVariate[CircularUnitaryMatrixDistribution[d]];

  (* Compute variance of spectral radius at dimension d *)
  spectralRadiusVariance[d_, trials_: sampleSize] := Module[
    {radii = Table[
       Module[{U, M},
         U = Table[HaarUnitary[d], {g}];
         M = Total[Table[KroneckerProduct[A[[j]], U[[j]]], {j, 1, g}]];
         Max[Abs[Eigenvalues[M]]]
       ],
       {trials}]
    },
    Variance[radii]
  ];

  (* Compute variance data for all d = 1..maxSize *)
  varianceData = Table[{d, spectralRadiusVariance[d]}, {d, minSize, maxSize}];

  (* Plot variance data *)
  varPlot = ListLinePlot[
    varianceData,
    AxesLabel -> {"d", "Variance of Spectral Radius"},
    PlotMarkers -> Automatic,
    PlotStyle -> Red,
    GridLines -> Automatic,
    PlotLabel -> "Variance of Spectral Radius vs Dimension d"
  ];

  (* Fit variance to power law decay model a/d^b + c *)
  fit = FindFit[
    varianceData,
    a/d^b + c,
    {{a, 1}, {b, 1}, {c, 0}},
    d
  ];
  fitEquation = a/d^b + c /. fit;
  exponentB = b /. fit;

  (* Plot fitted curve *)
  fitPlot = Plot[
    fitEquation,
    {d, minSize, maxSize},
    PlotStyle -> {Dashed, Black},
    PlotLegends -> {"Fit: a/d^b + c"},
    PlotRange -> All
  ];

  (* Show combined plot with annotation *)
  Show[
    varPlot,
    fitPlot,
    Epilog -> {
      Inset[
        Style[
          Row[{"Fit exponent b \[TildeTilde] ", NumberForm[exponentB, {3, 2}]}],
          Medium,
          Blue
        ],
        Scaled[{0.6, 0.8}] (* Adjust position as desired *)
      ]
    }
  ]
]

(* 
  SpectralRadiusFitPlot[A_, sampleSize_, maxSize_]
  --------------------------------------------------------------------
  Computes and plots the average spectral radius of the matrix:
     A1 \[CircleTimes] U1 + A2 \[CircleTimes] U2
  where U1, U2 are Haar-random unitary matrices of size d\[Times]d, and A1, A2
  are fixed k\[Times]k Hermitian matrices (passed as a list A = {A1, A2}).

  Parameters:
    - A: List of fixed k\[Times]k matrices {A1, A2}
    - sampleSize: Number of Haar samples per dimension d
    - maxSize: Maximum matrix dimension d to compute up to

  Output:
    - A plot of the average spectral radius vs dimension d

  Dependencies:
    - Requires Wolfram Language 13+ for Haar unitary distribution

*)
LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusFitPlot[A_, sampleSize_,minSize_, maxSize_] := Module[
  {
    g = Length[A], (* Number of unitaries (should be 2) *)
    averageSpectralRadius, HaarUnitary,
    data, fit, fitExpr, exponentB,SumMax
  },

  (* Haar-random unitary *)
  HaarUnitary[d_] := RandomVariate[CircularUnitaryMatrixDistribution[d]];

  (* Compute average spectral radius at dimension d *)
  averageSpectralRadius[d_, trials_: sampleSize] := Module[
    {radii = Table[
       Module[{U, M},
         U = Table[HaarUnitary[d], {g}];
         M = Total[Table[KroneckerProduct[A[[j]], U[[j]]], {j, 1, g}]];
         Max[Abs[Eigenvalues[M]]]
       ],
       {trials}]
    },
    Mean[radii]
  ];

  (* Compute average spectral radius data *)
  data = Table[{d, averageSpectralRadius[d]}, {d, minSize, maxSize}];
	SumMax = N[Abs[EstimateSpectralRadius[A, "outer",  Verbose -> False, Extremum -> "Max"]]];
	

  (* Plot data and fit together *)
  ListLinePlot[
    {
      data,
       Table[{d, SumMax}, {d, minSize, maxSize}]
    },
    PlotStyle -> {Blue, {Thin, Black}},
    PlotMarkers -> Automatic,
    AxesLabel -> {"d", "Average Spectral Radius"},
    GridLines -> Automatic,
    PlotLabel -> "Average Spectral Radius with Fit",
    PlotLegends -> {
      "Average Spectral Radius at size dxd",
      Row[{"Spectral Radius of Limiting Operator: ", NumberForm[SumMax]}]
    }
  ]
]; 

LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusHistogramPlot[A_, sampleSize_, d_] := Module[
  {
    g = Length[A],
    HaarUnitary,
    radii, SumMax,
    hist, binEdges, counts, maxCount,
    minBin, maxBin, buffer = 0.01
  },

  HaarUnitary[n_] := RandomVariate[CircularUnitaryMatrixDistribution[n]];

  radii = Table[
    Module[{U, M},
      U = Table[HaarUnitary[d], {g}];
      M = Total[Table[KroneckerProduct[A[[j]], U[[j]]], {j, 1, g}]];
      Max[Abs[Eigenvalues[M]]]
    ],
    {sampleSize}
  ];

  SumMax = N[Abs[EstimateSpectralRadius[A, "outer", Verbose -> False, Extremum -> "Max"]]];

  (* Create histogram bins manually *)
  hist = HistogramList[radii, 20]; (* { {bin edges}, {counts} } *)
  binEdges = hist[[1]];
  counts = hist[[2]];
  maxCount = Max[counts];

  minBin = Min[Min[binEdges], SumMax] - buffer;
  maxBin = Max[Max[binEdges], SumMax] + buffer;

  Show[
    Histogram[radii, {binEdges},
      "Count",
      ChartStyle -> LightBlue,
      AxesLabel -> {"Spectral Radius", "Count"},
      PlotLabel -> Row[{
        "Histogram of Spectral Radii at d = ", d,
        " (n = ", sampleSize, ")"
      }],
      PlotRange -> {{minBin, maxBin}, {0, maxCount * 1.1}},
      ImageSize -> Large
    ],
    Graphics[{Red, Dashed, Thick, Line[{{SumMax, 0}, {SumMax, maxCount * 1.1}}]}],
    Epilog -> {
      Text[Style[
        Row[{"SumMax \[TildeTilde] ", NumberForm[SumMax, {5, 4}]}], Red, 12],
        {SumMax, maxCount * 1.05}, {-1, 0}
      ]}
  ]
];
  LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusHistogramPlotWithFit[A_, sampleSize_, d_] := Module[  {
    g = Length[A],
    HaarUnitary,
    radii, SumMax,
    hist, binEdges, counts, maxCount,
    minBin, maxBin, buffer = 0.01,
    distFit, pdfPoints, binWidth, pdfScale,
    labelX, labelYStart, labelSpacing, labels
  },

  (* Haar-random unitary *)
  HaarUnitary[n_] := RandomVariate[CircularUnitaryMatrixDistribution[n]];

  (* Spectral radius samples *)
  radii = Table[
    Module[{U, M},
      U = Table[HaarUnitary[d], {g}];
      M = Total[Table[KroneckerProduct[A[[j]], U[[j]]], {j, 1, g}]];
      Max[Abs[Eigenvalues[M]]]
    ],
    {sampleSize}
  ];

  (* Limiting operator spectral radius *)
  SumMax = N[Abs[EstimateSpectralRadius[A, "outer", Verbose -> False, Extremum -> "Max"]]];

  (* Histogram and bin settings *)
  hist = HistogramList[radii, 20];
  binEdges = hist[[1]];
  counts = hist[[2]];
  maxCount = Max[counts];
  binWidth = binEdges[[2]] - binEdges[[1]];

  minBin = Min[Min[binEdges], SumMax] - buffer;
  maxBin = Max[Max[binEdges], SumMax] + buffer;

  (* Distribution fit *)
  distFit = Quiet@FindDistribution[radii];
  pdfScale = sampleSize * binWidth;

  (* PDF curve points *)
  pdfPoints = Table[
    {x, pdfScale * PDF[distFit, x]},
    {x, minBin, maxBin, (maxBin - minBin)/300}
  ];

  (* Label positioning *)
  labelX = maxBin - 0.02 (maxBin - minBin);
  labelYStart = maxCount * 1.05;
  labelSpacing = maxCount * 0.05;

  labels = {
    Inset[
      Style[Row[{"Fit: ", distFit}], Orange, 12, Italic],
      {labelX, labelYStart}, {1, 0}
    ],
    Inset[
      Style[Row[{"SumMax \[TildeTilde] ", NumberForm[SumMax, {5, 4}]}], Red, 12],
      {labelX, labelYStart - labelSpacing}, {1, 0}
    ]
  };

  (* Final plot *)
  Show[
    {
      Histogram[radii, {binEdges}, "Count",
        ChartStyle -> LightBlue,
        AxesLabel -> {"Spectral Radius", "Count"},
        PlotLabel -> Row[{
          "Spectral Radii at d = ", d,
          " (n = ", sampleSize, ")"
        }],
        PlotRange -> {{minBin, maxBin}, {0, maxCount * 1.2}},
        ImageSize -> Large
      ],
      Graphics[
        {
          Orange, Thick, Line[pdfPoints],
          Red, Dashed, Thick, Line[{{SumMax, 0}, {SumMax, maxCount * 1.1}}],
          labels
        }
      ]
    }
  ]
];

(* This is just for the Gausian Random Matrices *)

LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusHistogramPlot[A_, sampleSize_, d_] := Module[
  {
    g = Length[A],
    HaarUnitary,
    radii, SumMax,
    hist, binEdges, counts, maxCount,
    minBin, maxBin, buffer = 0.01
  },

  HaarUnitary[n_] := RandomVariate[CircularUnitaryMatrixDistribution[n]];

  radii = Table[
    Module[{U, M},
      U = Table[HaarUnitary[d], {g}];
      M = Total[Table[KroneckerProduct[A[[j]], U[[j]]], {j, 1, g}]];
      Max[Abs[Eigenvalues[M]]]
    ],
    {sampleSize}
  ];

  SumMax = N[Abs[EstimateSpectralRadius[A, "outer", Verbose -> False, Extremum -> "Max"]]];

  (* Create histogram bins manually *)
  hist = HistogramList[radii, 20]; (* { {bin edges}, {counts} } *)
  binEdges = hist[[1]];
  counts = hist[[2]];
  maxCount = Max[counts];

  minBin = Min[Min[binEdges], SumMax] - buffer;
  maxBin = Max[Max[binEdges], SumMax] + buffer;

  Show[
    Histogram[radii, {binEdges},
      "Count",
      ChartStyle -> LightBlue,
      AxesLabel -> {"Spectral Radius", "Count"},
      PlotLabel -> Row[{
        "Histogram of Spectral Radii at d = ", d,
        " (n = ", sampleSize, ")"
      }],
      PlotRange -> {{minBin, maxBin}, {0, maxCount * 1.1}},
      ImageSize -> Large
    ],
    Graphics[{Red, Dashed, Thick, Line[{{SumMax, 0}, {SumMax, maxCount * 1.1}}]}],
    Epilog -> {
      Text[Style[
        Row[{"SumMax \[TildeTilde] ", NumberForm[SumMax, {5, 4}]}], Red, 12],
        {SumMax, maxCount * 1.05}, {-1, 0}
      ]}
  ]
];
  LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusHistogramPlotWithFitGUE[A_, sampleSize_, d_] := Module[  {
    g = Length[A],
    radii, SumMax,
    hist, binEdges, counts, maxCount,
    minBin, maxBin, buffer = 0.01,
    distFit, pdfPoints, binWidth, pdfScale,
    labelX, labelYStart, labelSpacing, labels,GOutside,approxSpecRad
  },



  (* Spectral radius samples *)
  radii = Table[
    Module[{G, M},
      G = Table[(1/Sqrt[d])*RandomVariate[GaussianUnitaryMatrixDistribution[d]], {g}];
      M = Total[Table[KroneckerProduct[A[[j]], G[[j]]], {j, 1, g}]];
      Max[Abs[Eigenvalues[M]]]
    ],
    {sampleSize}
  ];
	GOutside = Table[(1/Sqrt[1000])*RandomVariate[GaussianUnitaryMatrixDistribution[1000]], {g}];
	approxSpecRad = Max[Abs[Eigenvalues[Total[Table[KroneckerProduct[A[[j]], GOutside[[j]]], {j, 1, g}]]] ]];
  (* Limiting operator spectral radius *)
  SumMax = N[Abs[EstimateSpectralRadius[A, "outer", Verbose -> False, Extremum -> "Max"]]];

  (* Histogram and bin settings *)
  hist = HistogramList[radii, 20];
  binEdges = hist[[1]];
  counts = hist[[2]];
  maxCount = Max[counts];
  binWidth = binEdges[[2]] - binEdges[[1]];

  minBin = Min[Min[binEdges], SumMax] - buffer;
  maxBin = Max[Max[binEdges], SumMax] + buffer;

  (* Distribution fit *)
  distFit = Quiet@FindDistribution[radii];
  pdfScale = sampleSize * binWidth;

  (* PDF curve points *)
  pdfPoints = Table[
    {x, pdfScale * PDF[distFit, x]},
    {x, minBin, maxBin, (maxBin - minBin)/300}
  ];

  (* Label positioning *)
  labelX = maxBin - 0.02 (maxBin - minBin);
  labelYStart = maxCount * 1.05;
  labelSpacing = maxCount * 0.05;

  labels = {
    Inset[
      Style[Row[{"Fit: ", distFit}], Orange, 12, Italic],
      {labelX, labelYStart}, {1, 0}
    ],
    Inset[
      Style[Row[{"Approx Mean \[TildeTilde] ", NumberForm[approxSpecRad, {5, 4}]}], Red, 12],
      {labelX, labelYStart - labelSpacing}, {1, 0}
    ]
  };

  (* Final plot *)
  Show[
    {
      Histogram[radii, {binEdges}, "Count",
        ChartStyle -> LightBlue,
        AxesLabel -> {"Spectral Radius", "Count"},
        PlotLabel -> Row[{
          "Spectral Radii at d = ", d,
          " (n = ", sampleSize, ")"
        }],
        PlotRange -> {{minBin, maxBin}, {0, maxCount * 1.2}},
        ImageSize -> Large
      ],
      Graphics[
        {
          Orange, Thick, Line[pdfPoints],
          Red, Dashed, Thick, Line[{{approxSpecRad, 0}, {approxSpecRad, maxCount * 1.1}}],
          labels
        }
      ]
    }
  ]
];

(* 
SpectralRadiusConvExpectation[A_, sampleSize_, minSize_, maxSize_]

Purpose:
  Computes and visualizes the convergence behavior of the spectral radius of 
  random matrix ensembles of the form \[CapitalSigma] A_j \[CircleTimes] U_j, where:
    - A_j are fixed coefficient matrices,
    - U_j are independently sampled Haar-random unitary matrices.

Inputs:
  A         \[LongDash] A list of coefficient matrices {A\:2081, A\:2082, ..., A_g}, where each A_j is a d\[Times]d matrix.
  sampleSize \[LongDash] Number of random realizations to average over at each matrix size.
  minSize    \[LongDash] Minimum matrix dimension (for U_j) to evaluate.
  maxSize    \[LongDash] Maximum matrix dimension (for U_j) to evaluate.

Behavior:
  - For each dimension d in [minSize, maxSize], samples `sampleSize` instances of 
    random matrices M_d = \[CapitalSigma] A_j \[CircleTimes] U_j, where each U_j is d\[Times]d Haar unitary.
  - Computes the average absolute difference between the spectral radius of M_d 
    and the estimated limiting spectral radius of the operator associated with A.
  - Plots the convergence curve E[|rad_d - rad|] as a function of d.

Output:
  A ListLinePlot showing convergence of the average spectral radius to its 
  limiting value, along with descriptive legends.

Dependencies:
  - Requires `EstimateSpectralRadius` function (assumed defined elsewhere).
  - Uses `CircularUnitaryMatrixDistribution` for Haar-random unitary sampling.
*)
LodewykJansenvanRensburg`ApproximateSpectra`SpectralRadiusConvExpectation[A_, sampleSize_,minSize_, maxSize_] := Module[
  {
    g = Length[A], (* Number of unitaries (should be 2) *)
    averageDifference, HaarUnitary,
    data, fit, fitExpr, exponentB,SumMax,labeledData
  },

  (* Haar-random unitary *)
  HaarUnitary[d_] := RandomVariate[CircularUnitaryMatrixDistribution[d]];
	SumMax = N[Abs[EstimateSpectralRadius[A, "outer",  Verbose -> False, Extremum -> "Max"]]];
  (* Compute average spectral radius at dimension d *)
  averageDifference[d_, trials_: sampleSize] := Module[
    {radii = Table[
       Module[{U, M},
         U = Table[HaarUnitary[d], {g}];
         M = Total[Table[KroneckerProduct[A[[j]], U[[j]]], {j, 1, g}]];
         Abs[SumMax-Max[Abs[Eigenvalues[M]]]]
       ],
       {trials}]
    },
    Mean[radii]
  ];

  (* Compute average spectral radius data *)
  data = Table[{d, averageDifference[d]}, {d, minSize, maxSize}];

	
  (* Wrap data with custom legend using Legended *)
  labeledData = {
    Legended[
      data,
      Column[{
        "Average Error to Limiting Spectral Radius",
        Row[{"Spectral Radius of Limiting Operator: ", NumberForm[SumMax]}],
        Row[{"Number of Realizations: ", sampleSize}]
      }]
    ]
  };

  (* Plot *)
  ListLinePlot[
    labeledData,
    PlotStyle -> Blue,
    PlotMarkers -> Automatic,
    AxesLabel -> {"d", "Average Error"},
    GridLines -> Automatic,
    PlotLabel -> Row[{"E[| rad_", Style["d", Italic], " - rad |]"}]
  ]
];

(*
  AnalyzeCharPoly

  Description:
    Given a list of two square matrices A = {A1, A2}, this function computes and analyzes
    the characteristic polynomials of the matrices:
      - AiSumAiBar = Conjugate[A1] \[CircleTimes] A1 + Conjugate[A2] \[CircleTimes] A2
      - AiMinusAiBar = Conjugate[A1] \[CircleTimes] A1 - Conjugate[A2] \[CircleTimes] A2

    It calculates their characteristic polynomials, roots, multiplicities, and
    plots these polynomials together with their roots. It also produces tables
    summarizing roots and multiplicities.

  Inputs:
    - A : List of two square matrices {A1, A2}
    - showExtraInfo (optional, default True) : Boolean controlling whether to
      print additional detailed info and extra plots.

  Outputs:
    - Prints combined plot of polynomials and roots with tables side-by-side.
    - Optionally prints detailed polynomial and matrix info with additional plots.
    - Returns an Association containing key computed objects for further use.

  Usage:
    AnalyzeCharPoly[{A1, A2}, True]  (* with extra info *)
    AnalyzeCharPoly[{A1, A2}, False] (* without extra info *)
*)
LodewykJansenvanRensburg`ApproximateSpectra`AnalyzeCharPoly[A_List, showExtraInfo_: True] := Module[
  {
    AiSumAiBar, AiMinusAiBar,
    SpecRad,
    charPoly, charPolyMinus,
    roots, rootsMinus,
    mults, multsMinus,
    multsWithSqrt, multsWithSqrtMinus,
    display, displayMinus,
    RootTableSum, RootTableMins,
    RootsSumPlotItem, RootsMinusPlotItem,
    zmin, zmax,
    plot
  },
  
  Print["============================================================================================"];
  Print["Analyzing Characteristic Polynomial of A1bar \[CircleTimes] A1 + A2bar \[CircleTimes] A2 and A1bar \[CircleTimes] A1 - A2bar \[CircleTimes] A2 : "];
  Print["==========================================================================================\n"];
  
  (* Compute AiSumAiBar and AiMinusAiBar *)
  AiSumAiBar = KroneckerProduct[Conjugate[A[[1]]], A[[1]]] + KroneckerProduct[Conjugate[A[[2]]], A[[2]]];
  AiMinusAiBar = KroneckerProduct[Conjugate[A[[1]]], A[[1]]] - KroneckerProduct[Conjugate[A[[2]]], A[[2]]];
  
  (* Spectral radius *)
  SpecRad = N[Max[Abs[Eigenvalues[AiSumAiBar]]]];
  
  (* Characteristic polynomials *)
  charPoly = CharacteristicPolynomial[AiSumAiBar, z];
  charPolyMinus = CharacteristicPolynomial[AiMinusAiBar, z];
  
  (* Roots and multiplicities *)
  roots = z /. NSolve[charPoly == 0, z];
  rootsMinus = z /. NSolve[charPolyMinus == 0, z];
  
  mults = Tally[roots];
  multsMinus = Tally[rootsMinus];
  
  multsWithSqrt = mults /. {r_?NumericQ, m_} :> {r, m, Sqrt[Abs[r]]};
  multsWithSqrtMinus = multsMinus /. {r_?NumericQ, m_} :> {r, m, Sqrt[Abs[r]]};
  
  display = multsWithSqrt /. {r_, m_, s_} :> {N[r, 12], m, N[s, 12]};
  displayMinus = multsWithSqrtMinus /. {r_, m_, s_} :> {N[r, 12], m, N[s, 12]};
  
  (* Tables *)
  RootTableSum = Style[
    TableForm[display, TableHeadings -> {None, {"Root", "Multp", "Sqrt(Abs of Root)"}}], 
    FontSize -> 10
  ];
  
  RootTableMins = Style[
    TableForm[displayMinus, TableHeadings -> {None, {"Root", "Multp", "Sqrt(Abs of Root)"}}],
    FontSize -> 10
  ];
  
  (* Roots plots *)
  RootsSumPlotItem = ListPlot[
    Table[{Re[r], 0}, {r, roots}],
    PlotStyle -> {Purple, PointSize[Large]}
  ];
  
  RootsMinusPlotItem = ListPlot[
    Table[{Re[r], 0}, {r, rootsMinus}],
    PlotStyle -> {Red, PointSize[Large]}
  ];
  
  (* Plot range *)
  zmin = -SpecRad - 0.5;
  zmax = SpecRad + 0.5;
  
  (* Combined polynomial plot *)
  plot = Plot[
    {charPoly, charPolyMinus},
    {z, zmin, zmax},
    PlotRange -> All,
    PlotStyle -> {Blue, Green},
    AxesLabel -> {"z", "P(z)"},
    PlotLegends -> Placed[{"sumAibarAi", "minAibarAi"}, Above],
    PlotLabels -> None
  ];
  
  (* Show combined plot with roots, and tables side by side *)
  Print[
    Row[{
      Show[plot, RootsSumPlotItem, RootsMinusPlotItem, ImageSize -> 800],
      Column[{RootTableSum, RootTableMins}]
    }]
  ];
  
  (* Optional extra info *)
  If[TrueQ[showExtraInfo],
    Print["========= Additional Info ============= "];
    
    Print["Characteristic Polynomial of A1bar \[CircleTimes] A1 + A2bar \[CircleTimes] A2: ", charPoly];
    Print["Matrix: A1bar \[CircleTimes] A1 + A2bar \[CircleTimes] A2 = ", MatrixForm[AiSumAiBar]];
    
    If[VectorQ[CoefficientList[charPoly, z], Im[#] == 0 &],
      Show[
        Plot[charPoly, {z, -SpecRad - 0.05, SpecRad + 0.05}, PlotRange -> All,
          AxesLabel -> {"z", "P(z)"},
          PlotStyle -> Thick
        ],
        RootsSumPlotItem
      ],
      ListPlot[
        Table[{Re[r], Im[r]}, {r, roots}],
        PlotStyle -> {Red, PointSize[Large]},
        AxesLabel -> {"Re", "Im"},
        AspectRatio -> 1,
        GridLines -> Automatic
      ]
    ];
    
    Print["\nCharacteristic Polynomial of A1bar \[CircleTimes] A1 - A2bar \[CircleTimes] A2: ", charPolyMinus];
    Print["Matrix: A1bar \[CircleTimes] A1 - A2bar \[CircleTimes] A2 = ", MatrixForm[AiMinusAiBar]];
    
    If[VectorQ[CoefficientList[charPolyMinus, z], Im[#] == 0 &],
      Show[
        Plot[charPolyMinus, {z, -SpecRad - 0.05, SpecRad + 0.05}, PlotRange -> All,
          AxesLabel -> {"z", "P(z)"},
          PlotStyle -> Thick
        ],
        RootsMinusPlotItem
      ],
      ListPlot[
        Table[{Re[r], Im[r]}, {r, rootsMinus}],
        PlotStyle -> {Red, PointSize[Large]},
        AxesLabel -> {"Re", "Im"},
        AspectRatio -> 1,
        GridLines -> Automatic
      ]
    ];
    
    Print["========= End ============= "];
  ];
  
  (* Return main combined output as result (optional) *)
  <|
    "AiSumAiBar" -> AiSumAiBar,
    "AiMinusAiBar" -> AiMinusAiBar,
    "CharPolySum" -> charPoly,
    "CharPolyMinus" -> charPolyMinus,
    "RootsSum" -> roots,
    "RootsMinus" -> rootsMinus,
    "RootTableSum" -> RootTableSum,
    "RootTableMins" -> RootTableMins,
    "Plot" -> plot
  |>
]


End[] (* End `Private` *)

EndPackage[]
