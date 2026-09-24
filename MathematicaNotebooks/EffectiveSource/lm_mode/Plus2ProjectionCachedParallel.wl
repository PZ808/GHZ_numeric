(* ::Package:: *)

(*
  Plus2ProjectionCachedParallel.wl

  Clean projection layer for the spin +2 Teukolsky effective source.

  Expected infrastructure already defined in the parent notebook:

    centralRule
      Association with keys "yGrid" and "Weights".  The entries are the
      angular quadrature nodes y = cos(theta) and their dy weights.

    omegaOfM[m]
      Numerical circular-orbit mode frequency.

    sphValsOnRule[s, m, omega, ellList, rule]
      Cache-aware spheroidal-table accessor.  It must return an Association
      mapping ell -> values on rule["yGrid"].  If the disk-cache accessor has
      a different name, pass it through the option "SphValuesProvider".

    plus2TfunOfM[m]
      Builder returning a numerical source function Tfun[X, Y].  The second
      coordinate is the mixed coordinate Y = -r0val y.

  Main speedups relative to the earlier notebook implementation:

    1. Cached spheroidal values are loaded once per m as a single ell-by-y
       matrix.  No spheroidal functions are built on auxiliary kernels.

    2. The Teukolsky source is sampled once at each (X,Y) point.  All ell
       projections are then formed simultaneously by matrix multiplication.
       The previous ell-task layout recomputed the same source once per ell.

    3. Parallel tasks are radial chunks, not individual ell modes.

    4. Each completed m mode is atomically saved, so interrupted runs resume.
*)

(* ---------------------------------------------------------------------- *)
(* 0. Generic utilities                                                    *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  plus2EnsureDirectory,
  plus2ShortHash,
  plus2NumberTag,
  plus2SpinTag,
  plus2RuleTag,
  plus2ParameterTag,
  plus2PatchTag,
  plus2ConventionTag,
  plus2RunTag,
  plus2DefaultCacheRoot,
  plus2ProjectionDirectory,
  plus2ProjectionMFile,
  plus2AtomicDumpSave,
  plus2ReadMXObject,
  $Plus2MXObject
];

plus2EnsureDirectory[dir_String] :=
  If[! DirectoryQ[dir],
    CreateDirectory[dir, CreateIntermediateDirectories -> True]
  ];

plus2ShortHash[expr_] := Module[{hex},
  hex = IntegerString[Hash[expr, "SHA256"], 16];
  StringTake[hex, Min[16, StringLength[hex]]]
];

plus2NumberTag[x_?NumericQ] := Module[{s},
  s = ToString @ NumberForm[
    Chop[N[x, 16]],
    {Infinity, 12},
    NumberPoint -> ".",
    NumberPadding -> {"", ""},
    ExponentFunction -> (Null &)
  ];
  s = StringTrim[s];
  If[StringContainsQ[s, "."],
    s = StringReplace[s, RegularExpression["0+$"] -> ""];
    s = StringReplace[s, RegularExpression["\\.$"] -> ""];
  ];
  If[s === "-0", s = "0"];
  StringReplace[s, {"-" -> "m", "." -> "p"}]
];

plus2SpinTag[] := Module[{aa},
  If[! NumericQ[aKerr],
    Print["aKerr must be numerical before the projection cache is used."];
    Abort[];
  ];
  aa = Chop[N[aKerr, 16]];
  If[TrueQ[aa == 0], "a0p0", "a" <> plus2NumberTag[aa]]
];

plus2RuleTag[rule_Association] := Module[{y, w},
  y = N[rule["yGrid"], MachinePrecision];
  w = N[rule["Weights"], MachinePrecision];
  "grid" <> plus2ShortHash[{y, w}]
];

plus2ParameterTag[] := Module[{pars},
  pars = {
    Quiet @ Check[N[aKerr, 16], Missing["aKerr"]],
    Quiet @ Check[N[MKerr, 16], Missing["MKerr"]],
    Quiet @ Check[N[r0val, 16], Missing["r0val"]],
    Quiet @ Check[N[\[ScriptCapitalB]0val, 16], Missing["B0val"]]
  };
  "params" <> plus2ShortHash[pars]
];

plus2PatchTag[patchAssoc_Association] :=
  "xgrid" <> plus2ShortHash[
    Normal @ Map[N[#, MachinePrecision] &, patchAssoc]
  ];

plus2ConventionTag[
  basisSpin_Integer,
  phiFactor_,
  conjugateBasis_
] := "conv" <> plus2ShortHash[{
  basisSpin,
  N[phiFactor, 16],
  TrueQ[conjugateBasis]
}];

plus2RunTag[
  patchAssoc_Association,
  basisSpin_Integer,
  phiFactor_,
  conjugateBasis_
] := StringRiffle[{
  plus2PatchTag[patchAssoc],
  plus2ConventionTag[basisSpin, phiFactor, conjugateBasis]
}, "_"];

plus2DefaultCacheRoot[] := Quiet @ Check[
  FileNameJoin[{NotebookDirectory[], "TeffProjectionCache"}],
  FileNameJoin[{$TemporaryDirectory, "TeffProjectionCache"}]
];

plus2ProjectionDirectory[
  root_String,
  ellMax_Integer,
  rule_Association,
  runTag_String
] := FileNameJoin[{
  root,
  "Plus2Projections_" <> plus2SpinTag[] <> "_" <> plus2ParameterTag[],
  "ellMax" <> ToString[ellMax],
  plus2RuleTag[rule],
  runTag
}];

plus2ProjectionMFile[
  root_String,
  ellMax_Integer,
  rule_Association,
  runTag_String,
  m_Integer
] := FileNameJoin[{
  plus2ProjectionDirectory[root, ellMax, rule, runTag],
  "m" <> ToString[m] <> ".mx"
}];

plus2AtomicDumpSave[file_String, obj_] := Module[{dir, tmp},
  dir = DirectoryName[file];
  plus2EnsureDirectory[dir];
  tmp = FileNameJoin[{
    dir,
    "." <> FileBaseName[file] <> "-" <> CreateUUID[] <> ".mx"
  }];
  Clear[$Plus2MXObject];
  $Plus2MXObject = obj;
  DumpSave[tmp, $Plus2MXObject];
  If[FileExistsQ[file], DeleteFile[file]];
  RenameFile[tmp, file];
  Clear[$Plus2MXObject];
  file
];

plus2ReadMXObject[file_String] := Module[{},
  If[! FileExistsQ[file], Return[$Failed]];
  Clear[$Plus2MXObject];
  Get[file];
  $Plus2MXObject
];

(* ---------------------------------------------------------------------- *)
(* 1. Cached angular projector                                             *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  $Plus2AngularProjectorMemory,
  clearPlus2AngularProjectorMemory,
  plus2EllList,
  plus2ResolveSphValuesProvider,
  plus2DefaultSphValuesProvider,
  plus2BuildAngularProjector,
  plus2CachedAngularProjector,
  plus2CheckAngularProjector
];

$Plus2AngularProjectorMemory = <||>;

clearPlus2AngularProjectorMemory[] :=
  ($Plus2AngularProjectorMemory = <||>;);

plus2EllList[sourceSpin_Integer, ellMax_Integer, m_Integer] :=
  Range[Max[Abs[sourceSpin], Abs[m]], ellMax];

plus2ResolveSphValuesProvider[Automatic] := plus2DefaultSphValuesProvider;
plus2ResolveSphValuesProvider[f_] := f;

plus2DefaultSphValuesProvider[
  s_Integer,
  m_Integer,
  omega_?NumericQ,
  ellList_List,
  rule_Association
] := Module[{vals},
  If[DownValues[sphValsOnRule] === {},
    Print[
      "No definition was found for sphValsOnRule.  Load the spheroidal " <>
      "table-cache section first, or pass its accessor through " <>
      "\"SphValuesProvider\"."
    ];
    Abort[];
  ];
  vals = sphValsOnRule[s, m, omega, ellList, rule];
  If[! AssociationQ[vals],
    Print["The spheroidal value provider did not return an Association."];
    Print[Short[vals, 3]];
    Abort[];
  ];
  vals
];

Options[plus2BuildAngularProjector] = {
  "BasisSpin" -> -2,
  "PhiFactor" -> 2 Pi,
  "ConjugateBasis" -> True,
  "MixedYMap" -> Automatic,
  "SphValuesProvider" -> Automatic
};

plus2BuildAngularProjector[
  m_Integer,
  ellMax_Integer,
  rule_Association,
  omega_?NumericQ,
  OptionsPattern[]
] := Module[
  {
    basisSpin, phiFactor, conjugateBasis, mixedYMap, provider,
    ellList, yGrid, weights, YGrid, sphAssoc, sphMat, bra, projMat
  },

  basisSpin = OptionValue["BasisSpin"];
  phiFactor = OptionValue["PhiFactor"];
  conjugateBasis = TrueQ[OptionValue["ConjugateBasis"]];
  provider = plus2ResolveSphValuesProvider[OptionValue["SphValuesProvider"]];

  ellList = plus2EllList[+2, ellMax, m];
  yGrid = Developer`ToPackedArray @ N[rule["yGrid"], MachinePrecision];
  weights = Developer`ToPackedArray @ N[rule["Weights"], MachinePrecision];

  If[Length[yGrid] =!= Length[weights],
    Print["The angular quadrature rule has inconsistent node and weight lengths."];
    Abort[];
  ];

  mixedYMap = Replace[
    OptionValue["MixedYMap"],
    Automatic :> Function[{yy}, N[-r0val yy, MachinePrecision]]
  ];
  YGrid = Developer`ToPackedArray @ N[mixedYMap /@ yGrid, MachinePrecision];

  sphAssoc = provider[basisSpin, m, omega, ellList, rule];
  If[! And @@ (KeyExistsQ[sphAssoc, #] & /@ ellList),
    Print["The cached spheroidal table is missing one or more requested ell modes."];
    Print["Requested ell list: ", ellList];
    Abort[];
  ];

  sphMat = Developer`ToPackedArray @ N[Lookup[sphAssoc, ellList], MachinePrecision];

  If[Dimensions[sphMat] =!= {Length[ellList], Length[yGrid]},
    Print["Unexpected cached spheroidal-table dimensions: ", Dimensions[sphMat]];
    Print["Expected: ", {Length[ellList], Length[yGrid]}];
    Abort[];
  ];

  bra = If[conjugateBasis, Conjugate[sphMat], sphMat];
  projMat = Developer`ToPackedArray @ N[
    phiFactor Map[# weights &, bra],
    MachinePrecision
  ];

  <|
    "SourceSpin" -> +2,
    "BasisSpin" -> basisSpin,
    "m" -> m,
    "omega" -> omega,
    "ellList" -> ellList,
    "yGrid" -> yGrid,
    "YGrid" -> YGrid,
    "Weights" -> weights,
    "PhiFactor" -> phiFactor,
    "ConjugateBasis" -> conjugateBasis,
    "BasisMatrix" -> sphMat,
    "ProjectionMatrix" -> projMat,
    "RuleTag" -> plus2RuleTag[rule]
  |>
];

Options[plus2CachedAngularProjector] = Options[plus2BuildAngularProjector];

plus2CachedAngularProjector[
  m_Integer,
  ellMax_Integer,
  rule_Association,
  omega_?NumericQ,
  OptionsPattern[]
] := Module[{key, projector, buildOptions},
  buildOptions = {
    "BasisSpin" -> OptionValue["BasisSpin"],
    "PhiFactor" -> OptionValue["PhiFactor"],
    "ConjugateBasis" -> OptionValue["ConjugateBasis"],
    "MixedYMap" -> OptionValue["MixedYMap"],
    "SphValuesProvider" -> OptionValue["SphValuesProvider"]
  };
  key = plus2ShortHash[{
    m,
    ellMax,
    N[omega, 16],
    plus2RuleTag[rule],
    buildOptions
  }];
  If[KeyExistsQ[$Plus2AngularProjectorMemory, key],
    Return[$Plus2AngularProjectorMemory[key]]
  ];
  projector = plus2BuildAngularProjector[
    m,
    ellMax,
    rule,
    omega,
    Sequence @@ buildOptions
  ];
  $Plus2AngularProjectorMemory[key] = projector;
  projector
];

plus2CheckAngularProjector[projector_Association] := Module[{gram, nEll},
  nEll = Length[projector["ellList"]];
  gram = projector["ProjectionMatrix"] . Transpose[projector["BasisMatrix"]];
  <|
    "ellList" -> projector["ellList"],
    "Dimensions" -> Dimensions[projector["ProjectionMatrix"]],
    "MaxIdentityError" -> Max[Abs[gram - IdentityMatrix[nEll]]],
    "GramMatrix" -> gram
  |>
];

(* ---------------------------------------------------------------------- *)
(* 2. Parallel X-chunk projection                                          *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  $Plus2CurrentSource,
  $Plus2CurrentYGrid,
  $Plus2CurrentProjectionMatrix,
  plus2ProjectCurrentSourceAtX,
  plus2ProjectCurrentChunk,
  plus2ChooseChunkSize,
  plus2MakeChunkTasks,
  plus2PrepareParallelState,
  plus2ClearParallelState,
  plus2AssemblePatchMatrices,
  plus2LegacyDataView
];

Clear[
  $Plus2CurrentSource,
  $Plus2CurrentYGrid,
  $Plus2CurrentProjectionMatrix
];

plus2ProjectCurrentSourceAtX[x_?NumericQ] := Module[{vals},
  vals = Developer`ToPackedArray @ Table[
    $Plus2CurrentSource[x, $Plus2CurrentYGrid[[j]]],
    {j, Length[$Plus2CurrentYGrid]}
  ];
  If[! VectorQ[vals, NumericQ],
    Print["Non-numeric source values on kernel ", $KernelID, " at X = ", x];
    Print[Short[vals, 3]];
    Abort[];
  ];
  Developer`ToPackedArray[$Plus2CurrentProjectionMatrix . vals]
];

plus2ProjectCurrentChunk[task_Association] := Module[{xGrid, values},
  xGrid = task["XGrid"];
  values = Developer`ToPackedArray @ Table[
    plus2ProjectCurrentSourceAtX[x],
    {x, xGrid}
  ];
  <|
    "Patch" -> task["Patch"],
    "Indices" -> task["Indices"],
    "Values" -> values
  |>
];

plus2ChooseChunkSize[patchAssoc_Association, Automatic] := Module[
  {nXMax, nKernels},
  LaunchKernels[];
  nXMax = Max[Length /@ Values[patchAssoc]];
  nKernels = Max[1, Length[Kernels[]]];
  Max[1, Ceiling[nXMax/(4 nKernels)]]
];

plus2ChooseChunkSize[_Association, n_Integer?Positive] := n;

plus2MakeChunkTasks[patchAssoc_Association, chunkSize_] := Module[{size},
  size = plus2ChooseChunkSize[patchAssoc, chunkSize];
  Flatten @ KeyValueMap[
    Function[{patchName, xGrid},
      Map[
        Function[idx,
          <|
            "Patch" -> patchName,
            "Indices" -> idx,
            "XGrid" -> xGrid[[idx]]
          |>
        ],
        Partition[Range[Length[xGrid]], UpTo[size]]
      ]
    ],
    patchAssoc
  ]
];

plus2PrepareParallelState[Tfun_, projector_Association] := Module[{},
  LaunchKernels[];
  ParallelEvaluate[
    Clear[
      $Plus2CurrentSource,
      $Plus2CurrentYGrid,
      $Plus2CurrentProjectionMatrix
    ];
    $HistoryLength = 0;
  ];
  $Plus2CurrentSource = Tfun;
  $Plus2CurrentYGrid = projector["YGrid"];
  $Plus2CurrentProjectionMatrix = projector["ProjectionMatrix"];
  DistributeDefinitions[
    $Plus2CurrentSource,
    $Plus2CurrentYGrid,
    $Plus2CurrentProjectionMatrix,
    plus2ProjectCurrentSourceAtX,
    plus2ProjectCurrentChunk
  ];
];

plus2ClearParallelState[] := Module[{},
  ParallelEvaluate[
    Clear[
      $Plus2CurrentSource,
      $Plus2CurrentYGrid,
      $Plus2CurrentProjectionMatrix
    ]
  ];
  Clear[
    $Plus2CurrentSource,
    $Plus2CurrentYGrid,
    $Plus2CurrentProjectionMatrix
  ];
];

plus2AssemblePatchMatrices[
  chunkResults_List,
  patchAssoc_Association
] := Association @ KeyValueMap[
  Function[{patchName, xGrid},
    patchName -> Module[{parts, matrix},
      parts = SortBy[
        Select[chunkResults, Lookup[#, "Patch"] === patchName &],
        First @ Lookup[#, "Indices"] &
      ];
      matrix = Developer`ToPackedArray @ Join @@ Lookup[parts, "Values"];
      If[Length[matrix] =!= Length[xGrid],
        Print["Incorrect row count while reassembling patch ", patchName, "."];
        Print["Expected ", Length[xGrid], "; obtained ", Length[matrix], "."];
        Abort[];
      ];
      matrix
    ]
  ],
  patchAssoc
];

plus2LegacyDataView[
  m_Integer,
  ellList_List,
  patchAssoc_Association,
  patchMatrices_Association,
  nY_Integer
] := Association @ Flatten[
  KeyValueMap[
    Function[{patchName, xGrid},
      MapIndexed[
        Function[{ell, index},
          {+2, m, ell, patchName} -> <|
            "s" -> +2,
            "ell" -> ell,
            "m" -> m,
            "Patch" -> patchName,
            "nY" -> nY,
            "XGrid" -> xGrid,
            "Values" -> patchMatrices[patchName][[All, First[index]]],
            "SourceBuild" -> "windowed-at-perturbation-level"
          |>
        ],
        ellList
      ]
    ],
    patchAssoc
  ],
  1
];

(* ---------------------------------------------------------------------- *)
(* 3. Single-m driver                                                      *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  plus2ProjectionMetadata,
  runPlus2ProjectionForM
];

Options[runPlus2ProjectionForM] = {
  "EllMax" -> 15,
  "OmegaFunction" -> omegaOfM,
  "BasisSpin" -> -2,
  "PhiFactor" -> 2 Pi,
  "ConjugateBasis" -> True,
  "MixedYMap" -> Automatic,
  "SphValuesProvider" -> Automatic,
  "ChunkSize" -> Automatic,
  "KeepLegacyDataView" -> True
};

plus2ProjectionMetadata[
  m_Integer,
  ellMax_Integer,
  projector_Association,
  patchAssoc_Association,
  elapsed_?NumericQ
] := <|
  "Sector" -> "Plus2",
  "SourceSpin" -> +2,
  "BasisSpin" -> projector["BasisSpin"],
  "m" -> m,
  "omega" -> projector["omega"],
  "ellMax" -> ellMax,
  "ellList" -> projector["ellList"],
  "nY" -> Length[projector["yGrid"]],
  "RuleTag" -> projector["RuleTag"],
  "PhiFactor" -> projector["PhiFactor"],
  "ConjugateBasis" -> projector["ConjugateBasis"],
  "PatchSizes" -> Map[Length, patchAssoc],
  "ProjectorIdentityError" -> plus2CheckAngularProjector[projector]["MaxIdentityError"],
  "ElapsedSeconds" -> elapsed,
  "aKerr" -> Quiet @ Check[aKerr, Missing["aKerr"]],
  "MKerr" -> Quiet @ Check[MKerr, Missing["MKerr"]],
  "r0val" -> Quiet @ Check[r0val, Missing["r0val"]],
  "B0val" -> Quiet @ Check[\[ScriptCapitalB]0val, Missing["B0val"]],
  "DateString" -> DateString[{"ISODateTime"}],
  "SystemID" -> $SystemID,
  "VersionNumber" -> $VersionNumber,
  "ObjectVersion" -> 3,
  "SourceBuild" -> "windowed-at-perturbation-level"
|>;

runPlus2ProjectionForM[
  m_Integer,
  patchAssoc_Association,
  Tfun_,
  rule_Association,
  OptionsPattern[]
] := Module[
  {
    ellMax, omegaFunction, omega, projectorOptions, projector,
    chunkSize, tasks, timing, chunkResults, patchMatrices,
    keepLegacy, legacyData, metadata
  },

  ellMax = OptionValue["EllMax"];
  omegaFunction = OptionValue["OmegaFunction"];
  omega = N[omegaFunction[m], MachinePrecision];

  If[! NumericQ[omega],
    Print["omegaOfM[", m, "] is not numerical: ", omega];
    Abort[];
  ];

  projectorOptions = {
    "BasisSpin" -> OptionValue["BasisSpin"],
    "PhiFactor" -> OptionValue["PhiFactor"],
    "ConjugateBasis" -> OptionValue["ConjugateBasis"],
    "MixedYMap" -> OptionValue["MixedYMap"],
    "SphValuesProvider" -> OptionValue["SphValuesProvider"]
  };

  Print["Loading cached spheroidal projector for m = ", m, "..."];
  projector = plus2CachedAngularProjector[
    m,
    ellMax,
    rule,
    omega,
    Sequence @@ projectorOptions
  ];

  Print["  ell range: ", projector["ellList"]];
  Print["  angular nodes: ", Length[projector["yGrid"]]];
  Print["  projector identity error: ", plus2CheckAngularProjector[projector]["MaxIdentityError"]];
  Print["  radial grid sizes: ", Map[Length, patchAssoc]];

  chunkSize = OptionValue["ChunkSize"];
  tasks = plus2MakeChunkTasks[patchAssoc, chunkSize];
  Print["  radial chunk tasks: ", Length[tasks]];

  plus2PrepareParallelState[Tfun, projector];
  CheckAbort[
    timing = AbsoluteTiming[
      chunkResults = ParallelMap[
        plus2ProjectCurrentChunk,
        tasks,
        Method -> "CoarsestGrained"
      ];
    ],
    plus2ClearParallelState[];
    Abort[]
  ];
  plus2ClearParallelState[];

  patchMatrices = plus2AssemblePatchMatrices[chunkResults, patchAssoc];
  keepLegacy = TrueQ[OptionValue["KeepLegacyDataView"]];
  legacyData = If[
    keepLegacy,
    plus2LegacyDataView[
      m,
      projector["ellList"],
      patchAssoc,
      patchMatrices,
      Length[projector["yGrid"]]
    ],
    Missing["NotStored"]
  ];

  metadata = plus2ProjectionMetadata[
    m,
    ellMax,
    projector,
    patchAssoc,
    timing[[1]]
  ];

  <|
    "Metadata" -> metadata,
    "PatchAssoc" -> patchAssoc,
    "ProjectionMatrices" -> patchMatrices,
    "Data" -> legacyData
  |>
];

(* ---------------------------------------------------------------------- *)
(* 4. Atomic per-m storage and resumable all-m driver                       *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  savePlus2ProjectionM,
  loadPlus2ProjectionM,
  plus2ClearModeMemoryIfAvailable,
  runPlus2ProjectionsAllMCached
];

savePlus2ProjectionM[
  root_String,
  ellMax_Integer,
  rule_Association,
  runTag_String,
  m_Integer,
  obj_Association
] := Module[{file},
  file = plus2ProjectionMFile[root, ellMax, rule, runTag, m];
  plus2AtomicDumpSave[file, obj];
  If[! FileExistsQ[file],
    Print["ERROR: projection save failed for m = ", m];
    Print["Target file: ", file];
    Abort[];
  ];
  file
];

loadPlus2ProjectionM[
  root_String,
  ellMax_Integer,
  rule_Association,
  runTag_String,
  m_Integer
] := plus2ReadMXObject[plus2ProjectionMFile[root, ellMax, rule, runTag, m]];

plus2ClearModeMemoryIfAvailable[m_Integer, clearSymbolic_] := Module[{},
  If[DownValues[clearPlus2ModeMemory] =!= {},
    clearPlus2ModeMemory[m, clearSymbolic]
  ];
];

Options[runPlus2ProjectionsAllMCached] = {
  "MList" -> Automatic,
  "CacheRoot" -> Automatic,
  "Overwrite" -> False,
  "ReturnResults" -> False,
  "ClearCompiledSourceAfterSave" -> True,
  "ClearSymbolicModeCacheAfterSave" -> False,
  "OmegaFunction" -> omegaOfM,
  "BasisSpin" -> -2,
  "PhiFactor" -> 2 Pi,
  "ConjugateBasis" -> True,
  "MixedYMap" -> Automatic,
  "SphValuesProvider" -> Automatic,
  "ChunkSize" -> Automatic,
  "KeepLegacyDataView" -> True
};

runPlus2ProjectionsAllMCached[
  ellMax_Integer,
  patchAssoc_Association,
  TfunBuilder_,
  rule_Association,
  OptionsPattern[]
] := Module[
  {
    mList, cacheRoot, overwrite, returnResults, clearCompiled, clearSymbolic,
    perMOptions, runTag, projDir, manifestFile, manifest, savedFiles = {},
    allResults = <||>, m, file, Tfun, obj, savedFile
  },

  mList = Replace[
    OptionValue["MList"],
    Automatic :> Range[-ellMax, ellMax]
  ];
  cacheRoot = Replace[
    OptionValue["CacheRoot"],
    Automatic :> plus2DefaultCacheRoot[]
  ];
  overwrite = TrueQ[OptionValue["Overwrite"]];
  returnResults = TrueQ[OptionValue["ReturnResults"]];
  clearCompiled = TrueQ[OptionValue["ClearCompiledSourceAfterSave"]];
  clearSymbolic = TrueQ[OptionValue["ClearSymbolicModeCacheAfterSave"]];

  perMOptions = {
    "EllMax" -> ellMax,
    "OmegaFunction" -> OptionValue["OmegaFunction"],
    "BasisSpin" -> OptionValue["BasisSpin"],
    "PhiFactor" -> OptionValue["PhiFactor"],
    "ConjugateBasis" -> OptionValue["ConjugateBasis"],
    "MixedYMap" -> OptionValue["MixedYMap"],
    "SphValuesProvider" -> OptionValue["SphValuesProvider"],
    "ChunkSize" -> OptionValue["ChunkSize"],
    "KeepLegacyDataView" -> OptionValue["KeepLegacyDataView"]
  };

  runTag = plus2RunTag[
    patchAssoc,
    OptionValue["BasisSpin"],
    OptionValue["PhiFactor"],
    OptionValue["ConjugateBasis"]
  ];
  projDir = plus2ProjectionDirectory[cacheRoot, ellMax, rule, runTag];
  plus2EnsureDirectory[projDir];
  manifestFile = FileNameJoin[{projDir, "manifest.wl"}];

  manifest = <|
    "Sector" -> "Plus2",
    "SourceSpin" -> +2,
    "BasisSpin" -> OptionValue["BasisSpin"],
    "ellMax" -> ellMax,
    "mList" -> mList,
    "nY" -> Length[rule["yGrid"]],
    "RuleTag" -> plus2RuleTag[rule],
    "ParameterTag" -> plus2ParameterTag[],
    "PatchTag" -> plus2PatchTag[patchAssoc],
    "ConventionTag" -> plus2ConventionTag[
      OptionValue["BasisSpin"],
      OptionValue["PhiFactor"],
      OptionValue["ConjugateBasis"]
    ],
    "RunTag" -> runTag,
    "DateStarted" -> DateString[{"ISODateTime"}],
    "SystemID" -> $SystemID,
    "VersionNumber" -> $VersionNumber,
    "ObjectVersion" -> 3,
    "Files" -> {}
  |>;
  Put[manifest, manifestFile];

  Print["Projection cache directory:"];
  Print[projDir];
  Print["Manifest file:"];
  Print[manifestFile];

  Do[
    file = plus2ProjectionMFile[cacheRoot, ellMax, rule, runTag, m];

    If[FileExistsQ[file] && ! overwrite,
      Print["Skipping cached plus-2 projection for m = ", m];
      Print["  existing file: ", file];
      AppendTo[savedFiles, file];
      If[returnResults,
        allResults[m] = loadPlus2ProjectionM[cacheRoot, ellMax, rule, runTag, m]
      ];
      manifest["Files"] = savedFiles;
      manifest["LastCompletedM"] = m;
      Put[manifest, manifestFile];
      Continue[];
    ];

    Print[""];
    Print["============================================================"];
    Print["Computing plus-2 projection for m = ", m];
    Print["============================================================"];

    Tfun = TfunBuilder[m];
    obj = runPlus2ProjectionForM[
      m,
      patchAssoc,
      Tfun,
      rule,
      Sequence @@ perMOptions
    ];

    savedFile = savePlus2ProjectionM[cacheRoot, ellMax, rule, runTag, m, obj];
    Print["Saved plus-2 projection for m = ", m];
    Print["  ", savedFile];
    AppendTo[savedFiles, savedFile];

    If[returnResults, allResults[m] = obj];
    If[clearCompiled,
      plus2ClearModeMemoryIfAvailable[m, clearSymbolic]
    ];
    Clear[Tfun, obj];
    ClearSystemCache[];

    manifest["Files"] = savedFiles;
    manifest["LastCompletedM"] = m;
    Put[manifest, manifestFile];
    ,
    {m, mList}
  ];

  manifest["DateFinished"] = DateString[{"ISODateTime"}];
  manifest["Files"] = savedFiles;
  Put[manifest, manifestFile];

  If[
    returnResults,
    <|"Manifest" -> manifest, "Results" -> allResults|>,
    manifest
  ]
];

(* ---------------------------------------------------------------------- *)
(* 5. Small preflight check                                                 *)
(* ---------------------------------------------------------------------- *)

ClearAll[preflightPlus2CachedProjector];

Options[preflightPlus2CachedProjector] = Options[plus2BuildAngularProjector];

preflightPlus2CachedProjector[
  m_Integer,
  ellMax_Integer,
  rule_Association,
  omegaFunction_ : omegaOfM,
  OptionsPattern[]
] := Module[{omega, projector},
  omega = N[omegaFunction[m], MachinePrecision];
  If[! NumericQ[omega],
    Print["The supplied frequency is not numerical: ", omega];
    Abort[];
  ];
  projector = plus2CachedAngularProjector[
    m,
    ellMax,
    rule,
    omega,
    "BasisSpin" -> OptionValue["BasisSpin"],
    "PhiFactor" -> OptionValue["PhiFactor"],
    "ConjugateBasis" -> OptionValue["ConjugateBasis"],
    "MixedYMap" -> OptionValue["MixedYMap"],
    "SphValuesProvider" -> OptionValue["SphValuesProvider"]
  ];
  plus2CheckAngularProjector[projector]
];
