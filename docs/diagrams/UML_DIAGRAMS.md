# System Diagrams (UML)

This document contains UML diagrams describing the AISNP Selection Pipeline architecture, data flow, and system components using PlantUML and Mermaid syntax.

---

## Table of Contents

1. [Use Case Diagram](#1-use-case-diagram)
2. [Component Diagram](#2-component-diagram)
3. [Activity Diagram - Main Pipeline](#3-activity-diagram---main-pipeline)
4. [Activity Diagram - Statistical Selection](#4-activity-diagram---statistical-selection)
5. [Sequence Diagram - ML Training](#5-sequence-diagram---ml-training)
6. [Class Diagram](#6-class-diagram)
7. [Deployment Diagram](#7-deployment-diagram)
8. [State Machine Diagram](#8-state-machine-diagram)
9. [Package Diagram](#9-package-diagram)
10. [Data Flow Diagram (DFD)](#10-data-flow-diagram-dfd)

---

## 1. Use Case Diagram

### PlantUML

```plantuml
@startuml Use Case Diagram
left to right direction
skinparam packageStyle rectangle

actor Researcher as R
actor "Data Scientist" as DS
actor "Bioinformatician" as B

rectangle "AISNP Selection System" {
    usecase "Load VCF Data" as UC1
    usecase "Apply Quality Filters" as UC2
    usecase "Run Statistical Tests" as UC3
    usecase "Train ML Models" as UC4
    usecase "Evaluate Models" as UC5
    usecase "Compare AISNP Panels" as UC6
    usecase "Generate Reports" as UC7
    usecase "Export Results" as UC8
    
    UC1 --> UC2 : <<include>>
    UC2 --> UC3 : <<include>>
    UC3 --> UC4 : <<include>>
    UC4 --> UC5 : <<include>>
    UC5 --> UC7 : <<include>>
    UC6 --> UC5 : <<extend>>
}

R --> UC1
R --> UC7
DS --> UC4
DS --> UC5
DS --> UC6
B --> UC2
B --> UC3
B --> UC8

@enduml
```

### Mermaid

```mermaid
flowchart LR
    subgraph Actors
        R[👤 Researcher]
        DS[👤 Data Scientist]
        B[👤 Bioinformatician]
    end
    
    subgraph "AISNP Selection System"
        UC1((Load VCF Data))
        UC2((Apply Quality Filters))
        UC3((Run Statistical Tests))
        UC4((Train ML Models))
        UC5((Evaluate Models))
        UC6((Compare AISNP Panels))
        UC7((Generate Reports))
        UC8((Export Results))
    end
    
    R --> UC1
    R --> UC7
    DS --> UC4
    DS --> UC5
    DS --> UC6
    B --> UC2
    B --> UC3
    B --> UC8
    
    UC1 -.->|include| UC2
    UC2 -.->|include| UC3
    UC3 -.->|include| UC4
    UC4 -.->|include| UC5
    UC6 -.->|extend| UC5
```

---

## 2. Component Diagram

### PlantUML

```plantuml
@startuml Component Diagram
skinparam componentStyle uml2

package "Data Layer" {
    component [VCF Files] as VCF
    component [PLINK Files] as PLINK
    component [CSV Data] as CSV
    component [BED Files] as BED
    database "1000 Genomes\nRepository" as DB
}

package "Processing Layer" {
    component [HardFilters] as HF
    component [SituationalFilters] as SF
    component [StatisticalTests] as ST
    component [FSTSelection] as FST
    
    interface "IFilter" as IF
    interface "IStatTest" as IST
    
    HF -down- IF
    SF -down- IF
    ST -down- IST
}

package "ML Layer" {
    component [MLTraining] as MLT
    component [ModelEvaluation] as ME
    component [PanelComparison] as PC
    
    interface "IClassifier" as IC
    interface "IEvaluator" as IE
    
    MLT -down- IC
    ME -down- IE
}

package "Utility Layer" {
    component [Config] as CFG
    component [Utils] as UTL
    component [Plotting] as PLT
}

package "External Tools" {
    component [PLINK2] as P2
    component [bcftools] as BCF
    component [Ensembl API] as ENS
    component [scikit-learn] as SKL
}

' Relationships
DB --> VCF : provides
VCF --> HF : input
PLINK --> SF : input
HF --> PLINK : output
SF --> ST : filtered data
ST --> MLT : selected SNPs
MLT --> ME : trained models
PC --> ME : uses

HF ..> P2 : uses
SF ..> P2 : uses
BED ..> ENS : query
MLT ..> SKL : uses

CFG --> HF
CFG --> SF
CFG --> MLT
UTL --> ST
PLT --> ME

@enduml
```

### Mermaid

```mermaid
graph TB
    subgraph "Data Layer"
        VCF[📁 VCF Files]
        PLINK[📁 PLINK Files]
        CSV[📁 CSV Data]
        BED[📁 BED Files]
        DB[(1000 Genomes)]
    end
    
    subgraph "Processing Layer"
        HF[HardFilters]
        SF[SituationalFilters]
        ST[StatisticalTests]
        FST[FSTSelection]
    end
    
    subgraph "ML Layer"
        MLT[MLTraining]
        ME[ModelEvaluation]
        PC[PanelComparison]
    end
    
    subgraph "Utility Layer"
        CFG[Config]
        UTL[Utils]
        PLT[Plotting]
    end
    
    subgraph "External Tools"
        P2[PLINK2]
        BCF[bcftools]
        ENS[Ensembl API]
        SKL[scikit-learn]
    end
    
    DB --> VCF
    VCF --> HF
    HF --> PLINK
    PLINK --> SF
    SF --> ST
    ST --> MLT
    MLT --> ME
    PC --> ME
    
    HF -.-> P2
    SF -.-> P2
    BED -.-> ENS
    MLT -.-> SKL
    
    CFG --> HF
    CFG --> SF
    CFG --> MLT
```

---

## 3. Activity Diagram - Main Pipeline

### PlantUML

```plantuml
@startuml Activity Diagram - Main Pipeline
start

:Load 1000 Genomes VCF;
note right: 84M variants

:Extract EAS Samples;
note right: CHB, JPT, KHV\n(n=306)

partition "Part 1: Statistical Selection" {
    :Apply Hard Filters;
    note right
        - SNP-only
        - Biallelic
        - MAF > 0.16%
        - Call rate > 95%
    end note
    
    :Apply Situational Filters;
    note right
        - HWE test
        - Unique IDs
        - LD pruning
    end note
    
    fork
        :Chi-Squared Test;
    fork again
        :Mutual Information;
    fork again
        :Information Gain;
    fork again
        :KL Divergence;
    end fork
    
    :Compute Intersection;
    note right: 37 Consensus SNPs
    
    :Train ML Models;
    
    :Evaluate with 5-Fold CV;
}

partition "Part 2: Panel Comparison" {
    :Load Known AISNP Panels;
    note right
        - cal_et_al (52)
        - seldin_128 (124)
        - forenseq (55)
        - kidd_55 (53)
        - hsiao_lin_hwa (125)
    end note
    
    :Convert rsID to BED;
    
    :Extract Genotypes;
    
    :Run ML Comparison;
}

:Generate Reports;

:Export Results;
note right
    - CSV files
    - Plots
    - Summary
end note

stop

@enduml
```

### Mermaid

```mermaid
flowchart TD
    A([Start]) --> B[Load 1000 Genomes VCF]
    B --> C[Extract EAS Samples]
    
    subgraph Part1["Part 1: Statistical Selection"]
        C --> D[Apply Hard Filters]
        D --> E[Apply Situational Filters]
        E --> F1[Chi-Squared Test]
        E --> F2[Mutual Information]
        E --> F3[Information Gain]
        E --> F4[KL Divergence]
        F1 --> G[Compute Intersection]
        F2 --> G
        F3 --> G
        F4 --> G
        G --> H[Train ML Models]
        H --> I[Evaluate with 5-Fold CV]
    end
    
    subgraph Part2["Part 2: Panel Comparison"]
        I --> J[Load Known AISNP Panels]
        J --> K[Convert rsID to BED]
        K --> L[Extract Genotypes]
        L --> M[Run ML Comparison]
    end
    
    M --> N[Generate Reports]
    N --> O[Export Results]
    O --> P([End])
```

---

## 4. Activity Diagram - Statistical Selection

### PlantUML

```plantuml
@startuml Activity Diagram - Statistical Selection
|Filtering|
start
:Receive QC-passed SNPs;
note right: N SNPs after QC

|Chi-Squared|
:Create 3x3 Contingency Table;
:Calculate χ² statistic;
:Apply FDR Correction;
if (q-value < 0.05?) then (yes)
    :Add to χ² set;
else (no)
    :Exclude;
endif

|Mutual Information|
:Calculate H(Population);
:Calculate H(Population|Genotype);
:Compute MI = H(P) - H(P|G);
:Rank all SNPs by MI;
:Select Top 500;

|Information Gain|
:Calculate Population Entropy;
:Calculate Conditional Entropy;
:Compute IG = H - Σ P(g)·H(P|g);
:Rank all SNPs by IG;
:Select Top 500;

|KL Divergence|
:Get Genotype Distribution per Pop;
:Calculate Pairwise KL(P||Q);
:Average Across Population Pairs;
:Rank all SNPs by Avg KL;
:Select Top 500;

|Consensus|
:Compute Set Intersection;
note right
    Consensus = χ² ∩ MI ∩ IG ∩ KL
end note

:Calculate Composite Rank;
note right: Average rank across tests

:Output 37 Consensus SNPs;
stop

@enduml
```

### Mermaid

```mermaid
flowchart TD
    subgraph Filtering
        A[Receive QC-passed SNPs]
    end
    
    A --> B & C & D & E
    
    subgraph ChiSquared["Chi-Squared Test"]
        B[Create 3x3 Contingency Table]
        B --> B1[Calculate χ² statistic]
        B1 --> B2[Apply FDR Correction]
        B2 --> B3{q < 0.05?}
        B3 -->|Yes| B4[Add to χ² set]
        B3 -->|No| B5[Exclude]
    end
    
    subgraph MI["Mutual Information"]
        C[Calculate H Population]
        C --> C1[Calculate H Pop given Genotype]
        C1 --> C2[Compute MI]
        C2 --> C3[Rank & Select Top 500]
    end
    
    subgraph IG["Information Gain"]
        D[Calculate Population Entropy]
        D --> D1[Calculate Conditional Entropy]
        D1 --> D2[Compute IG]
        D2 --> D3[Rank & Select Top 500]
    end
    
    subgraph KL["KL Divergence"]
        E[Get Genotype Dist per Pop]
        E --> E1[Calculate Pairwise KL]
        E1 --> E2[Average Across Pairs]
        E2 --> E3[Rank & Select Top 500]
    end
    
    B4 & C3 & D3 & E3 --> F[Compute Set Intersection]
    F --> G[Calculate Composite Rank]
    G --> H[Output 37 Consensus SNPs]
```

---

## 5. Sequence Diagram - ML Training

### PlantUML

```plantuml
@startuml Sequence Diagram - ML Training
skinparam sequenceMessageAlign center

actor User
participant "MLTraining" as MLT
participant "DataLoader" as DL
participant "Preprocessor" as PP
participant "StratifiedKFold" as SKF
participant "Classifier" as CLF
participant "Evaluator" as EVL
database "Results" as RES

User -> MLT: train_and_evaluate(data_path)
activate MLT

MLT -> DL: load_ml_matrix(path)
activate DL
DL --> MLT: DataFrame (306 x 38)
deactivate DL

MLT -> PP: encode_labels(df['pop'])
activate PP
PP --> MLT: encoded_labels
deactivate PP

MLT -> PP: scale_features(X)
activate PP
PP --> MLT: scaled_X
deactivate PP

MLT -> SKF: split(X, y, n_splits=5)
activate SKF
SKF --> MLT: fold_indices
deactivate SKF

loop for each classifier (7 models)
    loop for each fold (5 folds)
        MLT -> CLF: fit(X_train, y_train)
        activate CLF
        CLF --> MLT: trained_model
        deactivate CLF
        
        MLT -> CLF: predict(X_test)
        activate CLF
        CLF --> MLT: y_pred
        deactivate CLF
        
        MLT -> EVL: calculate_metrics(y_test, y_pred)
        activate EVL
        EVL --> MLT: {accuracy, f1, precision, recall}
        deactivate EVL
    end
    
    MLT -> EVL: aggregate_fold_results()
    activate EVL
    EVL --> MLT: {mean, std}
    deactivate EVL
end

MLT -> RES: save_results(results_df)
activate RES
RES --> MLT: saved
deactivate RES

MLT --> User: results_summary
deactivate MLT

@enduml
```

### Mermaid

```mermaid
sequenceDiagram
    actor User
    participant MLT as MLTraining
    participant DL as DataLoader
    participant PP as Preprocessor
    participant SKF as StratifiedKFold
    participant CLF as Classifier
    participant EVL as Evaluator
    participant RES as Results DB

    User->>MLT: train_and_evaluate(data_path)
    activate MLT
    
    MLT->>DL: load_ml_matrix(path)
    activate DL
    DL-->>MLT: DataFrame (306 x 38)
    deactivate DL
    
    MLT->>PP: encode_labels(df['pop'])
    activate PP
    PP-->>MLT: encoded_labels
    deactivate PP
    
    MLT->>PP: scale_features(X)
    activate PP
    PP-->>MLT: scaled_X
    deactivate PP
    
    MLT->>SKF: split(X, y, n_splits=5)
    activate SKF
    SKF-->>MLT: fold_indices
    deactivate SKF
    
    loop For each classifier (7 models)
        loop For each fold (5 folds)
            MLT->>CLF: fit(X_train, y_train)
            activate CLF
            CLF-->>MLT: trained_model
            deactivate CLF
            
            MLT->>CLF: predict(X_test)
            activate CLF
            CLF-->>MLT: y_pred
            deactivate CLF
            
            MLT->>EVL: calculate_metrics(y_test, y_pred)
            activate EVL
            EVL-->>MLT: metrics
            deactivate EVL
        end
        
        MLT->>EVL: aggregate_fold_results()
        activate EVL
        EVL-->>MLT: mean ± std
        deactivate EVL
    end
    
    MLT->>RES: save_results(results_df)
    activate RES
    RES-->>MLT: saved
    deactivate RES
    
    MLT-->>User: results_summary
    deactivate MLT
```

---

## 6. Class Diagram

### PlantUML

```plantuml
@startuml Class Diagram
skinparam classAttributeIconSize 0

package "config" {
    class PathConfig {
        +PROJECT_ROOT: Path
        +DATA_DIR: Path
        +OUTPUT_DIR: Path
        +GRAPHS_DIR: Path
        +get_path(name: str): Path
    }
    
    class FilterParams {
        +MAF_THRESHOLD: float = 0.0016
        +CALL_RATE: float = 0.05
        +HWE_PVALUE: float = 1e-6
        +LD_R2: float = 0.1
        +LD_WINDOW: int = 1000
    }
}

package "filters" {
    abstract class BaseFilter {
        +{abstract} apply(data: DataFrame): DataFrame
        +validate_input(data: DataFrame): bool
    }
    
    class HardFilter {
        -params: FilterParams
        +apply(data: DataFrame): DataFrame
        +filter_snps_only(): DataFrame
        +filter_biallelic(): DataFrame
        +filter_maf(): DataFrame
        +filter_call_rate(): DataFrame
    }
    
    class SituationalFilter {
        -params: FilterParams
        +apply(data: DataFrame): DataFrame
        +filter_hwe(): DataFrame
        +deduplicate_ids(): DataFrame
        +prune_ld(): DataFrame
    }
    
    BaseFilter <|-- HardFilter
    BaseFilter <|-- SituationalFilter
    HardFilter --> FilterParams
    SituationalFilter --> FilterParams
}

package "statistical" {
    abstract class StatisticalTest {
        +{abstract} compute(genotypes: ndarray, labels: ndarray): float
        +{abstract} select_top(scores: dict, n: int): list
    }
    
    class ChiSquaredTest {
        -alpha: float
        -correction: str
        +compute(genotypes, labels): float
        +create_contingency_table(): ndarray
        +apply_fdr_correction(): ndarray
    }
    
    class MutualInformation {
        +compute(genotypes, labels): float
        +calculate_entropy(X): float
        +calculate_conditional_entropy(X, Y): float
    }
    
    class InformationGain {
        +compute(genotypes, labels): float
    }
    
    class KLDivergence {
        +compute(genotypes, labels): float
        +pairwise_kl(P, Q): float
    }
    
    class ConsensusSelector {
        -tests: list[StatisticalTest]
        +run_all_tests(data): dict
        +compute_intersection(): set
        +rank_by_composite(): list
    }
    
    StatisticalTest <|-- ChiSquaredTest
    StatisticalTest <|-- MutualInformation
    StatisticalTest <|-- InformationGain
    StatisticalTest <|-- KLDivergence
    ConsensusSelector o-- StatisticalTest
}

package "ml" {
    abstract class BaseClassifier {
        +{abstract} fit(X, y): void
        +{abstract} predict(X): ndarray
        +{abstract} predict_proba(X): ndarray
    }
    
    class ClassifierFactory {
        +{static} get_classifiers(): dict
        +{static} create(name: str): BaseClassifier
    }
    
    class CrossValidator {
        -n_folds: int
        -stratified: bool
        +validate(clf, X, y): dict
        +aggregate_results(): DataFrame
    }
    
    class ModelEvaluator {
        +accuracy(y_true, y_pred): float
        +f1_score(y_true, y_pred): float
        +confusion_matrix(y_true, y_pred): ndarray
        +classification_report(): str
    }
    
    class MLPipeline {
        -classifiers: dict
        -cv: CrossValidator
        -evaluator: ModelEvaluator
        +run(data: DataFrame): DataFrame
        +compare_sources(sources: list): DataFrame
    }
    
    MLPipeline --> ClassifierFactory
    MLPipeline --> CrossValidator
    MLPipeline --> ModelEvaluator
    CrossValidator --> BaseClassifier
}

package "part2" {
    class RSIDConverter {
        -cache: dict
        +fetch_coordinates(rsid: str): tuple
        +batch_convert(rsids: list): DataFrame
        +save_bed(coords: DataFrame, path: str): void
    }
    
    class GenotypeExtractor {
        -plink_path: str
        +extract_snps(bed: str, pfile: str): DataFrame
        +to_ml_matrix(): DataFrame
    }
    
    class PanelComparator {
        -sources: dict
        -pipeline: MLPipeline
        +load_all_panels(): dict
        +compare_all(): DataFrame
        +generate_report(): str
    }
    
    PanelComparator --> MLPipeline
    PanelComparator --> GenotypeExtractor
}

@enduml
```

### Mermaid

```mermaid
classDiagram
    class PathConfig {
        +Path PROJECT_ROOT
        +Path DATA_DIR
        +Path OUTPUT_DIR
        +get_path(name) Path
    }
    
    class FilterParams {
        +float MAF_THRESHOLD
        +float CALL_RATE
        +float HWE_PVALUE
        +float LD_R2
    }
    
    class BaseFilter {
        <<abstract>>
        +apply(data)* DataFrame
        +validate_input(data) bool
    }
    
    class HardFilter {
        -FilterParams params
        +apply(data) DataFrame
        +filter_snps_only() DataFrame
        +filter_maf() DataFrame
    }
    
    class SituationalFilter {
        -FilterParams params
        +apply(data) DataFrame
        +filter_hwe() DataFrame
        +prune_ld() DataFrame
    }
    
    class StatisticalTest {
        <<abstract>>
        +compute(genotypes, labels)* float
        +select_top(scores, n)* list
    }
    
    class ChiSquaredTest {
        -float alpha
        +compute(genotypes, labels) float
        +apply_fdr_correction() ndarray
    }
    
    class MutualInformation {
        +compute(genotypes, labels) float
        +calculate_entropy(X) float
    }
    
    class ConsensusSelector {
        -list~StatisticalTest~ tests
        +run_all_tests(data) dict
        +compute_intersection() set
    }
    
    class MLPipeline {
        -dict classifiers
        -CrossValidator cv
        +run(data) DataFrame
        +compare_sources(sources) DataFrame
    }
    
    class CrossValidator {
        -int n_folds
        +validate(clf, X, y) dict
    }
    
    class ModelEvaluator {
        +accuracy(y_true, y_pred) float
        +f1_score(y_true, y_pred) float
        +confusion_matrix(y_true, y_pred) ndarray
    }
    
    BaseFilter <|-- HardFilter
    BaseFilter <|-- SituationalFilter
    HardFilter --> FilterParams
    StatisticalTest <|-- ChiSquaredTest
    StatisticalTest <|-- MutualInformation
    ConsensusSelector o-- StatisticalTest
    MLPipeline --> CrossValidator
    MLPipeline --> ModelEvaluator
```

---

## 7. Deployment Diagram

### PlantUML

```plantuml
@startuml Deployment Diagram
skinparam nodeStyle rectangle

node "Development Machine" as DM {
    node "Python Environment" as PY {
        component [Jupyter Notebooks] as JN
        component [Python Scripts] as PS
        
        artifact "pandas" as PD
        artifact "scikit-learn" as SKL
        artifact "numpy" as NP
        artifact "matplotlib" as MPL
        artifact "xgboost" as XGB
    }
    
    node "External Tools" as ET {
        component [PLINK2] as P2
        component [bcftools] as BCF
    }
    
    database "Local Storage" as LS {
        folder "1000genomes/" as G1K
        folder "output/" as OUT
        folder "graphs/" as GRP
        folder "docs/" as DOC
    }
}

cloud "External Services" as ES {
    component [Ensembl REST API] as ENS
    database "1000 Genomes FTP" as FTP
}

JN --> PS : imports
PS --> PD
PS --> SKL
PS --> NP
PS --> MPL
PS --> XGB

PS --> P2 : subprocess
PS --> BCF : subprocess

PS --> LS : read/write
PS --> ENS : HTTP requests

FTP --> G1K : download VCF

@enduml
```

### Mermaid

```mermaid
flowchart TB
    subgraph DM["Development Machine"]
        subgraph PY["Python Environment"]
            JN[Jupyter Notebooks]
            PS[Python Scripts]
            subgraph Libs["Libraries"]
                PD[pandas]
                SKL[scikit-learn]
                NP[numpy]
                XGB[xgboost]
            end
        end
        
        subgraph ET["External Tools"]
            P2[PLINK2]
            BCF[bcftools]
        end
        
        subgraph LS["Local Storage"]
            G1K[(1000genomes/)]
            OUT[(output/)]
            GRP[(graphs/)]
        end
    end
    
    subgraph ES["External Services"]
        ENS[Ensembl REST API]
        FTP[(1000 Genomes FTP)]
    end
    
    JN --> PS
    PS --> Libs
    PS --> P2
    PS --> BCF
    PS --> LS
    PS --> ENS
    FTP --> G1K
```

---

## 8. State Machine Diagram

### PlantUML

```plantuml
@startuml State Machine - SNP Processing
skinparam stateBackgroundColor LightBlue

[*] --> RawVariant : Load VCF

state "Quality Control" as QC {
    state RawVariant
    state SNPFiltered
    state BiallelicFiltered
    state MAFFiltered
    state CallRateFiltered
    
    RawVariant --> SNPFiltered : Apply SNP-only filter
    SNPFiltered --> BiallelicFiltered : Apply biallelic filter
    BiallelicFiltered --> MAFFiltered : Apply MAF filter
    MAFFiltered --> CallRateFiltered : Apply call rate filter
}

state "Population Filtering" as PF {
    state HWEFiltered
    state Deduplicated
    state LDPruned
    
    CallRateFiltered --> HWEFiltered : Apply HWE test
    HWEFiltered --> Deduplicated : Remove duplicates
    Deduplicated --> LDPruned : LD pruning
}

state "Statistical Selection" as SS {
    state Tested
    state Ranked
    state ConsensusSelected
    
    LDPruned --> Tested : Run 4 tests
    Tested --> Ranked : Calculate composite score
    Ranked --> ConsensusSelected : Select intersection
}

state "ML Validation" as ML {
    state Trained
    state Evaluated
    state Validated
    
    ConsensusSelected --> Trained : Train classifiers
    Trained --> Evaluated : Cross-validate
    Evaluated --> Validated : Final evaluation
}

Validated --> [*] : Export Results

@enduml
```

### Mermaid

```mermaid
stateDiagram-v2
    [*] --> RawVariant: Load VCF
    
    state "Quality Control" as QC {
        RawVariant --> SNPFiltered: SNP-only filter
        SNPFiltered --> BiallelicFiltered: Biallelic filter
        BiallelicFiltered --> MAFFiltered: MAF filter
        MAFFiltered --> CallRateFiltered: Call rate filter
    }
    
    state "Population Filtering" as PF {
        CallRateFiltered --> HWEFiltered: HWE test
        HWEFiltered --> Deduplicated: Remove duplicates
        Deduplicated --> LDPruned: LD pruning
    }
    
    state "Statistical Selection" as SS {
        LDPruned --> Tested: Run 4 tests
        Tested --> Ranked: Composite score
        Ranked --> ConsensusSelected: Intersection
    }
    
    state "ML Validation" as ML {
        ConsensusSelected --> Trained: Train classifiers
        Trained --> Evaluated: Cross-validate
        Evaluated --> Validated: Final eval
    }
    
    Validated --> [*]: Export Results
```

---

## 9. Package Diagram

### PlantUML

```plantuml
@startuml Package Diagram
skinparam packageStyle folder

package "AISNP Selection Pipeline" {
    
    package "scripts" <<Folder>> {
        package "config" {
            [PathConfig]
            [FilterParams]
            [MLParams]
        }
        
        package "filters" {
            [hard_filters.py]
            [situational_filters.py]
        }
        
        package "statistical" {
            [statistical_tests.py]
            [fst_selection.py]
        }
        
        package "ml" {
            [ml_training.py]
            [ml_evaluation.py]
        }
        
        package "utils" {
            [utils.py]
            [plotting.py]
        }
        
        package "part2" {
            [rsid_utils.py]
            [bed_to_matrix.py]
            [ml_comparison.py]
        }
        
        package "notebooks" {
            [01_hard_filtering.ipynb]
            [02_situational_filtering.ipynb]
            [02b_statistical_snp_selection.ipynb]
            [03_fst_and_pca.ipynb]
            [04c_ml_consensus_snps.ipynb]
            
            package "part2" as nb_p2 {
                [06_rsid_to_bed.ipynb]
                [07_bed_to_ml_matrix.ipynb]
                [08_known_aisnps_ml.ipynb]
            }
        }
    }
    
    package "data" <<Database>> {
        [1000genomes/]
        [known_aisnps/]
    }
    
    package "output" <<Folder>> {
        [ml_models/]
        [part2/]
    }
    
    package "docs" <<Folder>> {
        [PIPELINE.md]
        [RESULTS.md]
        [diagrams/]
        [slides/]
    }
}

' Dependencies
filters ..> config : uses
statistical ..> config : uses
ml ..> config : uses
part2 ..> ml : uses
notebooks ..> filters : imports
notebooks ..> statistical : imports
notebooks ..> ml : imports

@enduml
```

### Mermaid

```mermaid
flowchart TB
    subgraph scripts["📁 scripts"]
        subgraph config["config"]
            PC[PathConfig]
            FP[FilterParams]
        end
        
        subgraph filters["filters"]
            HF[hard_filters.py]
            SF[situational_filters.py]
        end
        
        subgraph statistical["statistical"]
            ST[statistical_tests.py]
            FS[fst_selection.py]
        end
        
        subgraph ml["ml"]
            MT[ml_training.py]
            ME[ml_evaluation.py]
        end
        
        subgraph part2["part2"]
            RU[rsid_utils.py]
            BM[bed_to_matrix.py]
            MC[ml_comparison.py]
        end
        
        subgraph notebooks["notebooks"]
            N01[01_hard_filtering]
            N02[02_situational_filtering]
            N02b[02b_statistical_snp]
            N04[04c_ml_consensus]
            subgraph nb_p2["part2"]
                N06[06_rsid_to_bed]
                N08[08_known_aisnps_ml]
            end
        end
    end
    
    subgraph data["📁 data"]
        G1K[(1000genomes/)]
        KA[(known_aisnps/)]
    end
    
    subgraph output["📁 output"]
        ML[(ml_models/)]
        P2[(part2/)]
    end
    
    filters -.-> config
    statistical -.-> config
    ml -.-> config
    part2 -.-> ml
```

---

## 10. Data Flow Diagram (DFD)

### PlantUML - Level 0 (Context)

```plantuml
@startuml DFD Level 0
skinparam rectangle {
    BackgroundColor LightYellow
}

actor "Researcher" as R

rectangle "AISNP Selection\nSystem" as SYS

database "1000 Genomes" as G1K
database "Known AISNP\nPanels" as KP
database "Results\nRepository" as RES

R --> SYS : Configuration\nParameters
G1K --> SYS : VCF Data\n(84M variants)
KP --> SYS : rsID Lists\n(6 sources)

SYS --> R : Performance\nReports
SYS --> RES : Selected SNPs\nTrained Models

@enduml
```

### PlantUML - Level 1 (Main Processes)

```plantuml
@startuml DFD Level 1
skinparam rectangle {
    BackgroundColor LightBlue
}

database "1000 Genomes\nVCF" as VCF
database "Known AISNP\nPanels" as KP

rectangle "1.0\nQuality\nFiltering" as P1
rectangle "2.0\nStatistical\nSelection" as P2
rectangle "3.0\nML\nTraining" as P3
rectangle "4.0\nPanel\nComparison" as P4
rectangle "5.0\nReport\nGeneration" as P5

database "Filtered\nPLINK Files" as D1
database "Consensus\nSNPs" as D2
database "ML\nResults" as D3
database "Comparison\nResults" as D4

VCF --> P1 : Raw variants
P1 --> D1 : QC-passed variants
D1 --> P2 : Filtered data
P2 --> D2 : 37 SNPs
D2 --> P3 : Selected SNPs
P3 --> D3 : Model metrics

KP --> P4 : rsID lists
D1 --> P4 : Genotype source
P4 --> D4 : Panel results

D3 --> P5 : Part 1 results
D4 --> P5 : Part 2 results
P5 --> Reports : Final reports

@enduml
```

### Mermaid - DFD Level 1

```mermaid
flowchart LR
    VCF[(1000 Genomes VCF)]
    KP[(Known AISNP Panels)]
    
    P1[/"1.0 Quality Filtering"/]
    P2[/"2.0 Statistical Selection"/]
    P3[/"3.0 ML Training"/]
    P4[/"4.0 Panel Comparison"/]
    P5[/"5.0 Report Generation"/]
    
    D1[(Filtered PLINK Files)]
    D2[(Consensus SNPs)]
    D3[(ML Results)]
    D4[(Comparison Results)]
    Reports[(Final Reports)]
    
    VCF --> P1
    P1 --> D1
    D1 --> P2
    P2 --> D2
    D2 --> P3
    P3 --> D3
    
    KP --> P4
    D1 --> P4
    P4 --> D4
    
    D3 --> P5
    D4 --> P5
    P5 --> Reports
```

---

## Diagram Rendering Tools

### PlantUML Rendering
- **Online**: [PlantUML Web Server](http://www.plantuml.com/plantuml)
- **VS Code**: PlantUML extension
- **IntelliJ**: PlantUML integration plugin
- **Command Line**: `java -jar plantuml.jar diagram.puml`

### Mermaid Rendering
- **GitHub**: Native support in Markdown
- **VS Code**: Markdown Preview Mermaid Support extension
- **Online**: [Mermaid Live Editor](https://mermaid.live)
- **Documentation**: [Mermaid Docs](https://mermaid.js.org/)

---

## Export Formats

| Format | PlantUML | Mermaid |
|--------|----------|---------|
| PNG | ✅ | ✅ |
| SVG | ✅ | ✅ |
| PDF | ✅ | ✅ |
| ASCII | ✅ | ❌ |
| LaTeX | ✅ | ❌ |

---

*UML Diagrams for AISNP Selection Pipeline*  
*Version 2.0 - January 2026*
