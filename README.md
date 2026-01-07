# pymfinder_network_motif_search

This repository provides an optimized implementation of the **pymfinder** framework, specifically adapted for identifying network motifs and statistically analyzing its composition.

## Key Features

* **Motif Composition Analysis**: Unlike original pymfinder that only provide total counts of motifs, this implementation output the specific node names for each motif instance.
* **Parallel Computing Randomization**: Employs `multiprocessing.Pool` for high-performance randomization, significantly reducing the computational time required for large networks across multiple CPU cores.


## Installation

### 1. Prerequisites

Follow the installation instructions for the core **pymfinder** library at the [official pymfinder repository](https://github.com/stoufferlab/pymfinder).

### 2. Environment Setup

The underlying `mfinder` C-extension is best supported in a Python 2.x environment. We recommend using **Conda**:

```bash
# Create and activate the environment
conda create -n pymfinder_motif python=2.7
conda activate pymfinder_motif

# Install required dependencies
pip install numpy pandas

```

### 3. Repository Setup

Clone this repository and ensure that `main.py` is located within your project's working directory.


## Core Functions

### *get_motif*

This function performs analysis based on the real observed network. It identifies all unique motif instances and records the identities of the nodes involved.

**Usage:**

```python
from main import get_motif

results = get_motif(
    input_path="network.txt", 
    output_path="real_participation.csv", 
    motifsize=3
)

```

**Output CSV Structure:**
| Column | Description |
| :--- | :--- |
| **MotifID** | The unique numerical identifier for the structural motif. The ID is based on the binary integer representation of the motif's adjacency matrix, following the original `mfinder` convention. |
| **NodeIDs** | A list of internal integer IDs assigned to the nodes forming the motif. |
| **Types** | Original node names joined by '-' participating in the motif. |


### *motif_random*

This function generates a null model by creating multiple randomized versions of the network while preserving the degree distribution. It utilizes parallel processing to quantify how frequently specific motifs occur by chance.

**Usage:**

```python
from main import motif_random

motif_random(
    input_path="network.txt", 
    output_path="random_participation.csv", 
    n_randomizations=100, 
    num_cores=40, 
    motifsize=3
)

```

**Output CSV Structure:**
| Column | Description |
| :--- | :--- |
| **Iteration** | The index of the randomization run. |
| **MotifID** | The unique numerical identifier for the structural motif identified in the randomized network. |
| **NodeIDs** | A list of internal integer IDs for the nodes forming the motif in the randomized graph. |
| **Types** | Original node names joined by '-' participating in the motif. |


## Input File Format

The input should be a space-separated or tab-separated `.txt` file representing an **edge list**. Node names must not contain spaces.

**Format:**

```text
<source_node> <target_node>

```

For more detailed configurations (weighted networks, Metropolis algorithms, etc.), please refer to the [pymfinder documentation](https://github.com/stoufferlab/pymfinder/blob/master/documentation/pymfinder_manual.pdf).

---
