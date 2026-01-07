# pymfinder_network_motif_search


This repository provides a specialized implementation of the **pymfinder** framework, adapted for identifying and statistically analyzing motifs.

## Installation

### 1. Prerequisites

First, follow the installation instructions for the core **pymfinder** library at the [pymfinder official repository](https://github.com/stoufferlab/pymfinder).

### 2. Environment Setup

We recommend using a conda environment with python2.x.

```bash
# Create and activate the environment
conda create -n pymfinder_network_motif_search python=2.7
conda activate pymfinder_network_motif_search

# Install dependencies
pip install numpy pandas

```

### 3. Repository Setup

Clone this repository and ensure `main.py` is located in your project's working directory.

---

## Core Functions

### `get_motif`

This function analyzes the **real observed network** and identifies all unique motif instances. It records the node names participating in each motif.

**Usage:**

```python
from main import get_motif

results = get_motif(
    input_path="spatial_network.txt", 
    output_path="real_participation.csv", 
    motifsize=3
)

```

**Output CSV Structure:**
| Column | Description |
| :--- | :--- |
| **MotifID** | The unique numerical identifier for the structural motif (The id is based on the representation of the adjacency matrix of the motif as a binary integer, following the original mfinder). |
| **NodeIDs** | A list of internal integer IDs for the nodes forming the motif. |
| **Types** | A string of the original node names. |

---

### `motif_random`

This function generates a **null model** by creating multiple randomized versions of your network (preserving the degree distribution). It uses parallel processing to quantify how often motifs occur by chance.

**Usage:**

```python
from main import motif_random

motif_random(
    input_path="spatial_network.txt", 
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
| **MotifID** | The unique numerical identifier for the structural motif (The id is based on the representation of the adjacency matrix of the motif as a binary integer, following the original mfinder). |
| **NodeIDs** | A list of internal integer IDs for the nodes forming the motif. |
| **Types** | A string of the original node names. |

---

## Input File Format

The input should be a space-separated or tab-separated `.txt` file representing an **edge list**. Node names should not contain spaces.

**Format:**

```text
<source_node> <target_node>

```


For more details configurations (weighted networks, metropolis algorithms, etc.), please refer to the [pymfinder](https://github.com/stoufferlab/pymfinder/blob/master/documentation/pymfinder_manual.pdf).

---

## Modified Functions
* **Motif compostion analysis**: 
* **Parallel Computing randomization**: Utilizes Python's `multiprocessing.Pool` for high-performance randomization.
