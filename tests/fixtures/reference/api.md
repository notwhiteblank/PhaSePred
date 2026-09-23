# PhaSePred API Documentation

PhaSePred is a comprehensive resource for predicting liquid-liquid phase separation (LLPS)
related proteins. This API provides programmatic read-only access to all prediction results.

- Base URL: `http://predict.phasep.pro/api/v1/`
  (via public IP the same API is reachable at `http://47.88.20.47/phasepred/api/v1/`)
- Data version: `2022-02-11` — 116,806 reviewed proteins across 20 species
- Authentication: none (public read-only)
- Rate limit: 10 requests/second per IP with burst 20; exceeding it returns HTTP `503`
- All responses are JSON. Success: `{"data": ...}`. Error: `{"error": {"code": ..., "message": ...}}`

## Endpoints

| Endpoint | Method | Description |
|---|---|---|
| `/api/v1/` | GET | Service index (version, endpoint list) |
| `/api/v1/meta/` | GET | Dataset statistics, organism list, tool descriptions |
| `/api/v1/protein/<entry>/` | GET | Full prediction results for one UniProt entry |
| `/api/v1/search/` | GET | Keyword search |
| `/api/v1/proteins/batch/` | GET, POST | Batch retrieval (up to 100 entries) |
| `/api/v1/downloads/` | GET | Per-species bulk download links |

### GET `/api/v1/protein/<entry>/`

Full prediction results for a single protein. `<entry>` is a UniProt accession
(case-insensitive, e.g. `O14958`).

Response fields inside `data`:

| Field | Type | Description |
|---|---|---|
| `entry` | string | UniProt accession |
| `entry_name` | string | UniProt entry name (e.g. `CASQ2_HUMAN`) |
| `status` | string | UniProt review status (`reviewed`) |
| `gene_names` | string | Gene names |
| `organism` | string | Organism |
| `sequence` | string | Protein sequence |
| `domain` | object | Pfam domains: `PfamID[]`, `domain[]`, `start[]`, `end[]` |
| `phasepred` | object | PhaSePred meta-predictor scores: `SaPS-8fea`, `PdPS-8fea`, `SaPS-10fea`, `PdPS-10fea` and `*_rnk` ranks |
| `rank` | object | `rank[]` and `method[]` (12 methods) |
| `tools` | object | Per-tool results, see below |

`tools` sub-fields: `catgranule`, `plaac`, `pscore`, `espritz`, `hydropathy`, `deepcoil`,
`seg`, `charge`, `phos`, `deepphase`. Each contains residue-level score arrays
(per-residue values as JSON number arrays), summary scores (`single`/`NLLR`/`FCR`),
region boundaries (`start[]`, `end[]`) and organism rank (`rnk`).

Example:

```bash
curl -s http://predict.phasep.pro/api/v1/protein/O14958/
```

```json
{
  "data": {
    "entry": "O14958",
    "entry_name": "CASQ2_HUMAN",
    "organism": "Homo sapiens (Human)",
    "phasepred": {"SaPS-8fea": 0.43, "PdPS-8fea": 0.60, "...": "..."},
    "tools": {
      "catgranule": {"residue": [0.0, 0.0, 0.27, "..."], "single": 0.12, "start": [1], "end": [45], "rnk": 0.31},
      "...": "..."
    },
    "domain": {"PfamID": ["PF01216"], "domain": ["Calsequestrin"], "start": [2], "end": [381]},
    "sequence": "MALLHSAR..."
  }
}
```

Errors: `404 not_found` if the entry does not exist.

### GET /api/v1/search/

| Parameter | Required | Default | Description |
|---|---|---|---|
| `q` | yes | — | Keyword, at least 2 characters |
| `condition` | no | `Protein_Name` | `Protein_Name` (matches gene/entry names) or `Uniprot_ID` (matches UniProt accessions) |
| `limit` | no | `20` | Page size, 1–100 |
| `offset` | no | `0` | Pagination offset |

Response `data`: `q`, `condition`, `limit`, `offset`, `has_more` (boolean), and
`results[]` with `status`, `entry`, `entry_name`, `organism`, `gene_names`.

```bash
curl -s "http://predict.phasep.pro/api/v1/search/?q=CASQ2&condition=Protein_Name&limit=5"
```

### GET or POST /api/v1/proteins/batch/

Retrieve up to 100 proteins in one request.

- GET: `?entries=O14958,P12345`
- POST: JSON body `{"entries": ["O14958", "P12345"]}` with `Content-Type: application/json`

Entries are case-insensitive and de-duplicated. Response `data`: `count`,
`results[]` (same shape as the single-protein endpoint) and `not_found[]`.

```bash
curl -s "http://predict.phasep.pro/api/v1/proteins/batch/?entries=O14958,P12345"
curl -s -X POST -H "Content-Type: application/json" \
     -d '{"entries": ["O14958", "P12345"]}' \
     http://predict.phasep.pro/api/v1/proteins/batch/
```

Errors: `400 bad_request` if entries are missing, malformed, or exceed 100.

### GET /api/v1/meta/

Dataset statistics: `data_version`, `total_proteins`, `organisms[]` (name + protein
count, 20 species) and `tools` (description of every tool field).

### GET /api/v1/downloads/

Lists per-species bulk download packages. Each item: `file`, `organism`, `url`
(relative to the site root, e.g. `/static/phasepred/database/human_reviewed.zip`),
`size_bytes`. 18 packages are available.

### GET /api/v1/

Service index: name, API version, data version and the endpoint list.

## Error codes

| HTTP | `error.code` | Meaning |
|---|---|---|
| 400 | `bad_request` | Missing/invalid parameters or malformed JSON body |
| 404 | `not_found` | Entry does not exist |
| 405 | `method_not_allowed` | Wrong HTTP method |
| 503 | — | Rate limit exceeded (nginx) |

## Python example

```python
import requests

BASE = "http://predict.phasep.pro/api/v1"

r = requests.get(f"{BASE}/search/", params={"q": "CASQ2", "limit": 5})
entries = [item["entry"] for item in r.json()["data"]["results"]]

r = requests.post(f"{BASE}/proteins/batch/", json={"entries": entries})
for protein in r.json()["data"]["results"]:
    print(protein["entry"], protein["phasepred"]["SaPS-8fea"])
```

## Citation

Chen Z, Hou C, et al. Screening membraneless organelle participants with machine-learning
models that integrate multimodal features. PNAS 2022;119:e2115369119.

## Terms

Free for non-commercial use for academic, government and non-profit institutions.
