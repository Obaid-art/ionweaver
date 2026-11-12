<p align="center"><img src="web/mark.svg" width="64" alt=""></p>
<h1 align="center">IonWeaver</h1>
<p align="center"><em>Make new ions that keep their family.</em></p>

---

An edit can produce a perfectly valid molecule and still destroy the thing that
made the ion useful: the charged group that defines its family. IonWeaver finds
the charge-bearing core the members of a family share, protects every atom and
every bond inside it, and changes only the neutral body around it. Each edit is
decided on the parent graph before anything is built, and every product is
certified five ways before it is offered. Nothing is repaired, retried or
filtered afterwards.

This repository is the web application: a static page (`web/`) and a small
engine server (`server/`) that runs the submitted manuscript's engine
(`engine/iw/`), module for module.

## Use it

**In the browser.** Open the page. Three ways to begin:

1. **A published family.** Pick one of 38 families of curated ions, choose an
   ion, say how many new ions you want, and optionally set **rules** beforehand:
   which elements may appear outside the core, no double or triple bonds, no
   rings, no heteroatoms, at most so many atoms outside the core, and bounds on
   the whole ion (molecular weight, H-bond acceptors and donors, log P, polar
   surface area, rotatable bonds, rings). The engine keeps examining until that
   many ions obey every rule; the core itself is never touched. The same rules
   apply across families on the Across-families page, so a whole set can be
   "no halogens, saturated, acyclic, one donor at most" in one go. Filter and sort what came back, ask for more (the next ions
   further down the study's own examination order), look at each one with its
   protected core in blue and the edit in amber, and download the shown set as
   CSV, SDF or SMILES. The workbench also shows the ion's parts: the core on
   its own and each side group, drawn as the library holds them.
2. **Your own ion.** Paste one SMILES string. IonWeaver names its family, shows
   the core it will protect, and offers the other published cores it also
   matches.
3. **Your own family.** Upload a sheet with two columns, `SMILES` and `family`
   (two ions per family are enough). The sheet is curated exactly as the
   published inventory was, the shared core of each family is discovered with
   the study's own layered search, and every member can then be edited against
   the published library plus the sheet's own components.

A family page also generates directly: say how many new ions you want from
the family and a spread of its ions serve as parents, no need to pick one.
From any starting ion you can also **edit again and again** (sequential
editing, up to twelve rounds, under the study's three selection rules), and
from the **Across families** page you can generate a quota of distinct ions per
family over every published family, or over the families of your sheet, as one
CSV with structures, provenance, the five checks and the descriptors.

**The studio** builds an ion by hand, one certified edit at a time. Point at
an atom or a bond; the panel shows every move the engine offers there and
every one it does not: attach any of the library's components (Graft), append
an element at a single, double or triple bond (Insert), change an element
(Substitute), remove an atom (Delete), exchange a whole side group (Swap),
cut a tail (Prune), change a bond order (Rebond), splice a different tail
(Cross-over). A move the engine offers is built and certified and becomes the
next ion; a move it does not offer is explained from the rule that stops it:
a protected core atom, a ring bond, a valence the ledger refuses, an atom or
bond environment the inventory never displays. "Let the engine choose" takes
the next certified move for you, at random, close to the start or far from
it. The trail can be undone step by step, downloaded, or handed to the
workbench to generate N ions from where you stopped.

**Locally.**

```bash
pip install -r requirements.txt
uvicorn app:app --app-dir server --port 8770
```

then open <http://localhost:8770>. Everything runs on your own computer.

**Tests.** `pytest -q` checks that the engine files are the notebook's, that
the library rebuilt at start-up reproduces the archive, that per-operator
capacities match the published counts, and that products come out in the
published examination order.

## How it is built

```
engine/iw/        the 18 engine modules of the submitted notebook, byte for byte
engine/MANIFEST.json   their SHA-256 values; the server re-hashes them at start-up
data/archive/     the archive tables the server rebuilds the published state from
server/           the FastAPI application and the thin layer that asks the engine
web/              the page: index.html, style.css, app.js, config.js
tests/            the reproduction suite and its fixture
tools/vendor_notebook.py   re-vendors the engine from a notebook
```

**The chemistry is the published chemistry.** `engine/iw/` is extracted from the
notebook by `tools/vendor_notebook.py` with no change; `MANIFEST.json` records
the hashes and the Method page shows whether they still match. The server
layer (`server/ionweaver_web/`) never decides chemistry: definitions come from
`operators.definitions_for`, products from `build.realise`, certificates from
`build.certify`, the candidate order from `single.order_for`, core discovery
from `discovery.discover`, curation from `curation.curate`. The only loop
written here is the sequential walk, so that its depth can be chosen; every
step inside it is the engine's.

**The published state is rebuilt, not trusted.** At start-up the 2,431 source
entries are re-curated, the environment registries re-learned and the editable
set and fragment library rebuilt; each is compared with the archive's own table,
and a server that does not reproduce the run refuses to start. The
maximum-common-substructure search is the one stage read from the archive
rather than re-run (it takes 35 s and is deterministic).

**Distinctness is a guarantee.** A delivered product is never the parent,
never another product of the same request, and in a batch never a product of
another family's request. Candidates are examined in the engine's keyed order
under the study's own seed (fixed inside the server; there is no seed to set),
so a reference ion's first products are the published ones and "more" simply
walks further down the same order. Batch parents are chosen the same way every
time: a family's editable ions sorted by body size, taken at evenly spaced
positions from smallest to largest.

**The fragment library is on show.** The 255 side groups the curated ions
display, each with its attachment bond, occurrences, families and which
operators may use it, so a user can see exactly what Graft, Swap and Cross-over
draw on.

## Deploy

The page is static and the engine needs Python with RDKit, so the two live in
two places and `web/config.js` says where the engine is.

* **GitHub Pages** serves `web/` (the `deploy` workflow publishes it on every
  push to `main`).
* **A Hugging Face Space** (Docker SDK) runs the engine from the `Dockerfile`.
  Set the repository variable `HF_SPACE` (for example `obaid-art/ionweaver`)
  and the secret `HF_TOKEN`, and the same workflow pushes the repository to the
  Space. Put the Space's URL in `web/config.js`. The free CPU tier is enough:
  a request for 60 ions takes well under a second, a batch over every family at
  100 ions each under half a minute.

The Space also serves the page itself at its root, so either address works.

## What it does not claim

Each output is one cation or one anion, not an ionic liquid: pairing with a
counter-ion, property prediction and synthesis assessment all come afterwards.
An ion absent from the reference inventory is not thereby unknown in the wider
literature. Stereochemistry is not preserved or evaluated. Core discovery is
automatic, but the initial family grouping is not. The descriptors are RDKit's
atom-contribution estimates on a charged ion and are offered as filters, not as
predictions.

## Licence

MIT.
