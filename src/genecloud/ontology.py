"""Minimal Gene Ontology reader (go-basic.obo): names, namespaces and ancestors (is_a + part_of)."""
from __future__ import annotations

import gzip
from functools import lru_cache
from pathlib import Path


class Ontology:
    def __init__(self, terms: dict[str, dict], alt: dict[str, str]):
        self.terms = terms
        self.alt = alt

    @classmethod
    def from_obo(cls, path: str | Path) -> "Ontology":
        opener = gzip.open if str(path).endswith(".gz") else open
        terms, alt, cur = {}, {}, None
        with opener(path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "[Term]":
                    cur = {"parents": []}
                    continue
                if line.startswith("["):
                    cur = None
                    continue
                if cur is None or not line:
                    continue
                key, _, val = line.partition(": ")
                if key == "id":
                    cur["id"] = val
                    terms[val] = cur
                elif key == "name":
                    cur["name"] = val
                elif key == "namespace":
                    cur["namespace"] = val
                elif key == "is_a":
                    cur["parents"].append(val.split(" ")[0])
                elif key == "relationship" and val.startswith("part_of "):
                    cur["parents"].append(val.split(" ")[1])
                elif key == "alt_id":
                    alt[val] = cur.get("id")
                elif key == "is_obsolete" and val == "true":
                    cur["obsolete"] = True
        return cls(terms, alt)

    def resolve(self, go_id: str) -> str | None:
        go_id = self.alt.get(go_id, go_id)
        t = self.terms.get(go_id)
        return None if t is None or t.get("obsolete") else go_id

    def name(self, go_id: str) -> str:
        return self.terms.get(go_id, {}).get("name", go_id)

    def namespace(self, go_id: str) -> str:
        return self.terms.get(go_id, {}).get("namespace", "")

    def ancestors(self, go_id: str) -> frozenset[str]:
        return self._anc(go_id)

    @lru_cache(maxsize=None)
    def _anc(self, go_id: str) -> frozenset[str]:
        t = self.terms.get(go_id)
        if t is None:
            return frozenset()
        out = {go_id}
        for p in t["parents"]:
            out |= self._anc(p)
        return frozenset(out)
