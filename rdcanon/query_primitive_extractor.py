from __future__ import annotations

import sys
from pathlib import Path


_AROMATIC_ENCODED_ATOMIC_NUMS = {
    1005: 5,
    1006: 6,
    1007: 7,
    1008: 8,
    1015: 15,
    1016: 16,
}
_AROMATIC_ATOMIC_NUMBERS = {5, 6, 7, 8, 15, 16}
_ORGANIC_SYMBOLS = {
    1: "H",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    15: "P",
    16: "S",
    17: "Cl",
    35: "Br",
    53: "I",
}


def _querymol_parse_smarts():
    try:
        from querymol import parse_smarts
        return parse_smarts
    except Exception:
        pass

    root = Path(__file__).resolve().parents[2]
    smarts_viz = root / "smarts-viz"
    if smarts_viz.exists():
        smarts_viz_str = str(smarts_viz)
        if smarts_viz_str not in sys.path:
            sys.path.insert(0, smarts_viz_str)

    try:
        from querymol import parse_smarts
        return parse_smarts
    except Exception:
        raise ImportError("Unable to import querymol.parse_smarts")


def _charge_token(charge: int) -> str:
    if charge == 0:
        return "+0"
    sign = "+" if charge > 0 else "-"
    mag = abs(charge)
    return f"{sign}{mag}" if mag > 1 else sign


def _atom_token_from_value(value) -> str:
    raw = int(value)
    aromatic = False
    if raw in _AROMATIC_ENCODED_ATOMIC_NUMS:
        raw = _AROMATIC_ENCODED_ATOMIC_NUMS[raw]
        aromatic = True

    symbol = _ORGANIC_SYMBOLS.get(raw)
    if symbol is None:
        return f"#{raw}"
    if aromatic and raw in _AROMATIC_ATOMIC_NUMBERS:
        return symbol.lower()
    return symbol


def _leaf_token(node, effective_negation: bool) -> str:
    description = getattr(node, "description", "")
    value = getattr(node, "value", None)

    if description in ("AtomType", "AtomAtomicNum") and value is not None:
        token = _atom_token_from_value(value)
    elif description == "AtomFormalCharge":
        token = _charge_token(int(value if value is not None else 0))
    elif description == "AtomHCount":
        v = int(value if value is not None else 0)
        token = f"H{v}" if v != 1 else "H"
    elif description == "AtomExplicitDegree":
        v = int(value if value is not None else 0)
        token = f"D{v}" if v != 1 else "D"
    elif description == "AtomTotalDegree":
        v = int(value if value is not None else 0)
        token = f"X{v}" if v != 1 else "X"
    elif description == "AtomIsAromatic":
        token = "a"
    elif description == "AtomIsAliphatic":
        token = "A"
    elif description == "AtomInRing":
        if effective_negation:
            return "R0"
        token = "R"
    elif description == "AtomInNRings":
        if value is None or int(value) < 0:
            token = "R"
        else:
            token = f"R{int(value)}"
    elif description == "AtomMinRingSize":
        token = f"r{int(value)}"
    elif description == "AtomHybridization":
        token = f"^{int(value)}"
    elif description == "AtomTotalValence":
        token = f"v{int(value)}"
    elif description == "AtomNull":
        token = "*"
    elif description == "RecursiveSmarts":
        token = "$(...)"
    else:
        token = description

    if effective_negation and token and not token.startswith("!"):
        token = f"!{token}"
    return token


def _collect_primitives(node, inherited_negation: bool = False) -> list[str]:
    op_obj = getattr(node, "op", None)
    op = getattr(op_obj, "value", op_obj)
    node_neg = bool(getattr(node, "negated", False))
    base_negation = inherited_negation ^ node_neg

    if op == "not":
        out: list[str] = []
        for child in getattr(node, "children", []):
            out.extend(_collect_primitives(child, not base_negation))
        return out

    if op in ("leaf", "recursive") or not getattr(node, "children", []):
        token = _leaf_token(node, base_negation)
        return [token] if token else []

    out: list[str] = []
    for child in getattr(node, "children", []):
        out.extend(_collect_primitives(child, base_negation))
    return out


def extract_primitives_from_atom_query(atom_smarts: str) -> list[str]:
    if not isinstance(atom_smarts, str) or not atom_smarts:
        raise ValueError("atom_smarts must be a non-empty string")

    parse_smarts = _querymol_parse_smarts()
    token = atom_smarts if atom_smarts.startswith("[") else f"[{atom_smarts}]"
    qmol = parse_smarts(token)
    if qmol.GetNumAtoms() == 0:
        raise RuntimeError(f"QueryMol parse returned no atoms for token: {token}")
    node = qmol.GetAtomWithIdx(0).GetQuery()
    primitives = _collect_primitives(node)
    filtered = [
        p for p in primitives
        if isinstance(p, str) and p and "$(...)" not in p
    ]
    if not filtered:
        raise RuntimeError(f"No primitives extracted from token: {token}")
    return filtered
