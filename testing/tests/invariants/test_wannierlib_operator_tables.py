"""The three operator tables, the name union, and the schema must agree.

WANNIERLIB_OP_NAMES is the list the reader walks: for each name in it, it looks for an
element of that name under <operators> and reads the switches the exposure tables say
that name offers. Nothing computes it from the tables, because a Fortran PARAMETER here
cannot be built from other PARAMETERs, so the union is written out by hand -- and a hand
is what this guards.

The failure it exists for is silent. Add a row to WANNIERLIB_OPR and forget the union,
and the reader never looks for that element: the input validates, the run starts, rc is
0, and the operator is simply not written. Nothing says so. The same holds the other way
between the tables and the schema: a switch a table offers but the schema has no
attribute for reads its default and stays off for ever.

Source text only -- no binary, no run -- so it costs nothing and holds on any machine.
"""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
OPTABLE = ROOT / "src/libraries/fleurinput/types_melem_optable.f90"
READER = ROOT / "src/libraries/fleurinput/types_wannierlib.F90"
SCHEMA = ROOT / "src/fleur/io/xml/FleurInputSchema.xsd"

#> Which exposure table backs which switch. This mapping IS the contract between the
#> Fortran tables and the XML attributes, and it is the one thing stated in both places.
TABLE_OF_SWITCH = {"ks": "WANNIERLIB_OPK", "rs": "WANNIERLIB_OPR", "ki": "WANNIERLIB_INTERP"}


def _strip_comments(text):
    """Drop whole-line Fortran comments, which carry brackets and quotes of their own."""
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("!"))


def _array_body(text, name):
    """The text of a PARAMETER array initialiser, and the length it declares."""
    m = re.search(name + r"\((\d+)\)\s*=\s*\[", text)
    assert m, f"{name} is not declared as a sized PARAMETER array in {OPTABLE.name}"
    start = m.end()
    end = text.index("]", start)
    return text[start:end], int(m.group(1))


@pytest.fixture(scope="module")
def tables():
    text = _strip_comments(OPTABLE.read_text())
    out = {}
    for name in ("WANNIERLIB_INTERP", "WANNIERLIB_OPR", "WANNIERLIB_OPK"):
        body, declared = _array_body(text, name)
        names = re.findall(r"t_melem_exposed\(\s*'([^']*)'", body)
        assert len(names) == declared, (
            f"{name} declares {declared} entries and holds {len(names)}; the declared "
            f"length is what every loop over it trusts")
        out[name] = names
    body, declared = _array_body(text, "WANNIERLIB_OP_NAMES")
    union = re.findall(r"'([^']*)'", body)
    assert len(union) == declared, (
        f"WANNIERLIB_OP_NAMES declares {declared} names and holds {len(union)}")
    out["union"] = union
    return out


@pytest.fixture(scope="module")
def schema_operators():
    """{element name: {switches its type offers}} for the children of <operators>."""
    text = SCHEMA.read_text()
    m = re.search(r'<xsd:complexType name="WannierlibOperatorsType">(.*?)</xsd:complexType>',
                  text, re.S)
    assert m, "WannierlibOperatorsType is gone from the schema"
    children = re.findall(r'<xsd:element[^>]*\bname="(\w+)"[^>]*\btype="(\w+)"', m.group(1))
    assert children, "WannierlibOperatorsType declares no operator elements"
    out = {}
    for ename, tname in children:
        tm = re.search(r'<xsd:complexType name="%s">(.*?)</xsd:complexType>' % tname, text, re.S)
        assert tm, f"<{ename}> has type {tname}, which the schema does not define"
        attrs = set(re.findall(r'<xsd:attribute[^>]*\bname="(ks|rs|ki)"', tm.group(1)))
        assert attrs, f"{tname} offers no ks/rs/ki attribute, so <{ename}> can ask for nothing"
        out[ename] = attrs
    return out


@pytest.mark.wannierlib
def test_the_union_is_exactly_the_three_tables(tables):
    """A name in a table and not in the union is an operator the reader never looks for."""
    union = set(tables["union"])
    assert len(union) == len(tables["union"]), (
        f"WANNIERLIB_OP_NAMES repeats a name: {sorted(n for n in union if tables['union'].count(n) > 1)}")
    tabled = set(tables["WANNIERLIB_INTERP"]) | set(tables["WANNIERLIB_OPR"]) | set(tables["WANNIERLIB_OPK"])
    assert union == tabled, (
        f"the union and the tables disagree. In a table and unreachable from the input: "
        f"{sorted(tabled - union)}. In the union and in no table: {sorted(union - tabled)}")


@pytest.mark.wannierlib
def test_the_schema_offers_exactly_the_tabled_switches(tables, schema_operators):
    """Every legal (operator, basis) pair has an attribute, and no illegal pair has one."""
    assert set(schema_operators) == set(tables["union"]), (
        f"the schema and WANNIERLIB_OP_NAMES disagree on the operator list. Only in the "
        f"schema: {sorted(set(schema_operators) - set(tables['union']))}. Only in the "
        f"union: {sorted(set(tables['union']) - set(schema_operators))}")
    for name, attrs in sorted(schema_operators.items()):
        expected = {sw for sw, table in TABLE_OF_SWITCH.items() if name in tables[table]}
        assert attrs == expected, (
            f"<{name}>: the schema offers {sorted(attrs)} and the tables back "
            f"{sorted(expected)}. A switch with no table behind it reads its default and "
            f"stays off; a table with no switch cannot be asked for at all")


@pytest.mark.wannierlib
def test_every_export_consumer_in_the_schema_is_read(tables):
    """A consumer the schema accepts but the reader never asks about is accepted and ignored."""
    text = SCHEMA.read_text()
    m = re.search(r'<xsd:complexType name="WannierlibExportType">(.*?)</xsd:complexType>',
                  text, re.S)
    assert m, "WannierlibExportType is gone from the schema"
    declared = set(re.findall(r'<xsd:element[^>]*\bname="(\w+)"', m.group(1)))
    read = set(re.findall(r"export_on\(xml,\s*xPathA,\s*'(\w+)'\)", READER.read_text()))
    assert declared == read, (
        f"accepted by the schema and never read: {sorted(declared - read)}; "
        f"read and not in the schema: {sorted(read - declared)}")
