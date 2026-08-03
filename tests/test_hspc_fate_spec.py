import json
import re
from pathlib import Path


ROOT = Path(__file__).parents[1]
SPEC = ROOT / "research" / "hspc-fate" / "v1-spec.md"
SOURCES = ROOT / "research" / "hspc-fate" / "sources.json"

REQ = re.compile(r"^- \*\*\[REQ-(\d{3})\]\*\* (.+)$")
AC = re.compile(r"^- \*\*\[AC-(\d{3})\]\*\* \(References REQ-(\d{3})\)$")
EARS = (
    re.compile(r"^The HSPC fate workflow shall .+\.$"),
    re.compile(r"^WHEN .+, the HSPC fate workflow shall .+\.$"),
    re.compile(r"^WHILE .+, the HSPC fate workflow shall .+\.$"),
    re.compile(r"^IF .+, THEN the HSPC fate workflow shall .+\.$"),
    re.compile(r"^WHERE .+, the HSPC fate workflow shall .+\.$"),
)


def test_each_requirement_has_one_ears_pattern_and_one_acceptance_block():
    lines = SPEC.read_text().splitlines()
    requirement_matches = [match for line in lines if (match := REQ.match(line))]
    criteria_matches = [match for line in lines if (match := AC.match(line))]
    requirements = {match[1]: match[2] for match in requirement_matches}
    criteria = {match[1]: match[2] for match in criteria_matches}

    expected = [f"{number:03d}" for number in range(1, len(requirements) + 1)]
    assert len(requirement_matches) == len(requirements)
    assert len(criteria_matches) == len(criteria)
    assert list(requirements) == expected
    assert list(criteria) == expected
    assert set(criteria.values()) == set(requirements)

    for requirement in requirements.values():
        assert sum(bool(pattern.fullmatch(requirement)) for pattern in EARS) == 1

    requirement_lines = {index for index, line in enumerate(lines) if REQ.match(line)}
    assert all(index in requirement_lines for index, line in enumerate(lines) if " shall " in line)


def test_each_acceptance_block_has_exactly_given_when_then():
    lines = SPEC.read_text().splitlines()
    acceptance_indexes = [index for index, line in enumerate(lines) if AC.match(line)]

    for keyword in ("Given", "When", "Then"):
        assert sum(line.startswith(f"  - **{keyword}** ") for line in lines) == len(acceptance_indexes)

    for index in acceptance_indexes:
        assert lines[index + 1] == "  - **Given** " + lines[index + 1].removeprefix("  - **Given** ")
        assert lines[index + 2] == "  - **When** " + lines[index + 2].removeprefix("  - **When** ")
        assert lines[index + 3] == "  - **Then** " + lines[index + 3].removeprefix("  - **Then** ")
        assert lines[index + 1].removeprefix("  - **Given** ")
        assert lines[index + 2].removeprefix("  - **When** ")
        assert lines[index + 3].removeprefix("  - **Then** ")
        if index + 4 < len(lines):
            assert not lines[index + 4].startswith("  - **")


def test_source_ledger_keeps_required_open_inputs_and_exclusions_explicit():
    ledger = json.loads(SOURCES.read_text())
    sources = {source["accession"]: source for source in ledger["sources"]}

    assert ledger["module"] == "hspc-fate"
    assert ledger["schema_version"] == "2.0.0"
    assert ledger["policy"]["data_in_git"] is False
    assert len(sources) == len(ledger["sources"])
    assert len({source["id"] for source in ledger["sources"]}) == len(ledger["sources"])
    assert sources["S-BSST1258"]["admission"] == "required_discovery"
    assert sources["GSE133181"]["admission"] == "required_validation"
    assert sources["S-BSST1258"]["access"] == "open"
    assert sources["GSE133181"]["access"] == "open"

    for source in ledger["sources"]:
        assert source["landing_url"].startswith("https://")
        assert source["asset_urls"]
        assert source["limitations"]
        if source["admission"] == "excluded":
            assert source["eligible_for_active_release"] is False
            assert source["exclusion_reason"]
