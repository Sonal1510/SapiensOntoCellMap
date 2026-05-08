#!/usr/bin/env python3
"""
Tests: anatomically correct sex-specific organ display in the HTML visualizer.

Rules enforced:
  Female SVG  — breast shown (anat-sex-female), prostate hidden (anat-sex-male)
  Male SVG    — prostate shown (anat-sex-male), breast hidden (anat-sex-female)
  CSS per SVG — female SVG hides anat-sex-male; male SVG hides anat-sex-female

Run:
    python3 tests/test_anatomy_sex_correctness.py
    # or via pytest:
    pytest tests/test_anatomy_sex_correctness.py -v
"""

import os
import re
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TEMPLATE_PATH = os.path.join(PROJECT_ROOT, "src", "visualization", "visualizer_template.html")
OUTPUT_PATH   = os.path.join(PROJECT_ROOT, "data", "processed_combined_db", "sapiens_visualizer.html")


def _parse_svg_blocks(html: str) -> tuple[str, str]:
    """Return (female_block, male_block) by splitting on container div IDs."""
    try:
        fc = html.index('id="anatomy-container-female"')
        mc = html.index('id="anatomy-container-male"')
        me = html.index('</div><!-- /anatomy-container-male -->')
    except ValueError as e:
        raise AssertionError(f"Could not find anatomy container divs: {e}")
    return html[fc:mc], html[mc:me]


def _ellipse_sex_class(block: str, organ: str) -> list[str]:
    """Return list of CSS classes on ellipse elements for a given organ region."""
    hits = re.findall(
        r'<ellipse[^>]*region-' + organ + r'[^>]*/>',
        block
    )
    classes = []
    for h in hits:
        m = re.search(r'class="([^"]+)"', h)
        if m:
            classes.append(m.group(1))
    return classes


def _css_rule(block: str, selector: str) -> str | None:
    """Return the display value for a CSS class selector in a <style> block."""
    pattern = re.escape(selector) + r'\s*\{[^}]*display\s*:\s*(\w+)'
    m = re.search(pattern, block)
    return m.group(1) if m else None


def run_tests(html_path: str, label: str) -> list[str]:
    """Run all sex-correctness assertions. Returns list of failure messages."""
    failures = []

    if not os.path.exists(html_path):
        return [f"File not found: {html_path}"]

    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()

    try:
        female_block, male_block = _parse_svg_blocks(content)
    except AssertionError as e:
        return [str(e)]

    # ── CSS rules ────────────────────────────────────────────────────────────

    # Female container must hide anat-sex-male elements (scoped selector)
    css_f_hide_male = _css_rule(female_block, "#anatomy-container-female .anat-sex-male")
    if css_f_hide_male != "none":
        failures.append(
            f"[CSS] Female SVG: #anatomy-container-female .anat-sex-male should be display:none, got '{css_f_hide_male}'"
        )

    # Female container must show anat-sex-female elements (scoped selector)
    css_f_show_female = _css_rule(female_block, "#anatomy-container-female .anat-sex-female")
    if css_f_show_female != "inline":
        failures.append(
            f"[CSS] Female SVG: #anatomy-container-female .anat-sex-female should be display:inline, got '{css_f_show_female}'"
        )

    # Male container must hide anat-sex-female elements (scoped selector)
    css_m_hide_female = _css_rule(male_block, "#anatomy-container-male .anat-sex-female")
    if css_m_hide_female != "none":
        failures.append(
            f"[CSS] Male SVG: #anatomy-container-male .anat-sex-female should be display:none, got '{css_m_hide_female}'"
        )

    # Male container must show anat-sex-male elements (scoped selector)
    css_m_show_male = _css_rule(male_block, "#anatomy-container-male .anat-sex-male")
    if css_m_show_male != "inline":
        failures.append(
            f"[CSS] Male SVG: #anatomy-container-male .anat-sex-male should be display:inline, got '{css_m_show_male}'"
        )

    # ── Element classes ───────────────────────────────────────────────────────

    # Female SVG: breast must have anat-sex-female (shown in female)
    breast_f = _ellipse_sex_class(female_block, "breast")
    if not breast_f:
        failures.append("[ELEMENT] Female SVG: no breast ellipse found")
    elif not all("anat-sex-female" in c for c in breast_f):
        failures.append(
            f"[ELEMENT] Female SVG: breast ellipse missing anat-sex-female. Classes: {breast_f}"
        )

    # Female SVG: prostate must have anat-sex-male (hidden in female)
    prostate_f = _ellipse_sex_class(female_block, "prostate")
    if not prostate_f:
        failures.append("[ELEMENT] Female SVG: no prostate ellipse found")
    elif not all("anat-sex-male" in c for c in prostate_f):
        failures.append(
            f"[ELEMENT] Female SVG: prostate ellipse missing anat-sex-male (visible on female!). Classes: {prostate_f}"
        )

    # Male SVG: breast must have anat-sex-female (hidden in male)
    breast_m = _ellipse_sex_class(male_block, "breast")
    if not breast_m:
        failures.append("[ELEMENT] Male SVG: no breast ellipse found")
    elif not all("anat-sex-female" in c for c in breast_m):
        failures.append(
            f"[ELEMENT] Male SVG: breast ellipse missing anat-sex-female (visible on male!). Classes: {breast_m}"
        )

    # Male SVG: prostate must have anat-sex-male (shown in male)
    prostate_m = _ellipse_sex_class(male_block, "prostate")
    if not prostate_m:
        failures.append("[ELEMENT] Male SVG: no prostate ellipse found")
    elif not all("anat-sex-male" in c for c in prostate_m):
        failures.append(
            f"[ELEMENT] Male SVG: prostate ellipse missing anat-sex-male. Classes: {prostate_m}"
        )

    # ── Label text sex-gating ─────────────────────────────────────────────────

    # Female SVG: Breast label must have anat-sex-female class
    breast_label_f = re.findall(r'<text[^>]*anat-sex-([^"]+)"[^>]*>[^<]*[Bb]reast', female_block)
    if not breast_label_f or not all(s == "female" for s in breast_label_f):
        # Try alternate attribute order
        breast_label_f2 = re.findall(r'onclick="selectAnatomyRegion\(\'region-breast\'\)"[^>]*>Breast', female_block)
        breast_text_cls = re.findall(r'<text([^>]*region-breast[^>]*)>Breast', female_block)
        if not any("anat-sex-female" in t for t in breast_text_cls):
            failures.append("[LABEL] Female SVG: Breast label missing anat-sex-female class")

    # Male SVG: Prostate label must have anat-sex-male class
    prostate_labels_m = re.findall(r'<text([^>]*)>[Pp]rostate', male_block)
    if not any("anat-sex-male" in t for t in prostate_labels_m):
        failures.append("[LABEL] Male SVG: Prostate label missing anat-sex-male class")

    # Female SVG: Prostate label must NOT appear without sex gating
    prostate_labels_f = re.findall(r'<text([^>]*)>[Pp]rostate', female_block)
    if any("anat-sex-female" in t for t in prostate_labels_f):
        failures.append("[LABEL] Female SVG: Prostate label has anat-sex-female (would show on female!)")
    for t in prostate_labels_f:
        if "anat-sex-male" not in t:
            failures.append(f"[LABEL] Female SVG: Prostate label missing anat-sex-male gate: '{t.strip()[:80]}'")

    # Male SVG: Uterus label must NOT appear without sex gating
    uterus_labels_m = re.findall(r'<text([^>]*)>[Uu]terus', male_block)
    for t in uterus_labels_m:
        if "anat-sex-female" not in t:
            failures.append(f"[LABEL] Male SVG: Uterus label missing anat-sex-female gate: '{t.strip()[:80]}'")

    return failures


def main():
    all_passed = True

    for path, label in [(TEMPLATE_PATH, "TEMPLATE"), (OUTPUT_PATH, "GENERATED OUTPUT")]:
        print(f"\n{'='*60}")
        print(f"Testing: {label}")
        print(f"File:    {path}")
        print('='*60)
        failures = run_tests(path, label)
        if failures:
            all_passed = False
            print(f"FAILED — {len(failures)} issue(s):")
            for f in failures:
                print(f"  ✗ {f}")
        else:
            print("  ✓ Female SVG: breast shown (anat-sex-female), prostate hidden (anat-sex-male)")
            print("  ✓ Male SVG:   prostate shown (anat-sex-male), breast hidden (anat-sex-female)")
            print("  ✓ CSS rules:  each SVG correctly gates its sex-specific elements")
            print("  ✓ Labels:     sex-specific labels gated correctly")
            print(f"ALL PASS")

    if not all_passed:
        sys.exit(1)


# ── pytest-compatible test functions ─────────────────────────────────────────

def test_template_sex_correctness():
    failures = run_tests(TEMPLATE_PATH, "template")
    assert not failures, "\n".join(failures)

def test_generated_output_sex_correctness():
    if not os.path.exists(OUTPUT_PATH):
        import pytest
        pytest.skip("Generated output not found — run db_process.py first")
    failures = run_tests(OUTPUT_PATH, "generated output")
    assert not failures, "\n".join(failures)


if __name__ == "__main__":
    main()
