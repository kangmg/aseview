"""Validate JavaScript files have no syntax errors.

Uses Node.js (if available) to parse JS files, falling back to basic
bracket/brace balancing check.
"""

import os
import subprocess
import pytest

JS_DIR = os.path.join(os.path.dirname(__file__), "..", "aseview", "static", "js")
THEMES_DIR = os.path.join(os.path.dirname(__file__), "..", "aseview", "themes")
# Legacy fallback directory kept for backward compat
TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "..", "aseview", "templates")

TEMPLATE_NAMES = [
    "molecular_viewer.html",
    "normal_viewer.html",
    "overlay_viewer.html",
    "frag_selector.html",
]

def _available_theme_files():
    """Return (theme, filename) pairs for all installed themes."""
    result = []
    if not os.path.isdir(THEMES_DIR):
        return result
    for theme in sorted(os.listdir(THEMES_DIR)):
        theme_path = os.path.join(THEMES_DIR, theme)
        if not os.path.isdir(theme_path):
            continue
        for name in TEMPLATE_NAMES:
            path = os.path.join(theme_path, name)
            if os.path.exists(path):
                result.append((theme, name, path))
    return result


def _molecular_template_files():
    """Return all molecular viewer templates that should stay behaviorally aligned."""
    result = []
    base_path = os.path.join(TEMPLATE_DIR, "molecular_viewer.html")
    if os.path.exists(base_path):
        result.append(("templates", "molecular_viewer.html", base_path))
    result.extend(
        (theme, name, path)
        for theme, name, path in _available_theme_files()
        if name == "molecular_viewer.html"
    )
    return result


def _has_node():
    try:
        subprocess.run(["node", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def _extract_script_blocks(html_path):
    """Extract inline <script>...</script> blocks from an HTML file."""
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    blocks = []
    start = 0
    while True:
        idx = content.find("<script>", start)
        if idx == -1:
            break
        idx += len("<script>")
        end = content.find("</script>", idx)
        if end == -1:
            break
        blocks.append(content[idx:end])
        start = end + len("</script>")
    return blocks


def _check_bracket_balance(code, filepath=""):
    """Basic bracket/brace/paren balance check.

    Note: This is a simplified checker that skips regex literals and
    single-line comments. For authoritative syntax checking, use
    ``node --check`` instead.
    """
    stack = []
    pairs = {")": "(", "]": "[", "}": "{"}
    in_string = None
    escape_next = False
    i = 0
    n = len(code)

    while i < n:
        ch = code[i]

        if escape_next:
            escape_next = False
            i += 1
            continue

        if ch == "\\":
            escape_next = True
            i += 1
            continue

        # Handle string literals
        if in_string:
            if ch == in_string:
                in_string = None
            i += 1
            continue

        if ch in ('"', "'", "`"):
            in_string = ch
            i += 1
            continue

        # Skip single-line comments
        if ch == "/" and i + 1 < n and code[i + 1] == "/":
            end = code.find("\n", i)
            i = end + 1 if end != -1 else n
            continue

        # Skip block comments
        if ch == "/" and i + 1 < n and code[i + 1] == "*":
            end = code.find("*/", i + 2)
            i = end + 2 if end != -1 else n
            continue

        # Skip regex literals (heuristic: / after certain tokens)
        if ch == "/":
            # Look back to see if this could be a regex
            j = i - 1
            while j >= 0 and code[j] in " \t":
                j -= 1
            if j >= 0 and code[j] in "=(:,;[!&|?{}\n":
                # Likely a regex literal — skip to closing /
                k = i + 1
                while k < n and code[k] != "/" and code[k] != "\n":
                    if code[k] == "\\":
                        k += 1  # skip escaped char in regex
                    k += 1
                if k < n and code[k] == "/":
                    # Skip flags after /
                    k += 1
                    while k < n and code[k].isalpha():
                        k += 1
                    i = k
                    continue

        if ch in "([{":
            stack.append(ch)
        elif ch in ")]}":
            if not stack or stack[-1] != pairs[ch]:
                pytest.fail(
                    f"Unmatched '{ch}' at position {i} in {filepath}"
                )
            stack.pop()

        i += 1

    if stack:
        pytest.fail(f"Unclosed brackets {stack} in {filepath}")


# ── Tests ────────────────────────────────────────────────────


class TestJSFiles:
    @pytest.fixture(params=["styles.js", "aseview.js"])
    def js_file(self, request):
        path = os.path.join(JS_DIR, request.param)
        if not os.path.exists(path):
            pytest.skip(f"{request.param} not found")
        return path

    @pytest.mark.skipif(not _has_node(), reason="Node.js not available")
    def test_node_parse(self, js_file):
        """Parse JS file with Node.js to detect syntax errors."""
        result = subprocess.run(
            ["node", "--check", js_file],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"JS syntax error in {js_file}:\n{result.stderr}"
        )

    def test_bracket_balance(self, js_file):
        with open(js_file, "r", encoding="utf-8") as f:
            code = f.read()
        _check_bracket_balance(code, js_file)


class TestHTMLTemplateScripts:
    """Test all HTML templates across every installed theme."""

    @pytest.fixture(params=_available_theme_files(), ids=lambda t: f"{t[0]}/{t[1]}")
    def theme_template(self, request):
        return request.param  # (theme, name, path)

    def test_script_blocks_balanced(self, theme_template):
        """All inline <script> blocks should have balanced brackets."""
        _theme, _name, path = theme_template
        blocks = _extract_script_blocks(path)
        for i, block in enumerate(blocks):
            _check_bracket_balance(block, f"{path} <script> block #{i}")

    def test_performance_functions_present(self, theme_template):
        """Templates should reference performance utility functions."""
        _theme, _name, path = theme_template
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
        if "renderMolecule" in content:
            assert (
                "findBondsSpatialGrid" in content or "bondPairs" in content
            ), f"{path} should use spatial grid bond detection"


class TestMolecularConstraintHighlightWiring:
    """Keep fixed-constraint highlights wired through cinematic atom/bond styles."""

    @pytest.fixture(
        params=_molecular_template_files(),
        ids=lambda t: f"{t[0]}/{t[1]}",
    )
    def molecular_template(self, request):
        return request.param

    def test_fixed_constraints_pass_highlight_helpers(self, molecular_template):
        _theme, _name, path = molecular_template
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()

        required_snippets = [
            "const constraintHighlightColor = 0xffeb3b;",
            "highlightAtoms: fixedIndexSet",
            "highlightColor: constraintHighlightColor",
            "bondAtoms: [i, j]",
            "atomIndex: i",
            "styleName,\n                            bondHelpers",
            "styleName, atomHelpers",
        ]
        for snippet in required_snippets:
            assert snippet in content, f"{path} missing {snippet!r}"
