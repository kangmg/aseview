"""Validate JavaScript files have no syntax errors.

Uses Node.js (if available) to parse JS files, falling back to basic
bracket/brace balancing check.
"""

import os
import subprocess
import pytest

JS_DIR = os.path.join(os.path.dirname(__file__), "..", "aseview", "static", "js")
TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "..", "aseview", "templates")


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
    @pytest.fixture(
        params=[
            "molecular_viewer.html",
            "normal_viewer.html",
            "overlay_viewer.html",
            "frag_selector.html",
        ]
    )
    def template_file(self, request):
        path = os.path.join(TEMPLATE_DIR, request.param)
        if not os.path.exists(path):
            pytest.skip(f"{request.param} not found")
        return path

    def test_script_blocks_balanced(self, template_file):
        """All inline <script> blocks should have balanced brackets."""
        blocks = _extract_script_blocks(template_file)
        for i, block in enumerate(blocks):
            _check_bracket_balance(
                block, f"{template_file} <script> block #{i}"
            )

    def test_performance_functions_present(self, template_file):
        """Templates should reference performance utility functions."""
        with open(template_file, "r", encoding="utf-8") as f:
            content = f.read()
        # styles.js (which defines these functions) is inlined or linked
        # The template should reference at least the spatial grid function
        if "renderMolecule" in content:
            assert (
                "findBondsSpatialGrid" in content
                or "bondPairs" in content
            ), f"{template_file} should use spatial grid bond detection"
