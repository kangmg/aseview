from typing import Any, Dict, Optional

from .export import ImageExportError


_EXPORT_SCRIPT = (
    "async ({ format, options }) => {"
    "const fn = format === 'gif' ? window.saveAsGIF : window.saveAsPNG;"
    "if (typeof fn !== 'function') {"
    "throw new Error(`Export function missing for ${format}`);"
    "}"
    "return await fn(options);"
    "}"
)


def render_export_data_url(
    html: str,
    *,
    export_format: str,
    options: Dict[str, Any],
    timeout_ms: int,
    browser_executable_path: Optional[str],
) -> Dict[str, Any]:
    try:
        from playwright.sync_api import Error as PlaywrightError
        from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
        from playwright.sync_api import sync_playwright
    except ModuleNotFoundError as exc:
        raise ImageExportError(_missing_playwright_message()) from exc

    try:
        result = _run_export(
            sync_playwright=sync_playwright,
            html=html,
            export_format=export_format,
            options=options,
            timeout_ms=timeout_ms,
            browser_executable_path=browser_executable_path,
        )
    except (PlaywrightError, PlaywrightTimeoutError) as exc:
        raise ImageExportError(str(exc)) from exc

    if not isinstance(result, dict):
        raise ImageExportError("Viewer export returned an invalid result")
    if result.get("ok") is False:
        message = result.get("message") or "Viewer export failed"
        raise ImageExportError(str(message))
    return dict(result)


def _run_export(
    *,
    sync_playwright: Any,
    html: str,
    export_format: str,
    options: Dict[str, Any],
    timeout_ms: int,
    browser_executable_path: Optional[str],
) -> Dict[str, Any]:
    viewport_width = int(options.get("width") or 1024)
    viewport_height = int(options.get("height") or 768)
    launch_options: Dict[str, Any] = {"headless": True}
    if browser_executable_path:
        launch_options["executable_path"] = browser_executable_path

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(**launch_options)
        try:
            page = browser.new_page(
                viewport={
                    "width": max(1, viewport_width),
                    "height": max(1, viewport_height),
                }
            )
            page.set_content(html, wait_until="networkidle", timeout=timeout_ms)
            page.wait_for_function(
                "() => typeof window.saveAsPNG === 'function' "
                "&& typeof window.saveAsGIF === 'function'",
                timeout=timeout_ms,
            )
            page.wait_for_function(
                "() => document.querySelector('canvas') !== null",
                timeout=timeout_ms,
            )
            return page.evaluate(
                _EXPORT_SCRIPT,
                {"format": export_format, "options": options},
            )
        finally:
            browser.close()


def _missing_playwright_message() -> str:
    return (
        "save_png/save_gif require the optional export dependency. "
        "Install with `pip install 'aseview[export]'`, then run "
        "`python -m playwright install chromium`."
    )
