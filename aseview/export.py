import base64
import math
from pathlib import Path
from typing import Any, Dict, Optional, Protocol


class HtmlViewer(Protocol):
    def get_html(self) -> str:
        ...


class ImageExportError(RuntimeError):
    pass


_PNG_QUALITY_TO_SCALE = {
    "low": 1.0,
    "medium": 2.0,
    "high": 3.0,
}

_GIF_QUALITY_TO_SAMPLE_INTERVAL = {
    "low": 20,
    "medium": 10,
    "high": 5,
}


def save_png(
    viewer: HtmlViewer,
    filename: str,
    *,
    width: Optional[int] = None,
    height: Optional[int] = None,
    scale: Optional[float] = None,
    transparent: bool = True,
    background_color: Optional[str] = None,
    quality: Optional[str] = None,
    timeout_ms: int = 30000,
    browser_executable_path: Optional[str] = None,
) -> Dict[str, Any]:
    output_path = Path(filename)
    options = _base_options(
        output_path,
        width=width,
        height=height,
        scale=_resolve_png_scale(scale, quality),
        transparent=transparent,
        background_color=background_color,
    )
    return _save_image_export(
        viewer,
        output_path,
        export_format="png",
        options=options,
        timeout_ms=timeout_ms,
        browser_executable_path=browser_executable_path,
    )


def save_gif(
    viewer: HtmlViewer,
    filename: str,
    *,
    width: Optional[int] = None,
    height: Optional[int] = None,
    scale: Optional[float] = None,
    transparent: bool = True,
    background_color: Optional[str] = None,
    frames: Optional[int] = None,
    fps: Optional[float] = None,
    delay: Optional[float] = None,
    sample_interval: Optional[int] = None,
    quality: Optional[str] = None,
    timeout_ms: int = 30000,
    browser_executable_path: Optional[str] = None,
) -> Dict[str, Any]:
    output_path = Path(filename)
    options = _base_options(
        output_path,
        width=width,
        height=height,
        scale=scale if scale is not None else 1.0,
        transparent=transparent,
        background_color=background_color,
    )
    if frames is not None:
        options["frames"] = _positive_int(frames, "frames")
    if fps is not None:
        options["fps"] = _positive_float(fps, "fps")
    if delay is not None:
        options["delay"] = _positive_float(delay, "delay")
    resolved_sample_interval = _resolve_gif_sample_interval(
        sample_interval,
        quality,
    )
    if resolved_sample_interval is not None:
        options["sampleInterval"] = resolved_sample_interval

    return _save_image_export(
        viewer,
        output_path,
        export_format="gif",
        options=options,
        timeout_ms=timeout_ms,
        browser_executable_path=browser_executable_path,
    )


def _base_options(
    output_path: Path,
    *,
    width: Optional[int],
    height: Optional[int],
    scale: float,
    transparent: bool,
    background_color: Optional[str],
) -> Dict[str, Any]:
    options: Dict[str, Any] = {
        "download": False,
        "returnDataUrl": True,
        "filename": output_path.name,
        "scale": _positive_float(scale, "scale"),
        "transparent": bool(transparent),
    }
    if width is not None:
        options["width"] = _positive_int(width, "width")
    if height is not None:
        options["height"] = _positive_int(height, "height")
    if background_color is not None:
        options["backgroundColor"] = background_color
    return options


def _resolve_png_scale(
    scale: Optional[float],
    quality: Optional[str],
) -> float:
    if scale is not None:
        return _positive_float(scale, "scale")
    if quality is None:
        return 1.0
    key = quality.lower()
    if key not in _PNG_QUALITY_TO_SCALE:
        raise ValueError("PNG quality must be one of: low, medium, high")
    return _PNG_QUALITY_TO_SCALE[key]


def _resolve_gif_sample_interval(
    sample_interval: Optional[int],
    quality: Optional[str],
) -> Optional[int]:
    if sample_interval is not None:
        return _positive_int(sample_interval, "sample_interval")
    if quality is None:
        return None
    key = quality.lower()
    if key not in _GIF_QUALITY_TO_SAMPLE_INTERVAL:
        raise ValueError("GIF quality must be one of: low, medium, high")
    return _GIF_QUALITY_TO_SAMPLE_INTERVAL[key]


def _positive_int(value: int, name: str) -> int:
    integer = int(value)
    if integer <= 0 or integer != value:
        raise ValueError(f"{name} must be a positive integer")
    return integer


def _positive_float(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be a positive number")
    return number


def _save_image_export(
    viewer: HtmlViewer,
    output_path: Path,
    *,
    export_format: str,
    options: Dict[str, Any],
    timeout_ms: int,
    browser_executable_path: Optional[str],
) -> Dict[str, Any]:
    if timeout_ms <= 0:
        raise ValueError("timeout_ms must be a positive integer")

    result = _render_export_data_url(
        viewer.get_html(),
        export_format=export_format,
        options=options,
        timeout_ms=timeout_ms,
        browser_executable_path=browser_executable_path,
    )
    data_url = result.pop("dataUrl", None)
    if not isinstance(data_url, str):
        raise ImageExportError("Viewer export did not return image data")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    image_bytes = _decode_image_data_url(data_url, export_format)
    output_path.write_bytes(image_bytes)

    result["path"] = str(output_path)
    result["bytes"] = len(image_bytes)
    return result


def _render_export_data_url(
    html: str,
    *,
    export_format: str,
    options: Dict[str, Any],
    timeout_ms: int,
    browser_executable_path: Optional[str],
) -> Dict[str, Any]:
    from ._export_playwright import render_export_data_url

    return render_export_data_url(
        html,
        export_format=export_format,
        options=options,
        timeout_ms=timeout_ms,
        browser_executable_path=browser_executable_path,
    )


def _decode_image_data_url(data_url: str, export_format: str) -> bytes:
    expected_prefix = f"data:image/{export_format};base64,"
    if not data_url.startswith(expected_prefix):
        raise ImageExportError(
            f"Viewer export returned invalid {export_format.upper()} data"
        )
    encoded = data_url[len(expected_prefix):]
    return base64.b64decode(encoded)
