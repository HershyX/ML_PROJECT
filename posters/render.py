"""Render the posters to A3 PDFs and PNG previews with hidden Chrome.

Run: python posters/render.py [poster file names...]
"""
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
names = sys.argv[1:] or ["poster1_project.html", "poster2_frontend.html"]

with sync_playwright() as p:
    browser = p.chromium.launch(channel="chrome", headless=True)
    page = browser.new_page(viewport={"width": 1123, "height": 1588}, device_scale_factor=2)
    for name in names:
        src = HERE / name
        page.goto(src.as_uri(), wait_until="networkidle")
        page.evaluate("document.fonts.ready")
        page.wait_for_timeout(500)
        overflow = page.evaluate("""() => [...document.querySelectorAll('.card .inner, .poster')]
            .filter(el => el.scrollHeight > el.clientHeight + 2).map(el => el.className)""")
        page.screenshot(path=HERE / f"{src.stem}.png", full_page=False)
        page.pdf(path=HERE / f"{src.stem}.pdf", width="297mm", height="420mm", print_background=True,
                 prefer_css_page_size=True)
        print(f"{src.stem}: pdf and png written, overflowing: {overflow or 'none'}")
    browser.close()
