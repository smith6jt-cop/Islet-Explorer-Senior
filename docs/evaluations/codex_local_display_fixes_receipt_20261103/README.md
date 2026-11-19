# Latest CODEX display fixes: local source receipt

Received and read app commit `6d57be1a7fc5b0bff17a0caf68d1b7842d762f05` and font-build commit `7ac421ed7d4b5ab5425d0130e8ff7445f4e87240`, against `68c1cf6`. The [receipt](receipt.json) binds the five changed files. **The scope is local review only.** No port, implementation changes, tests, build, runtime action or deployment approval was performed here.

The Spatial change removes container-coordinate guide settings, uses default paper coordinates, measures shared guide and x-axis-title clearance, and positions guides after accounting for the shared margin. The tissue right margin increases from 18 to 24 points. These changes address the sender's reported behavior in pinned Plotly.js 2.11; runtime compatibility was not independently tested in this receipt.

The build change copies four DejaVu TTFs and their license from the existing version-pinned Matplotlib installation and registers them with fontconfig. It adds no package installation or dependency change. The installer records copied-font hashes rather than comparing against hardcoded individual font digests. Its checks cover four faces, distinct superscript signs, twelve distinct profile bitmaps in ragg/Cairo PNG, and embedded font/Unicode mappings in four Cairo PDFs. The PDF fixtures do not themselves test twelve distinct profile renderings.

Detailed [app source read](../../../outputs/codex_local_display_fixes_receipt_20261103/app_source_read.json) and [font source read](../../../outputs/codex_local_display_fixes_receipt_20261103/font_source_read.json) distinguish inspected code from execution evidence. The reported 601 R assertions, 35 parsed files, nine UMAP cases, four tissue modes and image-renderer results remain sender-reported. The new image's final evidence and runbook are forthcoming; receipt creation does not wait for or certify those results.

Earlier scientific outputs, receipts and finding-closure records remain unchanged. This new source pin is a separate candidate; earlier image or source validation does not automatically validate it. The sender reports public services unchanged, and this recipient performed no public-service action.
