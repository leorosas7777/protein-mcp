# Vendored assets

## 3Dmol-min.js

- **What:** the minified browser build of 3Dmol.js, embedded verbatim in every
  `console_<accession>.html` that `build_console_page` writes, so the page needs no network when
  it is opened.
- **Version:** 2.5.5.
- **Where it came from:** `https://unpkg.com/3dmol@2.5.5/build/3Dmol-min.js`.
- **SHA-256:** `f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016` (537 792 bytes).
  Verified on 2026-11-19 by downloading that URL again and comparing: byte for byte identical.
  The builds of 2.4.2, 2.5.2, 2.5.3 and 2.5.4 and the one served from `3dmol.csb.pitt.edu/build/`
  all differ, so the version is identified by the hash and not by a file name.
- **Licence:** BSD-3-Clause (the `license` field of the npm package metadata of 3dmol 2.5.5, and
  the text of `3Dmol-LICENSE.txt` below). It requires the copyright notice and the licence text
  to travel with any redistribution in binary form, and the page redistributes the bundle:
  `3Dmol-LICENSE.txt` is the licence of the package (`https://unpkg.com/3dmol@2.5.5/LICENSE`) and
  `3Dmol-min.js.LICENSE.txt` is the header file the bundle itself points to
  (`/*! For license information please see 3Dmol-min.js.LICENSE.txt */`). Both were downloaded on
  2026-11-19 from the same package version as the bundle.
- **How to update it:** download the new build and its two licence files from the same package
  version, record the new version and SHA-256 here, and regenerate the pages.
