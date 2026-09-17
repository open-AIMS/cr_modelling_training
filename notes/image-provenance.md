# Provenance of the images on the published pages

Written 2026-09-17. The repository is public and the site is served from it, so
every image in `vignettes/images/` is republished to anyone who opens a module.
`CLAUDE.md` section 5 records that several are screenshots of publisher-typeset
journal pages and that each has to be checked before the course. This note is
the list to check against, and it records the decision beside each file once one
is made.

The inventory was taken mechanically on 2026-09-17 by matching each filename in
`vignettes/images/` against the text of every `.qmd`. It says where a file is
used, not what is in it.

## The decision for each file

Three routes are available, in the order `CLAUDE.md` section 5 prefers them.

Redraw the figure from the quantities behind it. This is the only route that
removes the question entirely, and for a figure of a fitted curve it is usually
a short piece of code against data the course already ships.

Cite it and link to it rather than reproducing it. A reader following a link to
the publisher's page sees the figure under the publisher's own terms.

Reproduce it, where the article's licence permits reuse with attribution. An
open-access article under a Creative Commons licence usually does; a typeset
page from a subscription journal usually does not, and authorship of the article
does not by itself settle it, because copyright in the typeset version commonly
sits with the publisher. The licence has to be read for the specific article.

## Captures from published articles, still in use

Each of these is named in `CLAUDE.md` section 5 as a publisher capture, or looks
like one. The third column is what the module's own caption says it is.

| File | Used in | What the caption says |
|---|---|---|
| `Fisher_IEAM_Table1.png` | module 3, line 551 | a table of toxicity estimates, "Reproduced from @fisheretal2023" |
| `NSEC_ieam.jpg` | module 4, line 113; module 7, line 979 | the N(S)EC figure, "Reproduced from @fisheretal2023" |
| `ieam_head.jpg` | module 4, line 824 | the title block of @fisheretal2023, described there as open access |
| `etnc_fig1.jpg` | module 3, line 113 | a threshold and a smooth curve, "Reproduced from @fisherfox2023" |
| `etnc_fig2.jpg` | module 3, lines 47 and 529 | the four toxicity estimates, "Reproduced from @fisherfox2023" |
| `necmod_fox2010.jpg` | module 3, line 171 | panel A redrawn from Fox (2010) |
| `Ritz_etal2026.png` | module 5, line 57 | no caption; shown with `echo: false` |
| `modelave_ecol.jpg` | module 4, line 64 | no caption; shown with `echo: false` |
| `glmbooks.jpg` | module 5, line 165 | no caption; shown with `echo: false` |

Two things follow from the table itself, before any licence is read.

The four images with no caption state no source on the page. Whatever is decided
about reproducing them, a figure taken from somewhere has to name where, so each
needs a caption or has to go.

`@fisheretal2023` is open access, and module 4 says so in its own caption and
links the code repository for the paper. Read its licence first: if it is a
Creative Commons licence permitting reuse with attribution, three of the nine
files above are settled at once.

## Images no longer used

Eighteen files in `vignettes/images/` are referenced by no module as of
2026-09-17: `ametryn.JPG`, `bayesnec_landing.png`, `bayesnec_logo.jpg`,
`bnec_function.jpg`, `brm_args.jpg`, `dev_branch.jpg`, `drc_drmhelp.JPG`,
`drc_edcomp.JPG`, `drc_edhelp.JPG`, `drc_models.JPG`, `ecx_function.jpg`,
`etnc_header.jpg`, `nec_function.jpg`, `new_issue.jpg`, `nsec_function.jpg`,
`plot_args.jpg`, `rdata_object.jpg` and `website.jpg`. Most are screenshots of R
help pages and package websites taken for the 2023 learnr modules.

`etnc_header.jpg` is in this group. It is the ET&C title block of
@fisherfox2023, and it is the file `CLAUDE.md` section 5 confirms was checked on
2026-09-11, so the one capture whose identity is established is the one no page
uses any more.

Deleting a file here removes it from the working tree and not from the history,
so a capture that should never have been committed needs more than a delete.
Nothing in this group is worth that on its own account; the point of the list is
that the live question is nine files out of the fifty-four in the folder.

## Photographs and screenshots taken by the presenter

The remaining images in use are experimental photographs, Positron screenshots
and figures generated from the `bayesnec` vignettes. They raise no question of
this kind, with one qualification worth checking once: a photograph of a person,
a vessel or a field site taken during funded work may have its own attribution
requirement, and `experimental_setup.JPG`, `tank_offset.JPG`, `microtox.JPG`,
`dredging_effects.JPG` and `herbicides_climate.JPG` are the candidates.
