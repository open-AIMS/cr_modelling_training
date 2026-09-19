# Provenance of the images on the published pages

Written 2026-09-17. The repository is public and the site is served from it, so
every image in `vignettes/images/` is republished to anyone who opens a module.
`CLAUDE.md` section 5 records that several are screenshots of publisher-typeset
journal pages and that each has to be checked before the course. That check is
done. Every source was looked up, and on 2026-09-20 the presenter decided that
every capture stays; the section *The licence of each source* records both. This
note is now the provenance record rather than a list of work outstanding.

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

## The licence of each source

Each source was looked up on 2026-09-19 in the Crossref record for its DOI,
which states the licence the publisher registered for the version of record. The *Journal of
Statistical Software* registers no licence with Crossref, so its article page
and its DOAJ record were read instead.

| Source | DOI | Licence of the version of record |
|---|---|---|
| @fisheretal2023, IEAM | `10.1002/ieam.4809` | CC BY 4.0 |
| @fisherfox2023, ET&C | `10.1002/etc.5610` | CC BY-NC 4.0 |
| Fisher et al. (2024), JSS | `10.18637/jss.v110.i05` | CC BY 4.0 in the article metadata; the site's general licence block says CC BY 3.0 |
| @Fox2010, Ecotoxicology and Environmental Safety | `10.1016/j.ecoenv.2009.09.012` | none; Elsevier registers a text-and-data-mining licence only |
| @Ritz2026, Environmental and Ecological Statistics | `10.1007/s10651-025-00698-y` | none; Springer registers a text-and-data-mining licence only |

Six files are settled by this and need only attribution, which each caption now
gives as the source and the licence together.

`Fisher_IEAM_Table1.png`, `NSEC_ieam.jpg` and `ieam_head.jpg` are from
@fisheretal2023 under CC BY 4.0, which permits reproduction for any purpose with
attribution. `NSEC_ieam.jpg` was shown in module 7 with no caption at all and now
has one.

`etnc_fig1.jpg` and `etnc_fig2.jpg` are from @fisherfox2023 under CC BY-NC 4.0.
The non-commercial condition turns on whether the use is directed toward
commercial advantage or monetary compensation. The published site is free to read
and shows no advertising, and the presenter recorded on 2026-09-19 that the
workshop is free to attend, so the site, the deck and anything handed out on the
day are all within the licence.

`Fisher_etal2024_JSS.png` is the title block of the JSS description of
`bayesnec`, under CC BY. The slide that shows it names the authors, the journal
and the issue in the sentence above it, which is the attribution.

Three files have no reuse licence. `Fox2010.png` and `Fox2010_nec.png` are pages
of a subscription Elsevier article, and `Ritz_etal2026.png` is a page of a
Springer article. For the two Fox captures the route that would have removed the
question is the one module 1 already takes in the paragraph below the figure: the
model is set out there as Quarto mathematics, so the capture of the typeset
equations adds the article's typography and nothing else.

### The decision on the unlicensed captures

Every capture stays. RF decided this on 2026-09-20, with the licence of each
source stated and with the offer to delete the two Fox captures, which module 1
does not depend on. Nothing is to be removed or redrawn on this account, and the
inventory above is kept as the record of what each file is and where it came
from rather than as a list of work outstanding.

## Captures from published articles, still in use

Each of these is named in `CLAUDE.md` section 5 as a publisher capture, or looks
like one. The third column is what the module's own caption says it is.

| File | Used in | What the caption says |
|---|---|---|
| `Fisher_IEAM_Table1.png` | module 3, line 551 | a table of toxicity estimates, "Reproduced from @fisheretal2023 under CC BY 4.0" |
| `NSEC_ieam.jpg` | module 4, line 113; module 7, line 979 | the N(S)EC figure, "Reproduced from @fisheretal2023 under CC BY 4.0" in both |
| `ieam_head.jpg` | module 4, line 824 | the title block of @fisheretal2023, "open access under CC BY 4.0" |
| `etnc_fig1.jpg` | module 3, line 113 | a threshold and a smooth curve, "Reproduced from @fisherfox2023 under CC BY-NC 4.0" |
| `etnc_fig2.jpg` | module 3, lines 47 and 529 | the four toxicity estimates, "Reproduced from @fisherfox2023 under CC BY-NC 4.0" |
| `necmod_fox2010.jpg` | module 3, line 171 | panel A redrawn from Fox (2010) |
| `Ritz_etal2026.png` | module 5, line 979 | Ritz, Gerhard and Streibig (2026), "the source of the measurements in this section" |
| `modelave_ecol.jpg` | module 4, line 64 | no caption; shown with `echo: false` |
| `glmbooks.jpg` | module 5, line 165 | no caption; shown with `echo: false` |
| `Fox2010.png` | module 1, line 195 | the title block of @Fox2010 |
| `Fox2010_nec.png` | module 1, line 197 | the definition of the model in @Fox2010 |
| `Fisher_etal2024_JSS.png` | the opening deck, slide *The presenters* | no caption; the slide text above it names the paper |

The last three rows were added on 2026-09-19 and were not in the mechanical
inventory. The licence of each source is in the section above.

`Ritz_etal2026.png` was captioned after the inventory was taken and the row
above now gives that caption. `modelave_ecol.jpg` and `glmbooks.jpg` still state
no source on the page they appear on, and neither was matched to an article by
the licence check, because the inventory records where a file is used and not
what is in it. That is a question of what the figure is rather than of whether it
may be reproduced, which the decision above settles.

The line numbers in the table were taken on 2026-09-17 and several are now
wrong. Find a file by name rather than by line.

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
uses any more. Its source is CC BY-NC 4.0, so holding it in the repository raises
nothing that needs undoing.

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
